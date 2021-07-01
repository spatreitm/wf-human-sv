#!/usr/bin/env nextflow
import groovy.json.JsonBuilder

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
wf-human-sv

Usage:
    nextflow run epi2melabs/wf_human_sv [options]

Script Otions:
    --help
    --fastq                   DIR         FASTQ file (required)
    --reference               FILE        FASTA format reference sequence (required)
    --sample                  STR         Name for sample being analysed (default: $params.sample)
    --out_dir                 DIR         Path for output (default: $params.out_dir)
    --min_read_support        INT/STR     Minimum read support required to call a SV (default: auto)
    --target_bedfile          FILE        BED file, SVs will only be called in these regions (optional)
    --mode                    STR         Switch between standard or benchmark modes. (default: $params.mode)

Benchmarking:
    In benchmarking mode the calls made using the pipeline will
    be evaluated against the supplied truthset. The report will
    be populated with the results.
"""
}


// Workflow utilities
def getFastq(path) {
    reads = file("${path}{**,.}/*.{fastq,fastq.gz,fq,fq.gz}", type: 'file')
    // If we have found some files, return them
    if (reads) {
        println("--fastq: Detected folder with ${reads.size()} input file(s).")
        return reads
    }
    // If we haven't found a list of fastq, perhaps it is a file
    reads = file("${path}", type: 'file')
    if (reads.exists()) {
        println("--fastq: Single input file detected.")
        return reads
    }
    println("--fastq: File(s) not detected, check path.")
    return false
}


def getReference(path) {
    reference = file(params.reference, type: "file")
    if (!reference.exists()) {
        println("--reference: File not detected, check path.")
        return false
    }
    return reference
}


def checkMinReadSupport(min_read_support) {
    if (min_read_support.toString().isInteger()) {
        return min_read_support
    } else if (min_read_support !== 'auto') {
        println("--min_read_support: Must be integer or 'auto'.")
        return false
    } else {
        return min_read_support
    }
}


// Workflow processes
process indexLRA {
    label "wf_human_sv"
    cpus 1
    input:
        file reference
    output:
        path "${reference.simpleName}.fna", emit: ref
        path "${reference.simpleName}.fna.gli", emit: lra_index
        path "${reference.simpleName}.fna.mmi", emit: mmi_index
    script:
        def simpleRef = reference.simpleName + '.fna'
    """
    cp -L $reference $simpleRef
    if [[ $reference == *.gz ]]
    then
        rm $reference
        mv $simpleRef ${simpleRef}.gz
        gunzip ${simpleRef}.gz
    fi
    lra index -ONT $simpleRef
    """
}


process mapLRA {
    label "wf_human_sv"
    cpus params.threads
    input:
        file reference
        file lra_index
        file mmi_index
        file reads
    output:
        path "lra.bam", emit: bam
        path "lra.bam.bai", emit: bam_index
    """
    catfishq -r $reads --max_mbp $params.max_bp \
    | seqtk seq -A - \
    | lra align -ONT -t $task.cpus $reference - -p s \
    | samtools addreplacerg -r \"@RG\tID:$params.sample\tSM:$params.sample\" - \
    | samtools sort -@ $task.cpus -o lra.bam -m 2G -
    samtools index -@ $task.cpus lra.bam
    """
}


process cuteSV {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
        file reference
    output:
        path "cutesv.vcf", emit: vcf
    script:
        def simpleRef = reference.simpleName
    """
    cuteSV \
        --threads $task.cpus \
        --sample $params.sample \
        --retain_work_dir \
        --report_readid \
        --genotype \
        --min_size $params.min_sv_length \
        --max_size $params.max_sv_length \
        --min_support $params.min_read_support_limit \
        --max_cluster_bias_INS $params.max_cluster_bias_INS \
        --diff_ratio_merging_INS $params.diff_ratio_merging_INS \
        --max_cluster_bias_DEL $params.max_cluster_bias_DEL \
        --diff_ratio_merging_DEL $params.diff_ratio_merging_DEL \
        $bam \
        $reference \
        cutesv.vcf \
        .
	"""
}


process mosdepth {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
    output:
        path "depth.regions.bed.gz", emit: mosdepth_bed
    """
	mosdepth \
        -x \
        -t $task.cpus \
        -b 1000000 \
        depth \
        $bam
	"""
}


process filterCalls {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file mosdepth_bed
        file target_bed
    output:
        path "filtered.${vcf}", emit: vcf
    script:
        def sv_types_joined = params.sv_types.join(" ")
        def target_bed = target_bed.name != 'OPTIONAL_FILE' ? "--target_bedfile ${target_bed}" : ""
    """
    get_filter_calls_command.py \
        $target_bed \
        --vcf $vcf \
        --depth_bedfile $mosdepth_bed \
        --min_sv_length $params.min_sv_length \
        --max_sv_length $params.max_sv_length \
        --sv_types $sv_types_joined \
        --min_read_support $params.min_read_support \
        --min_read_support_limit $params.min_read_support_limit > filter.sh

    bash filter.sh > filtered.${vcf}
	"""
}


process sortVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "sorted.${vcf}", emit: vcf
    """
    vcfsort $vcf > sorted.${vcf}
    """
}


process indexVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "${vcf}.gz", emit: vcf_gz
        path "${vcf}.gz.tbi", emit: vcf_tbi
    """
    cat $vcf | bgziptabix ${vcf}.gz
    """
}


// The following processes are for the benchmarking
// pathway of this pipeline. Please see the documentation
// for further details. 
process downloadTruthset {
    label "wf_human_sv"
    cpus 1
    output:
        path "*.vcf.gz", emit: truthset_vcf_gz
        path "*.vcf.gz.tbi",  emit: truthset_vcf_tbi
        path "*.bed",  emit: truthset_bed
    """
    wget $params.truthset_vcf \
    && wget $params.truthset_index \
    && wget $params.truthset_bed
    """
}


process intersectCallsHighconf {
    label "wf_human_sv"
    cpus 1
    input:
        file calls_vcf
        file truthset_bed
    output:
        path "eval_highconf.bed", emit: calls_highconf_bed
    """
    bedtools intersect \
        -a $truthset_bed \
        -b $calls_vcf \
        -u > eval_highconf.bed
    """
}


process excludeNonIndels {
    label "wf_human_sv"
    cpus 1
    input:
        file calls_vcf
    output:
        path "indels_${calls_vcf}", emit: indels_only_vcf_gz
        path "indels_${calls_vcf}.tbi", emit: indels_only_vcf_tbi
    """
    zcat $calls_vcf \
    | sed 's/SVTYPE=DUP/SVTYPE=INS/g' \
    | bcftools view -i '(SVTYPE = \"INS\" || SVTYPE = \"DEL\")' \
    | bgziptabix indels_${calls_vcf}
    """
}


process truvari {
    label "wf_human_sv"
    cpus 1
    input:
        file reference
        file calls_vcf
        file calls_vcf_tbi
        file truthset_vcf
        file truthset_vcf_tbi
        file calls_highconf_bed
    output:
        path "summary.json", emit: truvari_json
    """
    truvari bench \
        --passonly \
        --pctsim 0 \
        -b $truthset_vcf \
        -c $calls_vcf \
        -f $reference \
        -o $params.sample \
        --includebed $calls_highconf_bed
    mv $params.sample/summary.txt summary.json
    """
}


process report {
    label "wf_human_sv"
    cpus 1
    input:
        file reads
        file lra_bam
        file lra_bam_index
        file eval_json
    output:
        path "wf-human-sv-report.html", emit: html
        path "nanoplot.tar.gz", emit: nanoplot
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
        def paramsMap = params.toMapString()
    """
    echo '$paramsJSON' > params.json
    fastcat --read samples_reads.txt $reads > /dev/null

    report.py \
        wf-human-sv-report.html \
        $params.sample \
        --reads_summary samples_reads.txt \
        --params params.json \
        --revision $workflow.revision \
        --commit $workflow.commitId \
        $evalResults 

    mkdir -p nanoplot
    NanoPlot \
        -t $task.cpus \
        -p nanoplot/${params.sample}_ \
        --bam $lra_bam \
        --raw \
        --N50 \
        --title $params.sample \
        --downsample 100000
    tar -czvf nanoplot.tar.gz nanoplot
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wf_human_sv"

    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: { 
        f -> params.sample ? "${params.sample}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// Workflow main pipeline
workflow pipeline {
    take:
        reference
        reads
        target
    main:
        indexLRA(reference)
        mapLRA(indexLRA.out.ref, indexLRA.out.lra_index, indexLRA.out.mmi_index, reads)
        cuteSV(mapLRA.out.bam, mapLRA.out.bam_index, indexLRA.out.ref)
        mosdepth(mapLRA.out.bam, mapLRA.out.bam_index)
        filterCalls(cuteSV.out.vcf, mosdepth.out.mosdepth_bed, target)
        sortVCF(filterCalls.out.vcf)
        indexVCF(sortVCF.out.vcf)
    emit:
        ref = indexLRA.out.ref
        vcf = indexVCF.out.vcf_gz
        vcf_index = indexVCF.out.vcf_tbi
        bam = mapLRA.out.bam
        bam_index = mapLRA.out.bam_index
}


workflow standard {
    take:
        reference
        reads
        target
        optional_file
    main:
        println("================================")
        println("Running workflow: standard mode.")
        standard = pipeline(reference, reads, target)
        report(
            reads, 
            standard.bam, 
            standard.bam_index, 
            optional_file)
        results = report.out.html.concat(
            report.out.nanoplot,
            standard.vcf, 
            standard.vcf_index, 
            standard.bam, 
            standard.bam_index)
    emit:
        results
}


workflow benchmark {
    take:
        reference
        reads
        target
    main:
        println("=================================")
        println("Running workflow: benchmark mode.")
        standard = pipeline(reference, reads, target)
        truthset = downloadTruthset()
        filtered = excludeNonIndels(standard.vcf)
        highconf = intersectCallsHighconf(
            standard.vcf, 
            truthset.truthset_bed)
        truvari(
            standard.ref,
            filtered.indels_only_vcf_gz,
            filtered.indels_only_vcf_tbi,
            truthset.truthset_vcf_gz,
            truthset.truthset_vcf_tbi,
            highconf.calls_highconf_bed)
        report(
            reads, 
            standard.bam, 
            standard.bam_index, 
            truvari.out.truvari_json)
        results = report.out.html.concat(
            report.out.nanoplot,
            standard.vcf,
            standard.vcf_index, 
            standard.bam, 
            standard.bam_index)
    emit:
        results
}


// workflow entrypoint
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    if (!params.reference) {
        helpMessage()
        println("")
        println("`--reference` is required")
        exit 1
    }

    // Ready the optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    // Checking user parameters
    println("=================================")
    println("Checking inputs")

    // Acquire reads and reference
    reads = getFastq(params.fastq)
    reference = getReference(params.reference)
    if (!reads || !reference) {
        exit 1
    }

    // Check for target bedfile
    target = file(params.target_bedfile)
    if (!target.exists()) {
        target = OPTIONAL
    }

    // Check for sample name
    if (!params.sample) {
        println("--sample: Not set, please supply a name")
        exit 1
    }

    // Check min_read_support
    if (!checkMinReadSupport(params.min_read_support)) {
        exit 1
    }

    // Print all params
    println("=================================")
    println("Summarising parameters")
    params.each { it -> println("> $it.key: $it.value") }

    // Execute workflow
    switch(params.mode) {
        case "standard": 
            results = standard(reference, reads, target, OPTIONAL)
            output(results)
            break;
        case "benchmark": 
            results = benchmark(reference, reads, target)
            output(results)
            break;
        default:
            results = standard(reference, reads, target)
            output(results)
            break;
    } 
}

