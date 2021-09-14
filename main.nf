#!/usr/bin/env nextflow
import groovy.json.JsonBuilder

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
wf-human-sv

Usage:
    nextflow run epi2melabs/wf_human_sv [options]

Script Options:
    --help
    --fastq                   DIR         FASTQ file (required)
    --reference               FILE        FASTA format reference sequence (required)
    --sample                  STR         Name for sample being analysed (optional)
    --out_dir                 DIR         Path for output (default: $params.out_dir)
    --mode                    STR         Switch between standard or benchmark modes. (default: $params.mode)
    --min_read_support        INT/STR     Minimum read support required to call a SV (default: auto)
    --target_bedfile          FILE        BED file, SVs will only be called in these regions (optional)
    --report_name             STR     Optional report suffix (default: $params.report_name)

Benchmarking:
    In benchmarking mode the calls made using the pipeline will
    be evaluated against the supplied truthset. The report will
    be populated with the results. Paths to truthset files are
    set within the config file, and may be URLs or absolute paths.
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
        path "${reference.simpleName}", emit: ref
        path "${reference.simpleName}.gli", emit: lra_index
        path "${reference.simpleName}.mmi", emit: mmi_index
    script:
        def simpleRef = reference.simpleName
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
    script:
        def sample = params.sample ? params.sample : 'sample'
    """
    catfishq -r $reads --max_mbp $params.max_bp \
    | seqtk seq -A - \
    | lra align -ONT -t $task.cpus $reference - -p s \
    | samtools addreplacerg -r \"@RG\tID:$sample\tSM:$sample\" - \
    | samtools sort -@ $task.cpus -o lra.bam -
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
        def sample = params.sample ? params.sample : 'sample'
    """
    cuteSV \
        --threads $task.cpus \
        --sample $sample \
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
        file target_bed
    output:
        path "depth.regions.bed.gz", emit: mosdepth_bed
    script:
        def target_bed = target_bed.name != 'OPTIONAL_FILE' ? "${target_bed}" : 1000000
    """
	mosdepth \
        -x \
        -t $task.cpus \
        -b $target_bed \
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
        path "${vcf.simpleName}_filtered.vcf", emit: vcf
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

    bash filter.sh > ${vcf.simpleName}_filtered.vcf
	"""
}


process sortVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "${vcf.simpleName}_sorted.vcf", emit: vcf
    """
    vcfsort $vcf > ${vcf.simpleName}_sorted.vcf
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


process fastcatQC {
    label "wf_human_sv"
    cpus 1
    input:
        file reads
    output:
        path "per_read_stats.txt", emit: stats
    """
    fastcat --read per_read_stats.txt $reads > /dev/null
    """
}


// The following processes are for the benchmarking
// pathway of this pipeline. Please see the documentation
// for further details. 
process getTruthset {
    label "wf_human_sv"
    cpus 1
    output:
        path "*.vcf.gz", emit: truthset_vcf_gz
        path "*.vcf.gz.tbi",  emit: truthset_vcf_tbi
        path "*.bed",  emit: truthset_bed
    """
    for item in $params.truthset_vcf $params.truthset_index $params.truthset_bed
    do
        if wget -q --method=HEAD \$item;
        then
            echo "Downloading \$item."
            wget \$item
        elif [ -f \$item ];
        then
            echo "Found \$item locally."
            cp \$item .
        else
            echo "\$item cannot be found."
            exit 1
        fi
    done
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

    if [ ! -s eval_highconf.bed ]
    then
        echo "No overlaps found between calls and truthset"
        echo "Chr names in your reference and truthset may differ"
        exit 1
    fi
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
    script:
        def sample = params.sample ? params.sample : 'sample'
    """
    TRUVARI=\$(which truvari)
    python \$TRUVARI bench \
        --passonly \
        --pctsim 0 \
        -b $truthset_vcf \
        -c $calls_vcf \
        -f $reference \
        -o $sample \
        --includebed $calls_highconf_bed
    mv $sample/summary.txt summary.json
    """
}


process report {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file bam
        file bam_index
        file read_stats
        file eval_json
    output:
        path "wf-human-sv-*.html", emit: html
        path "nanoplot.tar.gz", emit: nanoplot
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
        def paramsMap = params.toMapString()
        def sample = params.sample ? params.sample : 'sample'
        def report_name = "wf-human-sv-" + params.report_name + '.html'

    """
    # Explicitly get software versions
    TRUVARI=\$(which truvari)
    python \$TRUVARI version | sed 's/ /,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    cuteSV --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    echo `lra -v | head -n 2 | tail -1 | cut -d ':' -f 2 | sed 's/ /lra,/'` >> versions.txt
    echo `seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d ':' -f 2 | sed 's/ /seqtk,/'` >> versions.txt

    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json

    # Generate wf-human-sv aplanat static report
    report.py \
        $report_name \
        $sample \
        $vcf \
        --reads_summary $read_stats \
        --params params.json \
        --versions versions.txt \
        --revision $workflow.revision \
        --commit $workflow.commitId \
        $evalResults 

    # Generate supplemental nanoplot report
    mkdir -p nanoplot
    NanoPlot \
        -t $task.cpus \
        -p nanoplot/${sample}_ \
        --bam $bam \
        --raw \
        --N50 \
        --title $sample \
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
        mosdepth(mapLRA.out.bam, mapLRA.out.bam_index, target)
        filterCalls(cuteSV.out.vcf, mosdepth.out.mosdepth_bed, target)
        sortVCF(filterCalls.out.vcf)
        indexVCF(sortVCF.out.vcf)
        fastcatQC(reads)
    emit:
        ref = indexLRA.out.ref
        vcf = indexVCF.out.vcf_gz
        vcf_index = indexVCF.out.vcf_tbi
        bam = mapLRA.out.bam
        bam_index = mapLRA.out.bam_index
        read_stats = fastcatQC.out.stats
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
            standard.vcf,
            standard.bam, 
            standard.bam_index,
            standard.read_stats,
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
        truthset = getTruthset()
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
            standard.vcf,
            standard.bam,
            standard.bam_index,
            standard.read_stats,
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

