#!/usr/bin/env python
"""Create workflow report."""

import argparse
from datetime import datetime
import json
import sys

from aplanat import report
from aplanat.components import fastcat
import aplanat.graphics
from bokeh.layouts import layout
from bokeh.models import BasicTickFormatter, Range1d
import pandas as pd


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser(sys.argv[1:])
    parser.add_argument(
        "output",
        help="Report output file.")
    parser.add_argument(
        "sample_name")
    parser.add_argument(
        "--reads_summary",
        required=True)
    parser.add_argument(
        "--eval_results",
        required=False)
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions", required=True,
        help="directory contained CSVs containing name,version.")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "wf-human-sv report", "wf-human-sv",
        revision=args.revision, commit=args.commit)

    #
    # Front matter
    #
    section = report_doc.add_section()
    section.markdown(f"```Sample: {args.sample_name}```")
    section.markdown(f"```Date: {datetime.today().strftime('%Y-%m-%d')}```")

    #
    # Evaluation results
    #
    section.markdown("## Evaluation results")
    if args.eval_results:
        data = None
        with open(args.eval_results) as f:
            data = json.load(f)
        section = report_doc.add_section()
        section.markdown("This sections displays the truvari"
                         " evaluation metrics for your SV calls.")
        exec_summary = aplanat.graphics.InfoGraphItems()
        exec_summary.append(
            "TP-call",
            "{:.2f}".format(data['TP-call']),
            "chart-pie"
        )
        exec_summary.append(
            "FP",
            "{:.2f}".format(data['FP']),
            "chart-pie"
        )
        exec_summary.append(
            "FN",
            "{:.2f}".format(data['FN']),
            "chart-pie"
        )
        exec_summary.append(
            "F1",
            "{:.2f}".format(data['f1']),
            "chart-pie"
        )
        exec_summary.append(
            "Precision",
            "{:.2f}".format(data['precision']),
            "chart-pie"
        )
        exec_summary.append(
            "Recall",
            "{:.2f}".format(data['recall']),
            "chart-pie"
        )
        exec_plot = aplanat.graphics.infographic(
            exec_summary.values(), ncols=4)
        section.plot(exec_plot, key="exec-plot")
    else:
        section.markdown(
            "This report was generated without evaluation"
            " results. To see them, re-run the workflow with"
            " --mode benchmark set.")

    #
    # Input dataset QC
    #
    reads_summary = args.reads_summary
    reads_summary_df = pd.read_csv(reads_summary, sep='\t')
    read_qual = fastcat.read_quality_plot(reads_summary_df)
    read_length = fastcat.read_length_plot(reads_summary_df)
    read_length.x_range = Range1d(0, 100000)
    read_length.xaxis.formatter = BasicTickFormatter(use_scientific=False)
    section = report_doc.add_section()
    section.markdown("## Read Quality Control")
    section.markdown("This sections displays basic QC"
                     " metrics indicating read data quality.")
    section.plot(
        layout(
            [[read_length, read_qual]],
            sizing_mode="stretch_width")
    )

    #
    # Params reporting
    #
    section = report_doc.add_section()
    section.markdown("## Workflow parameters")
    section.markdown("The table below highlights values of"
                     " the main parameters used in this analysis.")
    params = []
    with open(args.params) as f:
        params_data = json.load(f)
    for key, value in params_data.items():
        params.append((key, value))
    df_params = pd.DataFrame(params, columns=['Key', 'Value'])
    section.table(df_params, sortable=False, paging=False,
                  index=False, searchable=False)

    #
    # Software versions
    #
    section = report_doc.add_section()
    section.markdown("## Software versions")
    section.markdown('''The table below highlights versions
                    of key software used within the analysis''')
    versions = list()
    if args.versions is not None:
        with open(args.versions) as fh:
            for line in fh.readlines():
                name, version = line.strip().split(',')
                versions.append((name, version))
    versions = pd.DataFrame(versions, columns=('Name', 'Version'))
    section.table(versions, index=False)

    #
    # write report
    #
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
