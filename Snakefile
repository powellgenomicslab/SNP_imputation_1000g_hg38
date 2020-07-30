#!/usr/local/envs/py36/bin python3

import os
import sys
import pandas as pd
from glob import glob
import datetime
from datetime import datetime


### Set the date and time of run starting
date_now = datetime.date(datetime.now())
time_now = datetime.time(datetime.now())
datetime_now = str(date_now) + "_" + str(time_now)

# Extract variables from configuration file for use within the rest of the pipeline
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["refs"]
popscle_dict = config["popscle"]
popscle_extra_dict = config["popscle_extra"]
souporcell_dict = config["souporcell"]
souporcell_extra_dict = config["souporcell_extra"]
DoubletDetection_dict = config["DoubletDetection"]
DoubletDetection_manual_dict = config["DoubletDetection_manual"]
DoubletDetection_extra_dict = config["DoubletDetection_extra"]
scrublet_dict = config["scrublet"]
scrublet_manual_dict = config["scrublet_manual"]
scrublet_extra_dict = config["scrublet_extra"]
scds_dict = config["scds"]
CombineResults_dict = config["CombineResults"]

sys.path.append(input_dict["pipeline_dir"] + '/mods') 
import prepareArguments

# Use prepareArguments.py script to retrieve exact directories of single cell files
scrnaseq_libs_df = prepareArguments.get_scrnaseq_dirs(config)
scrnaseq_libs_df.to_csv(os.path.join(output_dict["output_dir"],'file_directories.txt'), sep = "\t", index = False)


# Get list of pools to process
samples = pd.read_csv(input_dict["samplesheet_filepath"], sep = "\t")
samples.columns = ["Pool", "N"]


### If the scrublet_check output is present => all the contents are there that are needed to move past 
if os.path.exists(output_dict["output_dir"] + "/scrublet/scrublet_check.done"):
    scrublet_decisions = pd.read_csv(output_dict["output_dir"] + "/manual_selections/scrublet_gene_pctl.txt", sep = "\t")

# Includes
include: input_dict["pipeline_dir"] + "/includes/Snakefile_popscle.smk"
include: input_dict["pipeline_dir"] + "/includes/Snakefile_souporcell.smk"
include: input_dict["pipeline_dir"] + "/includes/Snakefile_scrublet.smk"
include: input_dict["pipeline_dir"] + "/includes/Snakefile_scds.smk"
include: input_dict["pipeline_dir"] + "/includes/Snakefile_DoubletDetection.smk"
include: input_dict["pipeline_dir"] + "/includes/Snakefile_CombineResults.smk"


demuxlet_files = []
demuxlet_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/demuxlet_results.txt",  pool=samples.Pool))

souporcell_files = []
souporcell_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/souporcell_results.txt", pool=samples.Pool))
souporcell_files.append(expand(output_dict["output_dir"] + "/{pool}/souporcell/Individual_genotypes_subset.vcf.gz", pool=samples.Pool))

scds_files = []
scds_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/scds_results.txt", pool=samples.Pool))

### the scrublet files that will be run are dependent on user inputs in the yaml file
scrublet_files = []
scrublet_files.append(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv")
if os.path.exists(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv"):
    scrublet_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv", sep = "\t")
    if scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection):
        scrublet_selection["scrublet_Percentile"] = scrublet_selection["scrublet_Percentile"].astype(int)
        scrublet_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/{pctl}_scrublet_results.txt", zip, pool=scrublet_selection.Pool, pctl = scrublet_selection.scrublet_Percentile))
    elif scrublet_selection["scrublet_Percentile"].count() != len(scrublet_selection):
        if scrublet_manual_dict["run_scrublet_manual"] == False:
            scrublet_files.append(expand(output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/default_run_variables.txt", pool = samples.Pool, pctl = scrublet_dict["percentile"]))
        elif scrublet_manual_dict["run_scrublet_manual"] == True:
            scrublet_files.append(expand(output_dict["output_dir"] + "/{pool}/scrublet_{pctl}/manual_rerun_variables_" + datetime_now + ".txt",zip, pool = scrublet_manual_dict["scrublet_manual_threshold_pools"], pctl = scrublet_manual_dict["scrublet_manual_threshold_percentiles"]))

DoubletDetection_files = []
DoubletDetection_files.append(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv")
if os.path.exists(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv"):
    DoubletDetection_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv", sep = "\t")
    if len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) == len(DoubletDetection_selection):
        DoubletDetection_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/DoubletDetection_results.txt", pool=DoubletDetection_selection.Pool))
    elif len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) != len(DoubletDetection_selection):
        if DoubletDetection_manual_dict["run_DoubletDetection_manual"] == False:
            DoubletDetection_files.append(expand(output_dict["output_dir"] + "/{pool}/DoubletDetection/default_run_variables.txt", pool = samples.Pool))
        elif DoubletDetection_manual_dict["run_DoubletDetection_manual"] == True:
            DoubletDetection_files.append(expand(output_dict["output_dir"] + "/{pool}/DoubletDetection/manual_rerun_variables_" + datetime_now + ".txt", pool = DoubletDetection_manual_dict["DoubletDetection_manual_pools"]))


combined_files = []
if os.path.exists(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv") and os.path.exists(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv"):
    scrublet_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/scrublet/scrublet_percentile_manual_selection.tsv", sep = "\t")
    DoubletDetection_selection = pd.read_csv(output_dict["output_dir"] + "/manual_selections/DoubletDetection/DoubletDetection_manual_selection.tsv", sep = "\t")
    if scrublet_selection["scrublet_Percentile"].count() == len(scrublet_selection) and len(DoubletDetection_selection[DoubletDetection_selection['DoubletDetection_PASS_FAIL'].astype(str).str.contains('PASS', na=False)]) == len(DoubletDetection_selection):
        combined_files.append(expand(output_dict["output_dir"] + "/{pool}/CombinedResults/Final_Assignments_demultiplexing_doublets.txt", pool = samples.Pool))
        combined_files.append(output_dict["output_dir"] + "/QC_figures/UMI_vs_Genes_QC_scatter.png")


rule all:
    input:
        demuxlet_files,
        souporcell_files,
        scds_files,
        scrublet_files,
        DoubletDetection_files,
        combined_files


