import os
import subprocess
import argparse
from time import sleep
from datetime import datetime

import json
import pandas as pd

# Downloads manifest file from gclound for lrwgs samples.
def download_manifest_file(manifest_file):
    print("Downloading the manifest file for lrwgs ...")
    subprocess.call("gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v7/wgs/long_read/manifest.csv .")
    subprocess.call("mv manifest.csv {}".format(manifest_file))

# Selects samples based on the counts and returs paths to bam files on gcloud.
def select_samples(sample_count, dataset):
    if dataset == "lrwgs":
        manifest_file = "manifest_lrwgs.csv"
        selected_columns = ["grch38-bam", "grch38-bai"]
        if not os.path.exists(manifest_file):
            print("Warning: The manifest file for long reads WGS is not present or not correctly named.")
            download_manifest_file(manifest_file)
    elif dataset == "srwgs":
        selected_columns = ["cram_uri", "cram_index_uri"]
        manifest_file = "manifest_srwgs.csv"


    files_df = pd.read_csv(manifest_file)
    num_individuals = len(files_df)
    print("Working with manifest file of {} columns and {} samples".format(
            len(files_df.columns), num_individuals))
    selected_samples = files_df.sample(n=sample_count, random_state=1)
    print("Selected {} random sample(s)".format(sample_count))

    selected_samples = selected_samples[selected_columns]
    selected_samples = selected_samples.rename(columns={
                            selected_columns[0] : "alignment_file",
                            selected_columns[1] : "alignment_index_file"})

    return selected_samples[["alignment_file", "alignment_index_file"]]

# Write a file indicating one target region per line.
def write_region_file(regions, filename):
    with open(filename, "w+") as region_file:
        for region in regions:
            region_file.write(region + "\n")

# Write a file indicating gcloud path to one bam file  per line.
def write_bam_list(bams, filename):
    with open(filename, "w+") as bam_list_file:
        for bam in bams:
            bam_list_file.write(bam + "\n")

# Write the input json file for when results from all batches are being merged.
def write_input_merge_json(input_json_filename,
                     batch_vcf_filenames,
                     batch_vcf_index_filenames,
                        ):
    data = {"run_merge_advntr.per_batch_vcfs": batch_vcf_filenames,
            "run_merge_advntr.per_batch_vcf_indexes": batch_vcf_index_filenames}

    with open(input_json_filename, "w+") as input_json_file:
        json.dump(data, input_json_file)

# Write input json file based on selected sample(s) bam files.
def write_input_json(input_json_filename,
                     samples_files,
                     region,
                     vntr_id,
                     google_project):
    try: # Running on AoU
        token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
            capture_output=True, check=True, encoding='utf-8')
        token = str.strip(token_fetch_command.stdout)
    except (TypeError, AttributeError, FileNotFoundError): # Running locally
        token = ""


    # TODO: Write a file indicating all bam files for genotyping. This will be passed to the WDL task.
    #bam_list = "./bam_list.txt"
    #write_bam_list(bams = [sample_df['grch38-bam']], filename=bam_list)

    # TODO: Write a file indicating all regions for genotyping. This will be passed to the WDL task.
    #region_file = "./regions.txt"
    #target_regions = ["chr15:88855424-88857434"]
    #write_region_file(regions=target_regions, filename=region_file)

    data = {"run_advntr.bam_files": samples_files,
            "run_advntr.region": region,
            "run_advntr.vntr_id": vntr_id,
            "run_advntr.gcloud_token": token,
            "run_advntr.google_project": google_project}

    with open(input_json_filename, "w+") as input_json_file:
        json.dump(data, input_json_file)

# Write options file including output directory in the bucket.
def write_options_json(options_json_filename, output_path_gcloud):
    data = {"jes_gcs_root": output_path_gcloud,
            "workflow_failure_mode": "ContinueWhilePossible",
            #"useDockerImageCache": true,
            #"default_runtime_attributes": {"bootDiskSizeGb": 3},
            }
    with open(options_json_filename, "w+") as options_json_file:
        json.dump(data, options_json_file)


def run_single_command(command):
    process = subprocess.run(command)
    if process.stdout is not None:
        with open("wdl_stdout.txt", "w+") as wdl_out:
            wdl_out.write(process.stdout)
    if process.stderr is not None:
        with open("wdl_stderr.txt", "w+") as wdl_err:
            wdl_out.write(process.stderr)

# Get the actual file path in gcloud from the template for vcf.
def get_file_path_from_template(template, google_project):
    process = subprocess.run(["gsutil", "-u", google_project, "ls", template],
                         stdout=subprocess.PIPE)
    filename = process.stdout.decode('utf-8').strip()
    return filename

def run_merge_command(num_batches, output_name, output_parent_dir):

    # Write options file including output directory in the bucket.
    options_json = "aou_merge_options.json"
    bucket = os.getenv("WORKSPACE_BUCKET")
    google_project = os.getenv("GOOGLE_PROJECT")
    output_path_gcloud = os.path.join(bucket, "saraj", output_parent_dir, "merged_outputs")

    # Config file includes gcloud file system setup.
    #config_file = "/home/jupyter/cromwell.conf"
    # Number of concurrent jobs is set to 25 up from 10.
    config_file = "cromwell.conf"

    # Write input json file based on selected sample(s) bam files.
    input_json = "aou_merge_inputs.json"
    # Write output directory in options file.
    write_options_json(options_json_filename=options_json,
                   output_path_gcloud=output_path_gcloud)
    # Add a template for each merged file resulting from one batch.
    # The wildcard is placed so that we do not have to know the
    # cromwell workflow root directory.
    batch_vcf_filenames = []
    batch_vcf_index_filenames = []
    for batch_idx in range(num_batches):
        # Get the filename for vcf file
        template = "{}".format(bucket) + \
                "/saraj/{}/{}_{}/".format(output_parent_dir, output_name, batch_idx) + \
                "run_advntr/*/call-sort_index/merged_samples.sorted.vcf.gz"
        filename = get_file_path_from_template(template=template,
                google_project=google_project)
        if len(filename) > 0:
            batch_vcf_filenames.append(filename)
        else:
            print("VCF file not added to the input file because file does not exist for template ",
                    template)
        # Get the filename for vcf index file
        template = "{}".format(bucket) + \
                "/saraj/{}/{}_{}/".format(output_parent_dir, output_name, batch_idx) + \
                "run_advntr/*/call-sort_index/merged_samples.sorted.vcf.gz.tbi"
        filename = get_file_path_from_template(template=template,
                google_project=google_project)
        if len(filename) > 0:
            batch_vcf_index_filenames.append(filename)
        else:
            print("VCF Index file not added to the input file because file does not exist for template ",
                    template)
    # Gcloud token is updated every time write_input_json is called.
    write_input_merge_json(input_json_filename=input_json,
                     batch_vcf_filenames=batch_vcf_filenames,
                     batch_vcf_index_filenames=batch_vcf_index_filenames)

    workflow_command = [
               "java", "-jar",
               "-Dconfig.file={}".format(config_file),
              "cromwell-86.jar",
              "run",
              "advntr_merge.wdl",
              "--inputs", "{}".format(input_json),
              "--options", "{}".format(options_json),
              ]
    run_single_command(workflow_command)

def is_output_dir_empty(directory, google_project):
    base_dir = os.path.basename(directory)
    parent_directory = os.path.dirname(directory)
    # If the directory does not exist, we'll get an exception.
    # Therefore, we re-route the stderr to /dev/null to not see the exception
    # in the output which might cause confusion.
    #try:
    process = subprocess.run(["gsutil", "-u", google_project, "ls", parent_directory],
                         stdout=subprocess.PIPE)
    #except Exception as e:
    #    print("caught exception")
    #    pass
    filenames = process.stdout.decode('utf-8')
    for filename in filenames.split():
        # Removing the trailing slashes, so we can get the directory names 
        # using basename function later.
        filename = os.path.normpath(filename.strip())
        if os.path.basename(filename) == base_dir:
                return False
    return True

# Calls the wdl workflow based on input sample filenames
def run_genotype_command(target_samples_df, output_name, output_parent_dir, region, vntr_id):

    # Write options file including output directory in the bucket.
    options_json = "aou_genotype_options.json"
    bucket = os.getenv("WORKSPACE_BUCKET")
    google_project = os.getenv("GOOGLE_PROJECT")
    parent_path_gcloud = os.path.join(bucket, "saraj", output_parent_dir)
    output_path_gcloud = os.path.join(parent_path_gcloud, output_name)

    if not is_output_dir_empty(directory=parent_path_gcloud,
                               google_project=google_project):
        print("Error: Output directory should be empty before running the workflow "+ \
              "in order to merge the files from the latest run properly. " + \
              "i.e. make sure the vcf files from previous runs are not involved " + \
              "when merging vcf files from this run.")
        exit(1)
    # Config file includes gcloud file system setup.
    config_file = "/home/jupyter/cromwell.conf"

    # Write input json file based on selected sample(s) bam files.
    input_json = "aou_genotype_inputs.json"
    num_batches = get_num_batches(args)
    for batch_idx in range(num_batches):
        start_time = datetime.now()
        # Get the indexes of samples in the current batch.
        # Then call the wdl workflow for only one batch at a time.
        first_sample_idx = batch_idx * args.batch_size
        last_sample_idx = min((batch_idx + 1) * args.batch_size, args.sample_count)
        samples_files = list(target_samples_df['alignment_file'])[first_sample_idx:last_sample_idx]
        # Write output directory in options file.
        output_dir = output_path_gcloud + "_" + str(batch_idx)
        write_options_json(options_json_filename=options_json,
                       output_path_gcloud=output_dir)
        # Gcloud token is updated every time write_input_json is called.
        write_input_json(input_json_filename=input_json,
                         samples_files=samples_files,
                         region=region,
                         vntr_id=vntr_id,
                         google_project=google_project,
                         )

        workflow_command = [
                   "java", "-jar",
                   "-Dconfig.file={}".format(config_file),
                  "cromwell-86.jar",
                  "run",
                  "advntr.wdl",
                  "--inputs", "{}".format(input_json),
                  "--options", "{}".format(options_json),
                  ]
        run_single_command(workflow_command)
        duration = datetime.now() - start_time
        print("Running batch {} finished in time {}".format(batch_idx, duration))
        sleep(60)

def parse_input_args():
    parser = argparse.ArgumentParser(
                    prog='advntr_wdl_aou',
                    description='Runs adVNTR wdl workflow based on provided '+ \
                                'input files and output directory. ' + \
                                'Currently working with lrWGS samples.',
                    )
    parser.add_argument("--dataset", required=True,
                        choices=["lrwgs", "srwgs"],
                        type=str,
                        help="Which dataset to get the aligned reads from.")

    parser.add_argument("--sample-count",
                        required=True,
                        type=int,
                        help="The number of samples randomly selected from the " + \
                             "lrWGS data. Max value on v7 data is 1027.")

    parser.add_argument("--batch-size",
                        default=20,
                        type=int,
                        help="The batch size for running genotype on. " +\
                             "The batch size should be such that downloading " + \
                             "inputs can finish under 1 hour that is when gcloud " + \
                             "token expires [default=20].")
    parser.add_argument("--output-name",
                        required=True,
                        type=str,
                        help="The output directory name for this experiment.")
    parser.add_argument("--region",
                         type=str, required=True,
                         help="The region(s) where VNTRs are located +- 1000bp. " + \
                              "This is used to stream targeted region of the BAM file " + \
                              "from the original location. This adds to efficiency. ")
    parser.add_argument("--vntr-id",
                         type=str, required=True,
                         help="The VNTR id used for this analysis. " + \
                              "VNTR id should correspond to the region flag otherwise" + \
                              "there will be no spanning reads.")
    args = parser.parse_args()
    return args

def get_num_batches(args):
    from math import ceil
    return ceil(args.sample_count / args.batch_size)

if __name__ == "__main__":
    # Parse input arguments.
    args = parse_input_args()

    # Find the selected set of samples we want to run the workflow on.
    target_samples = select_samples(args.sample_count, dataset=args.dataset)
    # For test run on local server
    #target_samples = pd.DataFrame({'grch38-bam':["./HG00438.bam"]})
    output_parent_dir = "batch_genotyping/run_8"
    # Run WDL workflow based on input files and output name.
    run_genotype_command(target_samples_df=target_samples,
                    output_name=args.output_name,
                    output_parent_dir=output_parent_dir,
                    region=args.region,
                    vntr_id=args.vntr_id)
    num_batches = get_num_batches(args)
    run_merge_command(num_batches=num_batches,
                      output_parent_dir=output_parent_dir,
                      output_name=args.output_name)
