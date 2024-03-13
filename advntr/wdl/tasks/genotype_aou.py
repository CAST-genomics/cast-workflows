import os
import subprocess
import argparse

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

# Write input json file based on selected sample(s) bam files.
def write_input_json(input_json_filename,
                     samples_files,
                     region,
                     vntr_id,
                     output_dir,
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
    data = {"jes_gcs_root": output_path_gcloud}
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


# Calls the wdl workflow based on input sample filenames
def run_wdl_command(target_samples_df, output_name, region, vntr_id):

    cromwell_version = "86"

    # Write options file including output directory in the bucket.
    options_json = "aou_options.json"
    bucket = os.getenv("WORKSPACE_BUCKET")
    google_project = os.getenv("GOOGLE_PROJECT")
    output_path_gcloud = os.path.join(bucket, "saraj", output_name)
    write_options_json(options_json_filename=options_json,
                       output_path_gcloud=output_path_gcloud)

    # Config file includes gcloud file system setup.
    config_file = "/home/jupyter/cromwell.conf"

    # Get sample id from the BAM file.
    #sample_id = os.path.basename(target_sample["grch38-bam"]).replace(".bam", "").replace(".cram", "")

    # Write input json file based on selected sample(s) bam files.
    input_json = "aou_inputs.json"
    samples_files = list(target_samples_df['alignment_file'])
    write_input_json(input_json_filename=input_json,
                     samples_files=samples_files,
                     region=region,
                     vntr_id=vntr_id,
                     output_dir=output_path_gcloud,
                     google_project=google_project,
                     )

    workflow_command = [
               "java", "-jar",
               "-Dconfig.file={}".format(config_file),
              "cromwell-{}.jar".format(cromwell_version),
              "run",
              "advntr.wdl",
              "--inputs", "{}".format(input_json),
              "--options", "{}".format(options_json),
              ]
    run_single_command(workflow_command)

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

if __name__ == "__main__":
    # Parse input arguments.
    args = parse_input_args()

    # Find the selected set of samples we want to run the workflow on.
    target_samples = select_samples(args.sample_count, dataset=args.dataset)
    # For test run on local server
    #target_samples = pd.DataFrame({'grch38-bam':["./HG00438.bam"]})

    # Run WDL workflow based on input files and output name.
    run_wdl_command(target_samples_df=target_samples,
                    output_name=args.output_name,
                    region=args.region,
                    vntr_id=args.vntr_id)
