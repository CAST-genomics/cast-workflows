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
def select_samples(sample_count):
    manifest_file = "manifest_lrwgs.csv"
    if not os.path.exists(manifest_file):
        print("Warning: The manifest file for long reads WGS is not present or not correctly named.")
        download_manifest_file(manifest_file)

    files_df = pd.read_csv(manifest_file)
    num_individuals = len(files_df)
    print("Working with manifest file of {} columns and {} samples".format(
            len(files_df.columns), num_individuals))
    selected_samples = files_df.sample(n=sample_count, random_state=0)
    print("Selected {} random sample(s)".format(sample_count))
    return selected_samples[["grch38-bam", "grch38-bai"]]

# Write input json file based on selected sample(s) bam files.
def write_input_json(input_json_filename, sample_df):
    data = {"run_advntr.bam_file": sample_df["grch38-bam"],
            "run_advntr.bam_index": sample_df["grch38-bai"]}
    with open(input_json_filename, "w+") as input_json_file:
        json.dump(data, input_json_file)

# Write options file including output directory in the bucket.
def write_options_json(options_json_filename, output_path_gcloud):
    data = {"jes_gcs_root": output_path_gcloud}
    with open(options_json_filename, "w+") as options_json_file:
        json.dump(data, options_json_file)

# Calls the wdl workflow based on input sample filenames
def run_wdl_command(target_samples_df, output_name):

    cromwell_version = "86"

    # Write options file including output directory in the bucket.
    options_json = "aou_options.json"
    bucket = os.getenv("WORKSPACE_BUCKET")
    output_path_gcloud = os.path.join(bucket, "saraj", output_name)
    write_options_json(options_json_filename=options_json,
                       output_path_gcloud=output_path_gcloud)

    # Config file includes gcloud file system setup.
    config_file = "/home/jupyter/cromwell.conf"

    # TODO: edit the WDL file to work with a batch instead of a single file.

    for idx, target_sample in target_samples_df.iterrows():

        # Write input json file based on selected sample(s) bam files.
        input_json = "aou_inputs.json"
        write_input_json(input_json_filename=input_json,
                         sample_df=target_sample)
        command = ["java", "-jar",
                   "-Dconfig.file={}".format(config_file),
                  "cromwell-{}.jar".format(cromwell_version),
                  "run",
                  "advntr.wdl",
                  "--inputs", "{}".format(input_json),
                  "--options", "{}".format(options_json),
                  ]
        process = subprocess.run(command)
        if process.stdout is not None:
            with open("wdl_stdout.txt", "w+") as wdl_out:
                wdl_out.write(process.stdout)
        if process.stderr is not None:
            with open("wdl_stderr.txt", "w+") as wdl_err:
                wdl_out.write(process.stderr)

def parse_input_args():
    parser = argparse.ArgumentParser(
                    prog='advntr_wdl_aou',
                    description='Runs adVNTR wdl workflow based on provided '+ \
                                'input files and output directory. ' + \
                                'Currently working with lrWGS samples.',
                    )
    parser.add_argument("--sample-count",
                        required=True,
                        type=int,
                        help="The number of samples randomly selected from the " + \
                             "lrWGS data. Max value on v7 data is 1027.")

    parser.add_argument("--output-name",
                        required=True,
                        type=str,
                        help="The output directory name for this experiment.")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    # Parse input arguments.
    args = parse_input_args()

    # Find the selected set of samples we want to run the workflow on.
    target_samples = select_samples(args.sample_count)

    # Run WDL workflow based on input files and output name.
    run_wdl_command(target_samples_df=target_samples,
                    output_name=args.output_name)
