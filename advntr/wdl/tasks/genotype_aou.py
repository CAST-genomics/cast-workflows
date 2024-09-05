import os
import subprocess
import argparse
import shutil
from time import sleep
from datetime import datetime

import json
import pandas as pd
import tempfile


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
    parser.add_argument("--only-batch-index",
                         required=False,
                         type=int,
                         help="Only run one batch at a time and exit (0-indexed). " + \
                              "This is useful when running a large batch " + \
                              "and not wanting to run everything at the same time. " + \
                              "This feature is optional. If not provided all batches will run.")
    parser.add_argument("--output-name",
                        required=True,
                        type=str,
                        help="The output directory name for this experiment.")
    parser.add_argument("--region-file",
                         type=str, required=True,
                         help="The region(s) where VNTRs are located +- 1000bp. " + \
                              "This is used to stream targeted region of the BAM file " + \
                              "from the original location. This adds to efficiency. ")
    parser.add_argument("--vntr-id",
                         type=str, required=True,
                         help="The VNTR id used for this analysis. " + \
                              "VNTR id should correspond to the region flag otherwise" + \
                              "there will be no spanning reads. Use comma separated ids or " + \
                              "'ALL' for all VNTRs in the database.")
    parser.add_argument("--mem",
                         type=int, required=False, default=4,
                         help="Memory allocation on the wdl nodes[4].")
    parser.add_argument("--cromwell",
                         action="store_true", required=False,
                         help="Run with cromwell instead of the default cromshell.")
    args = parser.parse_args()
    return args

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


# Write the input json file for when results from all batches are being merged.
def write_input_merge_json(input_json_filename,
                     batch_vcf_filenames,
                     batch_vcf_index_filenames,
                     mem,
                        ):
    data = {"run_merge_advntr.per_batch_vcfs": batch_vcf_filenames,
            "run_merge_advntr.per_batch_vcf_indexes": batch_vcf_index_filenames,
            "run_merge_advntr.mem": mem}

    with open(input_json_filename, "w+") as input_json_file:
        json.dump(data, input_json_file)

# Write input json file based on selected sample(s) bam files.
def write_input_json(input_json_filename,
                     samples_files,
                     region_file,
                     vntr_id,
                     mem,
                     google_project):

    data = {"run_advntr.bam_files": samples_files,
            "run_advntr.region_file": region_file,
            "run_advntr.vntr_id": vntr_id,
            "run_advntr.mem": mem,
            "run_advntr.google_project": google_project}

    with open(input_json_filename, "w+") as input_json_file:
        json.dump(data, input_json_file)

# Write options file including output directory in the bucket.
def write_options_json(options_json_filename, output_path_gcloud):
    data = {"jes_gcs_root": output_path_gcloud,
            "workflow_failure_mode": "ContinueWhilePossible",
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

def run_workflow(wdl_file, json_file,
                  json_options_file, config_file,
                  wdl_dependencies_file,
                  cromwell, dryrun=False):
    """
    Run workflow on AoU

    Arguments
    ---------
    json_file : str
        JSON file path with input arguments
    json_options_file : str
        JSON with additional options for cromshell

    dryrun : bool
        Just print the command, don't actually run cromshell
    """
    if cromwell is False:
        cmd = "cromshell submit {wdl} {json} -op {options}".format(
                        wdl=wdl_file,
                        json=json_file, options=json_options_file)
        if wdl_dependencies_file and wdl_dependencies_file.strip() != "":
            cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
    else:
        cmd = "java -jar -Dconfig.file={} ".format(config_file) + \
                  "jar_files/cromwell-87.jar run {} ".format(wdl_file) + \
                  "--inputs {} --options {}".format(json_file, json_options_file)
    if dryrun:
        sys.stderr.write("Run: %s\n"%cmd)
        return
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    print(output.decode("utf-8"))


def zip_dependencies(wdl_dependencies_file, files):
    """
    Put all WDL dependencies into a zip file

    Arguments
    ---------
    wdl_dependencies_file : str
        Zip file to put other wdls in
    """
    dirname = tempfile.mkdtemp()
    for f in files:
        shutil.copyfile("./%s"%f, dirname+"/"+f)
    shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)


# Get the actual file path in gcloud from the template for vcf.
def get_file_path_from_template(template, google_project):
    process = subprocess.run(["gsutil", "-u", google_project, "ls", template],
                         stdout=subprocess.PIPE)
    filename = process.stdout.decode('utf-8').strip()
    return filename

def run_merge_command(num_batches, num_samples,
                      batch_size, only_batch_index,
                      output_name, output_parent_dir,
                      mem, cromwell):
    merge_batches_individually = False
    # Write options file including output directory in the bucket.
    options_json = "aou_merge_options.json"
    bucket = os.getenv("WORKSPACE_BUCKET")
    google_project = os.getenv("GOOGLE_PROJECT")
    output_path_gcloud = os.path.join(bucket, "saraj", output_parent_dir, "merged_outputs")

    if merge_batches_individually and not is_output_dir_empty(directory=output_path_gcloud,
                               google_project=google_project):
        print("Error: Output merge directory should be empty before running the workflow "+ \
              "in order to merge the files from the latest run properly. "
              )
        exit(1)

    # Config file includes gcloud file system setup.
    config_file = "/home/jupyter/cromwell.conf"
    # Number of concurrent jobs is set to 25 up from 10.
    #config_file = "cromwell.conf"

    # Write input json file based on selected sample(s) bam files.
    input_json = "aou_merge_inputs.json"
    wdl_file = "advntr_merge.wdl"
    # Write output directory in options file.
    write_options_json(options_json_filename=options_json,
                   output_path_gcloud=output_path_gcloud)
    # Add a template for each merged file resulting from one batch.
    # The wildcard is placed so that we do not have to know the
    # cromwell workflow root directory.
    batch_vcf_filenames = []
    batch_vcf_index_filenames = []
    if merge_batches_individually:
        for batch_idx in range(num_batches):
            if batch_idx != only_batch_index:
                continue
            # Get the filename for vcf file
            #template = "{}".format(bucket) + \
            #        "/saraj/{}/{}_{}/".format(output_parent_dir, output_name, batch_idx) + \
            #        "run_advntr/*/call-merge_sort/merged_samples.sorted.vcf.gz"
            for sample_idx_in_batch in range(batch_size):
                if batch_idx * batch_size + sample_idx_in_batch > num_samples:
                    break
                template = "{}".format(bucket) + \
                    "/saraj/{}/{}_{}/".format(output_parent_dir, output_name, batch_idx) + \
                    "run_advntr/*/call-advntr_single_sample/shard-{}/".format(sample_idx_in_batch) + \
                    "advntr_single_sample/*/call-sort_index/*.sorted.vcf.gz"
                filename = get_file_path_from_template(template=template,
                        google_project=google_project)
                if len(filename) > 0:
                    batch_vcf_filenames.append(filename)
                    batch_vcf_index_filenames.append(filename + ".tbi")
                else:
                    print("VCF file not added to the input file because file does not exist for template ",
                            template)
    else:
        # All batches are individually merged. Now merge at the next level
        for batch_idx in range(num_batches):
            template = "{}".format(bucket) + \
                    "/saraj/vntr_reference_panel/batches/" + \
                    "merged_samples_p_g_vntrs_batch_{}.sorted.vcf.gz".format(batch_idx)
            filename = get_file_path_from_template(template=template,
                    google_project=google_project)
            if len(filename) > 0:
                batch_vcf_filenames.append(filename)
                batch_vcf_index_filenames.append(filename + ".tbi")
            else:
                print("VCF file not added to the input file because file does not exist for template ",
                        template)

    # Gcloud token is updated every time write_input_json is called.
    write_input_merge_json(input_json_filename=input_json,
                     batch_vcf_filenames=batch_vcf_filenames,
                     batch_vcf_index_filenames=batch_vcf_index_filenames,
                     mem=mem)


    # Run workflow
    run_workflow(wdl_file=wdl_file,
                    json_file=input_json,
                    json_options_file=options_json,
                    config_file=config_file,
                    wdl_dependencies_file=None,
                    cromwell=cromwell)



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
def run_genotype_command(target_samples_df, output_name,
                         output_parent_dir, region_file,
                         vntr_id, cromwell,
                         batch_size, sample_count,
                         mem, only_batch_index,
                         ):

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
    wdl_file = "advntr.wdl"

    num_batches = get_num_batches(sample_count=sample_count,
                                   batch_size=batch_size)

    for batch_idx in range(num_batches):
        if only_batch_index is not None and batch_idx != only_batch_index:
            continue
        start_time = datetime.now()
        # Get the indexes of samples in the current batch.
        # Then call the wdl workflow for only one batch at a time.
        first_sample_idx = batch_idx * batch_size
        last_sample_idx = min((batch_idx + 1) * batch_size, sample_count)
        samples_files = list(target_samples_df['alignment_file'])[first_sample_idx:last_sample_idx]
        # Write output directory in options file.
        output_dir = output_path_gcloud + "_" + str(batch_idx)
        write_options_json(options_json_filename=options_json,
                       output_path_gcloud=output_dir)
        # Gcloud token is updated every time write_input_json is called.
        write_input_json(input_json_filename=input_json,
                         samples_files=samples_files,
                         region_file=region_file,
                         vntr_id=vntr_id,
                         mem=mem,
                         google_project=google_project,
                         )

        # Zip all the WDL depencies
        wdl_dependencies_file = output_name + "advntr-wdl.zip"
        zip_dependencies(wdl_dependencies_file=wdl_dependencies_file,
                        files=["advntr_single.wdl"])

        # Run workflow
        run_workflow(wdl_file=wdl_file,
                        json_file=input_json,
                        json_options_file=options_json,
                        config_file=config_file,
                        wdl_dependencies_file=wdl_dependencies_file,
                        cromwell=cromwell)

        duration = datetime.now() - start_time
        print("Running batch {} finished in time {}".format(batch_idx, duration))


def get_num_batches(sample_count, batch_size):
    from math import ceil
    return ceil(sample_count / batch_size)

if __name__ == "__main__":
    # Parse input arguments.
    args = parse_input_args()

    # Find the selected set of samples we want to run the workflow on.
    target_samples = select_samples(args.sample_count, dataset=args.dataset)
    # For test run on local server
    #output_parent_dir = "batch_genotyping/run_p_vntrs_g_vntrs"
    output_parent_dir = "batch_genotyping/{}".format(args.output_name)
    run_batches = False
    merge_batches = not run_batches
    # Run WDL workflow based on input files and output name.
    if run_batches:
        run_genotype_command(target_samples_df=target_samples,
                    output_name=args.output_name,
                    output_parent_dir=output_parent_dir,
                    region_file=args.region_file,
                    vntr_id=args.vntr_id,
                    cromwell=args.cromwell,
                    batch_size=args.batch_size,
                    sample_count=args.sample_count,
                    mem=args.mem,
                    only_batch_index=args.only_batch_index,
                    )
    num_batches = get_num_batches(sample_count=args.sample_count,
                                   batch_size=args.batch_size)
    if merge_batches:
        run_merge_command(num_batches=num_batches,
                      num_samples=args.sample_count,
                      batch_size=args.batch_size,
                      output_parent_dir=output_parent_dir,
                      output_name=args.output_name,
                      cromwell=args.cromwell,
                      only_batch_index=args.only_batch_index,
                      mem=args.mem)
