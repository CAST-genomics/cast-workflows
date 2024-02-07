import os
import subprocess
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
                      prog="extract_target_region",
                      description="Extract target regions from input bam file(s) over Gcloud." + \
                      " Provide updated index file." + \
                      " The target regions are later fed into AdVNTR for genotyping.")
    parser.add_argument("--google-project",
                        type=str,
                        required=True,
                        help="The GOOGLE_PROJECT environment variable obtained from the gcloud bucket.")
    parser.add_argument("--gcloud-token",
                        type=str,
                        required=True,
                        help="The gcloud token required to access data on the gcloud bucket.")
    parser.add_argument("--region",
                        type=str,
                        help="The region(s) where VNTRs are located.")
    parser.add_argument("--regions-file",
                        type=str,
                        help="The filename with one line per target region (VNTR). " + \
                        "If --region is set, this argument is ignored.")
    parser.add_argument("--input-bam",
                         type=str,
                         help="Path to the input bam file for one individual.")
    parser.add_argument("--input-bam-list",
                         type=str,
                         help="The filname with one path to a bam file per line. " + \
                         "If --input-bam is set, this argument is ignored.")
    args = parser.parse_args()

    # Check sample id is provided either for a single sample or as a file.
    if args.input_bam is None and args.input_bam_list is None:
        print("Error: Either --input-bam or --input-bam-list should be set to locate bam files.")
        exit(1)
    # Check region is provided either for a single region or as a file.
    if args.region is None and args.regions_file is None:
        print("Error: Either --region or --regions-file should be set to target regions from bam files.")
        exit(1)
    return args

def set_environment_variables(args):
    # Set flags for samtools to read over gcloud URL
    os.environ["HTSLIB_CONFIGURE_OPTIONS"] ="--enable-gs"
    os.environ["GCS_OAUTH_TOKEN"] = args.gcloud_token
    os.environ["GCS_REQUESTER_PAYS_PROJECT"] = args.google_project

def run_single_command(command):
    process = subprocess.call(command.split(" "))
    return
    if process.stdout is not None:
        with open("wdl_stdout.txt", "w+") as wdl_out:
            wdl_out.write(process.stdout)
    if process.stderr is not None:
        with open("wdl_stderr.txt", "w+") as wdl_err:
            wdl_out.write(process.stderr)

def get_single_target_bam_file(args, input_bam):
    sample_id = os.path.basename(input_bam).replace(".bam", "")
    command = ""
    if args.region is not None:
        command = "samtools view -h -b -o target_region_{}_unsorted.bam {} {}".format(
                    sample_id, input_bam, args.region)
    else:
        # if --region is not set, --regions-file should be set.
        # Otherwise the program ends with an error before this function is called.
        # Note -M -L FILE is equivalent to --use-index --target-file FILE
        # and is equivalent to --region-file FILE and --regions-file FILE.
        command = "samtools view -h -b -o target_region_{}_unsorted.bam -M -L {} {}".format(
                    sample_id, args.regions_file, input_bam)
    run_single_command(command)
    command = "samtools sort -o target_region_{}.bam target_region_{}_unsorted.bam".format(
                        sample_id, sample_id)
    run_single_command(command)
    command = "samtools index target_region_{}.bam".format(sample_id)
    run_single_command(command)

def get_target_bam_files(args):
    if args.input_bam is not None:
        # Working with a single bam file.
        get_single_target_bam_file(args, args.input_bam)
    else:
        # if --input-bam is not set, --input-bam-list should be set.
        # Otherwise the program ends with an error before this function is called.
        input_bam_list_filename = args.input_bam_list
        input_bams = []
        with open(input_bam_list_filename, "r") as input_bam_list_file:
            input_bams = input_bam_list_file.readlines()
            input_bams = [line.strip() for line in input_bams]
        for input_bam in input_bams:
            get_single_target_bam_file(args, input_bam)

if __name__ == "__main__":
    args = parse_args()
    set_environment_variables(args)
    get_target_bam_files(args)
