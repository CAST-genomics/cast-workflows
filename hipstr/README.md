# HipSTR Genotyping Workflow

# How to run HipSTR on multiple samples
  * Use hipstr_multi.wdl to run HipSTR on joint samples. This will generate a stutter_model.txt file.
  * Use stutter_model.txt generated from the hipstr_multi.wdl run as a input for hipstr_single.wdl run.

## Test workflow locally
  1. Download [Cromwell and womtool](https://github.com/broadinstitute/cromwell/releases/), [JDK 11](https://docs.oracle.com/en/java/javase/11/install/installation-jdk-macos.html#GUID-F575EB4A-70D3-4AB4-A20E-DBE95171AB5F).
  2. Clone the repository to your local computer.
  3. Run command `java -jar cromwell-84.jar run -i <path_to_input.json> <path_to_wdl_file.wdl>`,the default ouput directory for cromwell is `./cromwell-executions`. To change ouput directory, add `--options options.json` to the command line.
  * The options.json file can be used for change output directory.
