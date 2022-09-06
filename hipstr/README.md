# HipSTR Genotyping Workflow

# How to run HipSTR on multiple samples
  * Use hipstr_multi.wdl to run HipSTR on joint samples. This will generate a stutter_model.txt file.
  * Use stutter_model.txt generated from the hipstr_multi.wdl run as a input for hipstr_single.wdl run.
