#!/bin/bash

# run the following to check the syntax of a wdl file:
#java -jar ~/packages/wdl/womtool-47.jar validate hello.wdl

# Updated version as of 1/11/2024
java  -jar ~/packages/cromwell-86/cromwell-86.jar run advntr.wdl --inputs genomequery_inputs.json &> wdl_output.txt
