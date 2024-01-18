#!/bin/bash

# Updated version cromwell is v86 as of 1/11/2024
java  -jar -Dconfig.file=/home/jupyter/cromwell.conf cromwell-86.jar run advntr.wdl --inputs aou_inputs.json &> wdl_output.txt
