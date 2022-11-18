#!/bin/bash

# run the following to check the syntax of a wdl file:
#java -jar ~/packages/wdl/womtool-47.jar validate hello.wdl

java -jar ~/packages/wdl/cromwell-47.jar run advntr.wdl > wdl_output
