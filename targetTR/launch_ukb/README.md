# TargetTR on UKB RAP

## Setup WDL workflow

First, compile WDL workflows to run on DNA Nexus

```
# Note: Checks are more rigorous than those of womtool validate
java -jar "$DX_COMPILER_JAR" compile ../wdl/workflows/targetTR.wdl \
	-project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /TargetedSTR/ -streamFiles all -archive -separateOutputs

java -jar "$DX_COMPILER_JAR" compile ../wdl/workflows/merge_index_str.wdl \
  -project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /TargetedSTR/ -streamFiles all -archive 

java -jar "$DX_COMPILER_JAR" compile ../wdl/tasks/associatr.wdl \
  -project project-GG25fB8Jv7B928vqK7k6vYY6 -folder /STRPhewas/ -streamFiles all -archive
```

## Compile list of files to process

List all UKB cram/crai files
```
for i in $(seq 10 60)
do
  cmd="dx ls -l 'Bulk/Whole genome sequences/Whole genome CRAM files/${i}/'"
  sh -c "${cmd}"
done > ukb_cram_files_long.txt

./process_cram_list.py ukb_cram_files_long.txt > ukb_cram_and_index_files.txt
```

## Call python launcher 

The python launcher takes care of batching, which was difficult to do cleanly in WDL 1.0.

The code below shows a small example (3 batches of 25 samples each).
To run all files, increase `--batch-size` to a larger number (1000?) and remove the `--batch-num` argument.
```
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 25 \
  --batch-num 3 \
  --workflow-id workflow-GPfbXV8Jv7B27kpf6Y50QyQ9 \
  --file-list ukb_cram_and_index_files.txt
```

See `launch_scripts/` for full launches.

## Call phewas (under development...)

```
dx run -y -f GATM-EGFR-test.json --destination /STRPhewas/results/ workflow-GPbFzQ8Jv7BGF8fp4k4fbFfB
```

## associatr wishlist

* Take plink style phenotypes file as input, and don't require numerical sample IDs
* Change to not use positional arguments