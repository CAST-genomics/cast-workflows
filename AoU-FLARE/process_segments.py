#!/usr/bin/env python3

import argparse

def process_segments(input_file, output_file):
    # Initialize mappings
    cols_to_samp = {}
    samp_to_segments = {}

    # Open the input file
    with open(input_file, "r") as f:
        # Read the header line
        header = f.readline().strip().split("\t")
        
        # Parse the header for sample names and haplotypes
        for i in range(2, len(header)):
            sample, hap = header[i].split("-")
            cols_to_samp[i] = sample
            if sample not in samp_to_segments:
                samp_to_segments[sample] = [[], []]  # [segment labels, segment lengths]

        # Initialize previous state tracking variables
        prev_starts = [None] * len(header)
        prev_labels = [-1] * len(header)

        # Process each row in the file
        for line in f:
            row = line.strip().split("\t")
            current_pos = int(row[1])  # Current position (from the 'pos' column)
            
            for col in range(2, len(header)):  # Starting from the third column (the data columns)
                current_label = int(row[col])  # Current label for this column
                
                # If it's the first row, initialize the starts and labels
                if prev_starts[col] is None:
                    prev_starts[col] = current_pos
                    prev_labels[col] = current_label
                else:
                    # If we are still in the same segment, do nothing
                    if current_label == prev_labels[col]:
                        continue
                    else:
                        # Otherwise, we have a change in label; close the previous segment and start a new one
                        seg_start = prev_starts[col]
                        seg_end = current_pos - 1  # The segment ends at the previous position
                        seg_label = prev_labels[col]
                        sample = cols_to_samp[col]
                        
                        # Calculate segment length (ensure it's non-negative)
                        seg_length = seg_end - seg_start + 1
                        if seg_length < 0:
                            raise ValueError(f"Negative segment length detected: {seg_length} for sample {sample}")
                        
                        # Add the segment for the sample
                        samp_to_segments[sample][0].append(seg_label)
                        samp_to_segments[sample][1].append(seg_length)
                        
                        # Reset the previous states for the new segment
                        prev_starts[col] = current_pos
                        prev_labels[col] = current_label

    # At the end of the loop, add the last segment for each sample
    for col in range(2, len(header)):
        seg_label = prev_labels[col]
        sample = cols_to_samp[col]
        seg_start = prev_starts[col]
        seg_end = current_pos  # Use the last position from the loop
        seg_length = seg_end - seg_start + 1  # +1 to include both start and end positions
        
        # Add the last segment for the sample
        samp_to_segments[sample][0].append(seg_label)
        samp_to_segments[sample][1].append(seg_length)

    # Save the results to the output file
    with open(output_file, "w") as f:
        f.write("Sample\tLabels\tLengths\n")
        for sample, (labels, lengths) in samp_to_segments.items():
            f.write(f"{sample}\t{','.join(map(str, labels))}\t{','.join(map(str, lengths))}\n")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process ancestry segments from a TSV file.")
    parser.add_argument("input_file", help="Path to the input TSV file")
    parser.add_argument("output_file", help="Path to the output TSV file")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the file and save the results
    process_segments(args.input_file, args.output_file)
