import sys, os
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(prog = "collect_phenotypes",
                description="collect phenotype ids, name and category based on an association file.")

    parser.add_argument("--input-file", help="The input association .tab from where phenotypes are extracted.",
        type=str, default="outputs_min_cases_1000/phewas_results_TMCO1_165761972.csv")
    parser.add_argument("--output", help="The name of the output file.",
        type=str, default="phenotypes_min_1000_cases.csv")

    args = parser.parse_args()
    return args


def read_tab_file(name):
    phenotypes = []
    with open(name, "r") as in_file:
        for l_idx, line in enumerate(in_file.readlines()):
            if l_idx == 0:
                # Skip the header line
                continue
            phecode = line.split(",")[0]
            category = line.split(",")[-1].strip().replace("/", "_").replace(" ", "_")
            # There might be commas in the middle of the names.
            # In that case, there are likely inside a (['".
            comma_counter = 0
            string = ""
            for idx, char in enumerate(line):
                if char == ",":
                    comma_counter += 1
                if comma_counter == 12:
                    break
            # Move to the start of the next column
            idx += 1
            stack = []
            while idx < len(line):
                if line[idx] in ["[", "("]:
                    stack.append(line[idx])
                elif line[idx] in ["]", ")"]:
                    stack.pop()
                elif line[idx] == '"':
                    if len(stack) == 0 or stack[-1] != '"':
                        stack.append(line[idx])
                    else:
                        stack.pop()
                if len(stack) == 0 and line[idx] == ",":
                    # Moved to the next column
                    # Check the remainder of the line is the same length as the category as a sanity check
                    if idx + len(category) + 1 != len(line) - 1:
                        print("For line {} string {} idx {} len category {} len line {} does not match the length of the line".format(
                                l_idx, string, idx, len(category), len(line)))
                    break
                else:
                    string += line[idx]
                idx += 1
            if idx == len(line):
                print("String not correctly processed for line {} {}".format(l_idx, line))
            # Replace all special characters with underline so it can be easily processed
            string = string.replace("'", "_").replace('"', "_").\
                        replace('(', '_').replace(")", "_").\
                        replace('[', '_').replace("]", "_").\
                        replace('/', '_').replace(" ", "_").\
                        replace(',', '_').replace('*', '_')
            #print("category: {} \t string: {}".format(category, string))
            phenotypes.append((phecode, string, category))
    return phenotypes

def write_phenotypes(phenotypes, filename):
    with open(filename, "w") as out_file:
        out_file.write("phecode,phecode_string,phecode_category\n")
        for pheno_line in phenotypes:
            out_file.write(",".join(pheno_line) + "\n")
            

if __name__ == "__main__":
    args = parse_arguments()

    phenotypes = read_tab_file(args.input_file)
    write_phenotypes(phenotypes, args.output)
