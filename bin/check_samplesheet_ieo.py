#!/usr/bin/env python3
 
import os
import sys
import errno
import argparse
from pathlib import Path
import glob
import re

def parse_args(args=None):
    Description = "Reformat nf-core/chipseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet_ieo.py <FILE_IN> <FILE_OUT>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure: use csv comma separation:
    rid,sid,sample,replicate,path,lanes,is_input,which_input,antibody
    this design follow the logic of sample annotation at the IEO where each sample is assigned a unique SID id
    associated to each run of the specified sample is a RID
    """
    sample_mapping_dict = {}
    with open(file_in, "r", encoding='utf-8-sig') as fin:
        ## Check header
        MIN_COLS = 9
        HEADER = ["rid","sid","sample","replicate","path","lanes","is_input","which_input","antibody"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)
        ## Check sample entries else line=fin.readline()
        # to debug the file we can use the following:
        # line=fin.readline()
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]
            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                    "Line",
                    line,
                )
            ## Check sample name entries
            rid,sid,sample,replicate,path,lanes,is_input,which_input,antibody = lspl[: len(HEADER)]
            sample=sample + "_" + replicate  # create a unique identifie
            
            print(sample)
            
            if sample.find(" ") != -1:
                print(
                    f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                )
                sample = sample.replace(" ", "_")
            if sample.find("-") != -1:
                print(
                    f"WARNING: - have been replaced by underscores for sample: {sample}"
                )
                sample = sample.replace("-", "_")    
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)
                
            if which_input.find(" ") != -1:
                print(
                    f"WARNING: Spaces have been replaced by underscores for which_input: {which_input}"
                )
                which_input = which_input.replace(" ", "_")
            if which_input.find("-") != -1:
                print(
                    f"WARNING: - have been replaced by underscores for which_input: {which_input}"
                )
                which_input = which_input.replace("-", "_")    

            # According to IEO file path policy samples are folders :
            # "..._rid_sid/" within samples are with extension: ".fastq.gz" and again rid_sid_..._R1...fastq.gz and rid_sid_..._R2...fastq.gz
            # the default is paired end sequencing at IEO
            # Given the way the crisper screen library is designed R1 and R2 must be present! Check:
            p=Path(path)
            # here we consider the "lanes" in case issues associated to the NOVASeq lanes are present.
            # lanes can have values in any, or the L001, L002, ... 
            # Illumina outputs include: "_L002_R1_001.fastq.gz"
            ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, input ]]}
            lanes=lanes.strip().split(":")
            if any(x!='any' for (x) in lanes) and any(x.startswith('L00') for (x) in lanes) and any(x[-1].isdigit() for (x) in lanes):
                print("Proceed with specified lanes:")
                # in here we need to split by common partial match:
                # We need to get any files with the extension in either curent or any subfolder
                
                subfolders = list(p.glob("*"+rid+"*"+sid+"*"))
                
                search = '*fastq.gz'
                fastqs = list(subfolders[0].glob('**/' + search))
                if len(fastqs) == 0:
                    print_error(
                            f"The path provided does not contain the files as RID+SID",
                            "path:",
                            path
                        )
                
                for ll in lanes:
                    print(ll)
                    search = "*"+rid+"*"+sid+'*'+ll+'*R1*fastq.gz'
                    fastqs_1 = list(subfolders[0].glob('**/' + search))
                    search = "*"+rid+"*"+sid+'*'+ll+'*R2*fastq.gz'
                    fastqs_2 = list(subfolders[0].glob('**/' + search))
                    
                    # check if fastqs_1 is not empty
                    if len(fastqs_1) > 0:
                        print("Found", len(fastqs_1), "files matching the pattern.")
                    else:
                        print("No files found matching the pattern.")
                    # check if fastqs_2 is not empty
                    if len(fastqs_2) > 0:
                        print("Found", len(fastqs_2), "files matching the pattern.")
                    else:
                        print("No files found matching the pattern.")
                    
                    ## Auto-detect paired-end/single-end
                    sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
                    if sample and len(fastqs_1) !=0 and len(fastqs_2) !=0:  ## Paired-end short reads
                        sample_info = ["0", str(fastqs_1[0]), str(fastqs_2[0]), is_input,which_input,antibody]
                    elif sample and len(fastqs_1) !=0 and len(fastqs_2) ==0:  ## Single-end short reads
                        sample_info = ["1", str(fastqs_1[0]), str(fastqs_2[0]), is_input,which_input,antibody]
                    else:
                        print_error("Invalid combination of columns provided!", "Line", line)
                    
                    # for each we append:
                    if sample not in sample_mapping_dict:
                        sample_mapping_dict[sample] = [sample_info]
                        print(sample,sample_info)
                    else:
                        if sample_info in sample_mapping_dict[sample]: 
                            print_error("Samplesheet contains duplicate rows!", "Line", line)
                        else:
                            print("append")
                            sample_mapping_dict[sample].append(sample_info)
                            print(sample,sample_info)
                            
            elif all(x=='any' for (x) in lanes):
                print("Proceed with all lanes:")
                # we need to get to the unique elements of "Lanes"
                subfolders = list(p.glob("*"+rid+"*"+sid+"*"))
                search = '*R1*fastq.gz'
                fastqs = list(subfolders[0].glob('**/' + search)) # to find recursively in all subfolders use: '**/' + search
                fastqs = [ re.sub(r'^.*?_L', 'L', str(x)).split('_')[0] for x in fastqs]
                
                if len(fastqs) == 0:
                    print_error(
                            f"The path provided does not contain the files as RID+SID",
                            "path:",
                            path
                        )
                    
                for ll in fastqs:
                    print(ll)
                    search = "*"+rid+"*"+sid+'*'+ll+'*R1*fastq.gz'
                    fastqs_1 = list(subfolders[0].glob('**/' + search))
                    search = "*"+rid+"*"+sid+'*'+ll+'*R2*fastq.gz'
                    fastqs_2 = list(subfolders[0].glob('**/' + search))
                    
                    # check if fastqs_1 is not empty
                    if len(fastqs_1) > 0:
                        print("Found", len(fastqs_1), "files matching the pattern.")
                    else:
                        print("No files found matching the pattern.")
                    # check if fastqs_2 is not empty
                    if len(fastqs_2) > 0:
                        print("Found", len(fastqs_2), "files matching the pattern.")
                    else:
                        print("No files found matching the pattern.")
                    
                    ## Auto-detect paired-end/single-end
                    sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
                    if sample and len(fastqs_1) !=0 and len(fastqs_2) !=0:  ## Paired-end short reads
                        sample_info = ["0", str(fastqs_1[0]), str(fastqs_2[0]), is_input,which_input,antibody]
                    elif sample and len(fastqs_1) !=0 and len(fastqs_2) ==0:  ## Single-end short reads
                        sample_info = ["1", str(fastqs_1[0]), str(""), is_input,which_input,antibody]
                    else:
                        print_error("Invalid combination of columns provided!", "Line", line)
                    
                    # for each we append:
                    if sample not in sample_mapping_dict:
                        sample_mapping_dict[sample] = [sample_info]
                        print(sample,sample_info)
                    else:
                        if sample_info in sample_mapping_dict[sample]: 
                            print_error("Samplesheet contains duplicate rows!", "Line", line)
                        else:
                            print('append')
                            sample_mapping_dict[sample].append(sample_info)
                            print(sample,sample_info)
                
            else:
                print(lanes)
                print_error(
                    f"The lane  format must by either any or in the from L00* with multiple lanes separate by : check sample sheet!"
                )
    
    for key in sample_mapping_dict:
        print("Dict keys: ",key)
    
    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "single_end",
                        "fastq_1",
                        "fastq_2",
                        "is_input",
                        "which_input",
                        "antibody"
                     ]
                    )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample have identical attributes: should add a check for the control! that must match the ids
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error(
                        f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    control = val[-2] # the which_input
                    if control != 'false' and control not in sample_mapping_dict.keys():
                        print_error(
                            f"Control identifier has to match a provided sample identifier!",
                            "Control",
                            control
                        )
                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())