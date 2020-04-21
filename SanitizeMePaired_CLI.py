#!/usr/bin/env python

import sys
import os
import glob
import re
from colored import stylize, attr, fg
import argparse 

def main():
    cli = argparse.ArgumentParser()

    cli.add_argument('-i', '--InputFolder', help="Folder containing paired fq, fq.gz, fastq, and fastq.gz files. Program will recursively find paired reads", required=True)
    cli.add_argument('-r', '--Reference', help="Host Reference fasta or fasta.gz file", required=True)
    cli.add_argument('-o', '--OutputFolder', help="Output Folder. Default is ~/dehost_output/test", required=False, default='~/dehost_output/test')

    cli.add_argument('-t', '--threads', help="Number of threads. More is faster if your computer supports it", type=int, required=False, default=4)

    args = cli.parse_args()

    for_files = sorted([f for f in glob.glob(args.InputFolder+"/**", recursive = True) if re.search(r'(.*)_(R|)1(.*)\.((fastq|fq)(|\.gz))$', f)])
    rev_files = sorted([f for f in glob.glob(args.InputFolder+"/**", recursive = True) if re.search(r'(.*)_(R|)2(.*)\.((fastq|fq)(|\.gz))$', f)])
    OutputFolder = os.path.expanduser(args.OutputFolder)
    os.system(f"mkdir -p {OutputFolder}")
    f=open(f"{OutputFolder}/cmd.log", 'w+')

    if (len(for_files) != len(rev_files)):
        print(stylize(f"You have unequal numbers of forward and reverse files!", fg("red") + attr("bold")))
        raise Exception(stylize(f"You have {len(for_files)} forward files and {len(rev_files)} reverse files!", fg("red") + attr("bold")))

    for i in range(0, len(for_files)):

        #print(for_files[i])
        #print(rev_files[i])

        base = os.path.splitext(os.path.basename(for_files[i]))[0]
        #print(base)
        os.system(f"mkdir -p {OutputFolder}")
        minimap2_cmd = f"minimap2 -ax sr {args.Reference} {for_files[i]} {rev_files[i]} -t {args.threads} > {OutputFolder}/{base}.sam"
        f.write(minimap2_cmd+'\n')
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -u -f 4 {OutputFolder}/{base}.sam > {OutputFolder}/{base}_filtered.sam"
        f.write(samtools_cmd1+'\n')
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools bam2fq {OutputFolder}/{base}_filtered.sam > {OutputFolder}/{base}_filtered.fastq"
        f.write(samtools_cmd2+'\n')
        os.system(samtools_cmd2)

        split1_cmd = f"cat {OutputFolder}/{base}_filtered.fastq | grep '^@.*/1$' -A 3 --no-group-separator >{OutputFolder}/{base}_filtered_r1.fastq"
        split2_cmd = f"cat {OutputFolder}/{base}_filtered.fastq | grep '^@.*/2$' -A 3 --no-group-separator >{OutputFolder}/{base}_filtered_r2.fastq"

        os.system(split1_cmd)
        f.write(split1_cmd+'\n')
        os.system(split2_cmd)
        f.write(split2_cmd+'\n')

        delete_cmd1 = f"rm {OutputFolder}/{base}.sam"
        os.system(delete_cmd1)
        f.write(delete_cmd1+'\n')
        delete_cmd2 = f"rm {OutputFolder}/{base}_filtered.sam"
        f.write(delete_cmd2+'\n')
        os.system(delete_cmd2)
        delete_cmd3 = f"rm {OutputFolder}/{base}_filtered.fastq"
        os.system(delete_cmd3)
        f.write(delete_cmd3+'\n')

        print("progress: {}/{}".format(i+1, len(for_files)))

    f.close()

if __name__ == "__main__":
    sys.exit(main())
