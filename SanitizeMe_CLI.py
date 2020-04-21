#!/usr/bin/env python

import sys
import os
import glob
import re
import argparse 

def main():
    cli = argparse.ArgumentParser()
    
    cli.add_argument('-i', '--InputFolder', help="Folder containing fastq files. Only files ending in .fq, .fg.gz, .fastq, and .fastq.gz will be processed", required=True)
    cli.add_argument('-r', '--Reference', help="Host Reference fasta or fasta.gz file", required=True)
    cli.add_argument('-o', '--OutputFolder', help="Output Folder. Default is ~/dehost_output/test", required=False, default='~/dehost_output/test')

    cli.add_argument('-t', '--threads', help="Number of threads. Default is 4. More is faster if your computer supports it", type=int, required=False, default=4)
    method = cli.add_mutually_exclusive_group()
    method.add_argument('--Nanopore', help="Select if you used Nanopore Sequencing", action='store_const', dest='seq_method', const='map-ont', default='map-ont')
    method.add_argument('--PacBio', help="Select if you used PacBio Genonmic Reads", action='store_const', dest='seq_method', const='map-pb')
    method.add_argument('--PacBioCCS', help="Select if you used PacBio CCS Genomic Reads", action='store_const', dest='seq_method', const='asm20')
    method.add_argument('--ShortRead', help="Select if you have single end short reads (Illumina)", action='store_const', dest='seq_method', const='sr')
    args = cli.parse_args()

    files = sorted([f for f in glob.glob(args.InputFolder+"/*") if re.search(r'(.*)\.((fastq|fq)(|\.gz))$', f)])
    OutputFolder = os.path.expanduser(args.OutputFolder)
    os.system(f"mkdir -p {OutputFolder}")
    f=open(f"{OutputFolder}/cmd.log", 'w+')

    for j in range(0, len(files)):
        i = files[j]
        base = os.path.splitext(os.path.basename(i))[0]
        #print(base)
        os.system(f"mkdir -p {OutputFolder}")
        minimap2_cmd = f"minimap2 -ax {args.seq_method} {args.Reference} {i} -t {args.threads} > {OutputFolder}/{base}.sam"
        f.write(minimap2_cmd+'\n')
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -u -f 4 {OutputFolder}/{base}.sam > {OutputFolder}/{base}_filtered.sam"
        f.write(samtools_cmd1+'\n')
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools bam2fq {OutputFolder}/{base}_filtered.sam > {OutputFolder}/{base}_filtered.fastq"
        f.write(samtools_cmd2+'\n')
        os.system(samtools_cmd2)
        delete_cmd1 = f"rm {OutputFolder}/{base}.sam"
        os.system(delete_cmd1)
        f.write(delete_cmd1+'\n')
        delete_cmd2 = f"rm {OutputFolder}/{base}_filtered.sam"
        os.system(delete_cmd2)
        f.write(delete_cmd2+'\n')

        print("progress: {}/{}".format(j+1, len(files)))

    f.close()

if __name__ == "__main__":
    sys.exit(main())
