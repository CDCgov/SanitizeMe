#!/usr/bin/env python

import sys
import os
import glob
import re
from gooey import Gooey, GooeyParser

@Gooey(default_size=(750, 820), 
       progress_regex=r"^progress: (?P<current>\d+)/(?P<total>\d+)$",
       progress_expr="current / total * 100")
        

def main():
    cli = GooeyParser(description="Remove Host Reads from Long Read, Single End, Fastq Files")
    required_args = cli.add_argument_group("Input Output", gooey_options={'show_border': True, 'columns': 1})
    required_args.add_argument('--InputFolder', help="Folder containing fastq files. Only files ending in .fq, .fg.gz, .fastq, and .fastq.gz will be processed", required=True, widget='DirChooser')
    required_args.add_argument('--Reference', help="Host Reference fasta or fasta.gz file", required=True, widget='FileChooser')
    required_args.add_argument('--OutputFolder', help="Output Folder", required=False, default='~/dehost_output/test')

    parser = cli.add_argument_group("Options", gooey_options={'show_border': True,'columns': 1})
    parser.add_argument('--threads', help="Number of threads. More is faster if your computer supports it", type=int, required=False, default=4)

    method = parser.add_mutually_exclusive_group("Sequencing Method", gooey_options={'show_border': True})
    method.add_argument('--Nanopore', help = 'Data from Nanopore Sequencing', action='store_true')
    method.add_argument('--PacBio', help = 'Data from PacBio Sequencing', action='store_true')
    method.add_argument('--PacBioCCS', help = 'Data from PacBio CCS', action='store_true')
    method.add_argument('--ShortRead', help = 'Single end short reads (Legacy support for Illumina)', action='store_true')
    args = cli.parse_args()

    if(args.Nanopore):
        seq_method="map-ont"
    elif(args.PacBio):
        seq_method="map-pb"
    elif(args.PacBioCCS):
        seq_method="asm20"
    elif(args.ShortRead):
        seq_method="sr"

    files = sorted([f for f in glob.glob(args.InputFolder+"/*") if re.search(r'(.*)\.((fastq|fq)(|\.gz))$', f)])
    OutputFolder = os.path.expanduser(args.OutputFolder)
    os.system(f"mkdir -p {OutputFolder}")
    f=open(f"{OutputFolder}/cmd.log", 'w+')

    for j in range(0, len(files)):
        i = files[j]
        base = os.path.splitext(os.path.basename(i))[0]
        #print(base)
        os.system(f"mkdir -p {OutputFolder}")
        minimap2_cmd = f"minimap2 -ax {seq_method} {args.Reference} {i} -t {args.threads} > {OutputFolder}/{base}.sam"
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
