#! /bin/bash
set -e
set -u
set -o pipefail

root=`pwd`

cat 7292files.txt | while read -r fastq_file
do
	fastqc /home/RNA-Seq/data_mnt/ftp/LCS7292/rawdata/7292/${fastq_file}
done
