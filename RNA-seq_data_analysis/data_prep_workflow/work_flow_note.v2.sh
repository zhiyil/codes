## Merge ##
cat ../01mergedFastq/fq.txt | while read -r f1 f2 f3 f4; do prefix=${f1%%_*}; cat $f1 $f2 $f3 $f4 > ../01mergedFastq/${prefix}.merge.fq.gz; done

## Trim ##
for file in `cat trim_list.txt`; do prefix=${file%%.*}; cutadapt -j 0 -a ${seq} -o ../02trimmedFastq/${prefix}.trimmed.fq.gz ${file} 2> ../02trimmedFastq/${prefix}.trim.log; done

## Mapping ##
mkdir -p ../03hisat2Mapping
cd ../02trimmedFastq
ls *gz > files.txt

for file in `cat files.txt`
do
prefix=${file%%.*}
hisat2 -p8 -x /data/genome/hisat2_indexes/human/GRCh38/grch38_snp_tran/genome_snp_tran --dta --summary-file ../03hisat2Mapping/${prefix}.mapping.summary -U ${file} | samtools sort > ../03hisat2Mapping/${prefix}.sorted.bam 2> ../03hisat2Mapping/${prefix}.mapping.log
done

## I put a sample.list.txt file in .. folder (use prefix=${file$$.*} to generate) ##

# generate a mapping summary report for all samples
cd ../03hisat2Mapping

# form a txt file with all headers being written by echo
echo -e SampleID'\t'Total'\t'Unique'\t'%Unique'\t'Multiple'\t'%Multiple'\t'Overall > total_mapping_report.txt

# grab the mapping numbers from the *.summary files
for ID in `cat ../sample.list.txt`
do
awk -v OFS="\t"  'FNR==1 {total=$1}; FNR==4 {uni=$1; puni=$2}; FNR==5 {mul=$1; pmul=$2}; FNR==6 {all=$1}; END {print total,uni,puni,mul,pmul,all}' ${ID}.mapping.summary >> total_mapping_report.txt
done

#insert sample ID at the beginning of each line
i=1
for ID in `cat ../sample.list.txt`; do ((i++)); sed -i "${i}s/^/${ID}\t/" total_mapping_report.txt; done
# note1: using the -i flag with sed to perform a replacement in-place, overwriting the file
# note2: perform a substitution with s/PATTERN/REPLACEMENT/. In this example PATTERN is ^, the beginnng of the line, and REPLACEMENT is "$ID\t", from the loop variable. The s/// command is within double-quotes so that the shell can expand variables
# note3: use ++ or -- to do increment or decrement by 1, must expand () for it to work
# note4: i was initially set as 1, as we replace ^ with ${ID}\t starting at the 2nd line (the 1st line is header)

sed -i 's/[()]//g' total_mapping_report.txt # this is to get rid of the perentheses around the percentage numbers
# to check the resultant total_mapping_report.txt file, can use this to view:
column -t total_mapping_report.txt

## indexing the bam files ##
for ID in `cat ../sample.list.txt`
do
samtools index ${ID}.sorted.bam 2> ${ID}.index.log
done

## FastQC ##
cd ../02trimmedFastq
for ID in `cat ../sample.list.txt`
do
mkdir -p ../04fastqc/sample_${ID}
fastqc ${ID}.trimmed.fq.gz -o ../04fastqc/sample_${ID} 2> ../04fastqc/sample_${ID}/${ID}.fastqc.log
done

## StringTie assembly ##

# assemble individual gtf #
mkdir -p ../05assbl
cd ../03hisat2Mapping
for ID in `cat ../sample.list.txt`
do
stringtie ${ID}.sorted.bam -p 8 -o ../05assbl/${ID}.gtf -G /data/genome/hisat2_indexes/human/GRCh38/Homo_sapiens.GRCh38.99.gtf 2> ../05assbl/${ID}.gtf.log
done

# merge assemblies #
cd ../05assbl
ls -1 *gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf -G /data/genome/hisat2_indexes/human/GRCh38/Homo_sapiens.GRCh38.99.gtf assembly_GTF_list.txt

#  estimate transcript abundance #
for ID in `cat ../sample.list.txt`
do
mkdir -p ../06trx_abundance/${ID}/
stringtie ${ID}.sorted.bam -p 8 -e -B -o ../06trx_abundance/${ID}/${ID}.abundance.gtf -G ../05assbl/stringtie_merged.gtf -A ../06trx_abundance/${ID}/${ID}.gene_abund.tab 2> ../06trx_abundance/${ID}/${ID}.abundance.gtf.log
done

# Generate counts matrix #
# Prepare a txt file, `prepDE.list.txt` for prepDE.py to run on
mkdir ../07countsMtx
for ID in `cat ../sample.list.txt`; do filePath=`readlink -f ${ID}.abundance.gtf`; echo -e  ${ID}'\t'${filePath} >> ../prepDE.list.txt; done

cd ../07countsMtx
prepDE.py -i prepDE.list.txt
