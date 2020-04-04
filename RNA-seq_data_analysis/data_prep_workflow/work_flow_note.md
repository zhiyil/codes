## 1. Put all files in the same folder

Go to the toplevel dir of the tree containing the fastq.gz files, run:

`find . -name '*.gz' -exec cp {} /home/ftp/LCS7823/raw_data/ \;`

## 2. Trim adapters in the raw reads

`ls` the data files to `dataFile`, then run the following command:

`while read -r f; do cutadapt -j 0 -a ${seq} -o ../00cleanData/${f}.cut.fq.gz ${f} 2> ../00cleanData/${f}.log; done < dataFile`

## 3. Merge multiple fastq files for the same sample

The files of the raw data are like below:

```{bash}
Elmar-mRNASeq-n12-10_S7_L001_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-4_S8_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-10_S7_L002_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-4_S8_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-10_S7_L003_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-4_S8_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-10_S7_L004_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-4_S8_L004_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-11_S1_L001_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-5_S4_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-11_S1_L002_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-5_S4_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-11_S1_L003_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-5_S4_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-11_S1_L004_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-5_S4_L004_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-12_S10_L001_R1_001.fastq.gz.cut.fq.gz  Elmar-mRNASeq-n12-6_S9_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-12_S10_L002_R1_001.fastq.gz.cut.fq.gz  Elmar-mRNASeq-n12-6_S9_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-12_S10_L003_R1_001.fastq.gz.cut.fq.gz  Elmar-mRNASeq-n12-6_S9_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-12_S10_L004_R1_001.fastq.gz.cut.fq.gz  Elmar-mRNASeq-n12-6_S9_L004_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-1_S5_L001_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-7_S6_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-1_S5_L002_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-7_S6_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-1_S5_L003_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-7_S6_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-1_S5_L004_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-7_S6_L004_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-2_S2_L001_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-8_S11_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-2_S2_L002_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-8_S11_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-2_S2_L003_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-8_S11_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-2_S2_L004_R1_001.fastq.gz.cut.fq.gz    Elmar-mRNASeq-n12-8_S11_L004_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-3_S12_L001_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-9_S3_L001_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-3_S12_L002_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-9_S3_L002_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-3_S12_L003_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-9_S3_L003_R1_001.fastq.gz.cut.fq.gz
Elmar-mRNASeq-n12-3_S12_L004_R1_001.fastq.gz.cut.fq.gz   Elmar-mRNASeq-n12-9_S3_L004_R1_001.fastq.gz.cut.fq.gz
```

Loop over samples 1-12, and merge the fastq files for the same sample number:

`for n in {1..12}; do gunzip -c *-${n}_* > ${n}.cut.fq.gz 2> log; done`

or use `cat` to merge:
`cat ../01mergedFastq/fq.txt | while read -r f1 f2 f3 f4; do prefix=${f1%%_*}; cat $f1 $f2 $f3 $f4 > ../01mergedFastq/${prefix}.merge.fq.gz; done`

*Note: in the above example, each sample has 4 fastq files, therefore, assign the 4 files in $f1, $f2, $f3, and $f4, and use `cat` to merge.*

## 4. Run `fastqc`

```{bash}
for ID in {1..12}
do
mkdir -p ../02fastqc/sample_${ID}
fastqc ${ID}.cut.fq.gz -o ../02fastqc/sample_${ID} 2> ../02fastqc/sample_${ID}/log
done
```

## 5. Mapping the reads by hisat2

For each RNA-Seq sample, map the reads to the genome with HISAT2 using the `--dta` option. **It is highly recommended to use the reference annotation information when mapping the reads, which can be either embedded in the genome index (built with the `--ss` and `--exon` options, see HISAT2 manual)**, or provided separately at run time (using the --known-splicesite-infile option of HISAT2). The SAM output of each HISAT2 run must be **sorted and converted to BAM using samtools** as explained above.

```{bash}
cd ../01mergedFastq

mkdir -p ../03hisat2Mapping

for ID in {1..12}
do
hisat2 -p8 -x /data/genome/hisat2_indexes/human/GRCh38/grch38_snp_tran/genome_snp_tran --dta --summary-file ../03hisat2Mapping/${ID}.mapping.summary -U ${ID}.cut.fq.gz | samtools sort > ../03hisat2Mapping/${ID}.sorted.bam 2> ../03hisat2Mapping/${ID}.log
done
```

*Note: the mouse reference genome is located at `/data/genome/hisat2_indexes/mouse/GRCm38/grcm38_snp_tran/genome_snp_tran`*

## 6. Indexing the BAM files

```{bash}
cd ../03hista2Mapping

for ID in {1..12}
do
samtools index ${ID}.sorted.bam 2> ${ID}.index.log
done
```

## 7. Statistics of the mapping

```{bash}
mkdir -p ../04alnStat

for ID in {1..12}
do
echo sample_${ID}_flagStat >> ..//04alnStat/flag_stat.txt
samtools flagstat ${ID}.bam >> ..//04alnStat/flag_stat.txt
echo sample_${ID}_bam_stat >> ..//04alnStat/bam_stat.txt
bam_stat.py -i ${ID}.bam >> ..//04alnStat/bam_stat.txt
done
```

## 8. Assemble transcriptome by StringTie


### 8.1 First assembe transcripts for each sample using the sorted BAM files by StringTie

For each RNA-Seq sample, run StringTie to assemble the read alignments obtained in the previous step; it is recommended to run StringTie with the `-G` option if the reference annotation is available.

```{bash}
mkdir -p ../05assbl

cd ../03hisat2Mapping

for ID in {1..12}
do
stringtie ${ID}.sorted.bam -p 8 -o ../05assbl/${ID}.gtf -G /data/genome/hisat2_indexes/human/GRCh38/Homo_sapiens.GRCh38.99.gtf 2> ../05assbl/${ID}.gtf.log
done
```

*Note: the mouse gtf file location: `/data/genome/hisat2_indexes/mouse/GRCm38/Mus_musculus.GRCm38.99.gtf`*

### 8.2 Merge assemblies by StringTie

Run `StringTie` with `--merge` in order to generate a non-redundant set of transcripts observed in all the RNA-Seq samples assembled previously. The `stringtie --merge` mode takes as input a list of all the assembled transcripts files (in GTF format) previously obtained for each sample, as well as a reference annotation file (-G option) if available.

```{bash}
cd ../05assbl
ls -1 *gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -o stringtie_merged.gtf -G /data/genome/hisat2_indexes/human/GRCh38/Homo_sapiens.GRCh38.99.gtf assembly_GTF_list.txt
```

*The merged transcriptome is named `stringtie_merged.gtf`*

### 8.3 Using the merged transcriptome to estimate transcript abundance

Use the indiviual BAM files to generate `sample${ID}.abundance.gtf` in each subfolder of `/06trx_abundance/sample${ID}/`. Note: -e and -B flags are used to produce input files for downstream differential analysis by `edgeR` or `DEseq`, or `Ballgown` pipeline.

```{bash}
cd ../03hisat2Mapping
mkdir -p ../06trx_abundance/sample{1..12}/

for ID in {1..12}
do
stringtie ${ID}.sorted.bam -p 8 -e -B -o ../06trx_abundance/sample${ID}/sample${ID}.abundance.gtf -G ../05assbl/stringtie_merged.gtf -A ../06trx_abundance/sample${ID}/sample${ID}.gene_abund.tab 2> ../06trx_abundance/sample${ID}/${ID}.abundance.gtf.log
done
```

### 8.4 Generate counts matrix for `edgeR` analysis by `prepDE.py`

## 8.4.1 Prepare a txt file, `prepDE.list.txt`

```{bash}
mkdir ../07countsMtx

for ID in `cat ../sample.list.txt`; do filePath=`readlink -f ${ID}.abundance.gtf`; echo -e  ${ID}'\t'${filePath} >> ../07countsMtx/prepDE.list.txt; done
```

`cat prepDE.list.txt`

```{bash}
Sample1	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample1/sample1.abundance.gtf
Sample2	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample2/sample2.abundance.gtf
Sample3	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample3/sample3.abundance.gtf
Sample4	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample4/sample4.abundance.gtf
Sample5	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample5/sample5.abundance.gtf
Sample6	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample6/sample6.abundance.gtf
Sample7	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample7/sample7.abundance.gtf
Sample8	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample8/sample8.abundance.gtf
Sample9	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample9/sample9.abundance.gtf
Sample10	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample10/sample10.abundance.gtf
Sample11	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample11/sample11.abundance.gtf
Sample12	/home/RNA-Seq/data_mnt/ftp/LCS7823/06trx_abundance/sample12/sample12.abundance.gtf
```

```{bash}
cd ../07countsMtx
prepDE.py -i prepDE.list.txt
```

`ls`

```{bash}
gene_count_matrix.csv  sample_lst.txt  transcript_count_matrix.csv
```
