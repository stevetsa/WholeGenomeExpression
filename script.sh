#!/bin/bash
# Usage: sh script.sh SRR531311 SRR531315 ./RefDir 8 ./outDir > log &

if [ "$#" -ne 5 ]; then
       echo "Not Enough Arguments. Please enter sh script.sh [SRA ID #1] [SRA ID #2] [Reference Genome Dir] [CPUs] [Output Directory]"
       exit 1
fi

SRA1ID="$1"
SRA2ID="$2"
REFDIR=$(readlink -f "$3")
CORES="$4"
OUTDIR=$(readlink -f "$5")

## hg
# For human SRA, uncomment this section and  all sections marked ##hg
#printf "Getting host genome and create BLAST databases\n"
#date
#mkdir $REFDIR
#wget -P $REFDIR ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_genomic.fna.gz
#gunzip $REFDIR/GCF_000001405.37_GRCh38.p11_genomic.fna.gz
#makeblastdb -in $REFDIR/GCF_000001405.37_GRCh38.p11_genomic.fna -dbtype nucl -out $REFDIR/ref
#wget -P $REFDIR ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_genomic.gff.gz
#gunzip $REFDIR/GCF_000001405.37_GRCh38.p11_genomic.gff.gz

##mm
 mouse reference  - default
# 
printf "Getting Reference genome and create BLAST databases\n"
date
mkdir $REFDIR
wget -P $REFDIR `esearch -db assembly -query 'Mus musculus [orgn] AND latest [SB]' | efetch -format docsum | xtract -pattern DocumentSum$
                 awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'` 
gunzip $REFDIR/*.fna.gz
makeblastdb -in $REFDIR/*_genomic.fna -dbtype nucl -out $REFDIR/ref
wget -P $REFDIR ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff$
gunzip $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff.gz

mkdir $OUTDIR
printf "Creating SAM alignment files with magicblast\n"
date
printf "processing %s\n" "$SRA1ID"
date
magicblast -sra $SRA1ID -db $REFDIR/ref -out $OUTDIR/$SRA1ID.sam -num_threads $CORES -no_unaligned
printf "processing %s\n" "$SRA2ID"
date
magicblast -sra $SRA2ID -db $REFDIR/ref -out $OUTDIR/$SRA2ID.sam -num_threads $CORES -no_unaligned


printf "Convert SAM to BAM\n"
date
samtools view -bS $OUTDIR/$SRA1ID.sam | samtools sort -o $OUTDIR/$SRA1ID.bam
samtools view -bS $OUTDIR/$SRA2ID.sam | samtools sort -o $OUTDIR/$SRA2ID.bam
#For older samtools < 1.3
#samtools view -bS $OUTDIR/$SRA1ID.sam | samtools sort - $OUTDIR/$SRA1ID
#samtools view -bS $OUTDIR/$SRA2ID.sam | samtools sort - $OUTDIR/$SRA2ID

printf "Counting\n"
date
##mm
featureCounts -t exon -g gene -a $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff -o $OUTDIR/counts_gene1.txt $OUTDIR/$SRA1ID.bam
featureCounts -t exon -g gene -a $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff -o $OUTDIR/counts_gene2.txt $OUTDIR/$SRA2ID.bam

##hg
#featureCounts -t exon -g gene -a $REFDIR/GCF_000001405.37_GRCh38.p11_genomic.gff -o $OUTDIR/counts_gene1.txt $OUTDIR/$SRA1ID.bam
#featureCounts -t exon -g gene -a $REFDIR/GCF_000001405.37_GRCh38.p11_genomic.gff -o $OUTDIR/counts_gene2.txt $OUTDIR/$SRA2ID.bam

printf "Done\n"
date
