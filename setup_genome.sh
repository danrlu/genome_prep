# This script will download, index, and prep reference genomes.

# Prep N2 reference genome
RELEASE=WS276
PROJECT=PRJNA13758

# Download and prepare FASTA Genome
wget "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/${PROJECT}/sequence/genomic/c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz"
wget "ftp://ftp.wormbase.org/pub/wormbase/releases/${RELEASE}/species/c_elegans/${PROJECT}/c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz"

# If wget not working, can download manually from here:
# ftp://ftp.wormbase.org/pub/wormbase/releases/WS276/species/c_elegans/PRJNA13758/


bwa index c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz
gunzip -c c_elegans.${PROJECT}.${RELEASE}.genomic.fa.gz > c_elegans.${PROJECT}.${RELEASE}.genomic.fa # gzip on Quest doesn't have -k option

samtools faidx c_elegans.${PROJECT}.${RELEASE}.genomic.fa

gatk CreateSequenceDictionary -R c_elegans.${PROJECT}.${RELEASE}.genomic.fa -O c_elegans.${PROJECT}.${RELEASE}.genomic.dict


# to prep the gff3 for bcftools csq
zcat c_elegans.${PROJECT}.${RELEASE}.annotations.gff3.gz | grep -P "\tWormBase\t" > c_elegans.${PROJECT}.${RELEASE}.annotations.WormBase.gff3

Rscript setup_gff_CSQ.R c_elegans.${PROJECT}.${RELEASE}.annotations.WormBase.gff3
