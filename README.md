# Genome-annotation
Genome assemble and genomic annotation for the grassland caterpillar Gynaephora qinghaiensis
1. Genome assemble


2. Genome annotation
2.1 Repeat annotation

2.2 Protein coding gene annotation
2.2.1 ab initio prediction
Braker2 (v2.1.2)

2.2.2 homology-based prediction
Exonerate (v2.2.023) 

GenomeThreader (v1.7.1)

2.2.3 Transcriptome-based prediction
hisat2-build -p 28 genome.fa genome
hisat2 --dta -p 20 -x CYMC-genome -1 E1.R1.fastq.gz -2 E1.R2.fastq.gz -S E1.sam
samtools sort -@ 20 -o E1.bam E1.sam
stringtie -p 20 -o E1.gtf -l E1 E1.bam
stringtie --merge -p 50 -o stringtie_merge.gtf mergelist.txt
pwd_to_TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl stringtie_merge.gtf ../CYMC.genome.fasta.masked.NewID > transcripts.fasta
pwd_to_TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl stringtie_merge.gtf > stringtie_merge.gff3










