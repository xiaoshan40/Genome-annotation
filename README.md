Genome assemble and genomic annotation for the grassland caterpillar *Gynaephora qinghaiensis*
======

1.Genome assemble
-------
1.1 quality-filtering 
-------
fastp (v0.20.0) for BGI short-reads
```
fastp -I raw.R1.fq.gz -I raw.R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz -n 0 -f 5 -F 5 -t 5 -T 5
```
Guppy (v3.2.2+9fe0a7811) for Nanopore long-reads
```
guppy_basecaller -i fast5/ -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg -s fastq --device cuda:all:100% --num_callers 16
```

1.2 Genome size and heterozygosity estimation using kmerfreq
-------
```
kmerfreq -k 17 ngs_fastq.lst
```
ngs_fastq.lst is a library file that contains the path of all the input sequence files in fastq format

1.3 Genome assemble and polish
-------
*de novo* genome assemble using NextDenovo (v2.3.113)
```
nextDenovo nextDenovo.run.cfg
```
genome assembly polish using NextPolish (v1.3.014)
```
nextPolish nextPolish.run.cfg
```
remove redundant contigs purge_haplotigs
```
minimap2 -ax map-ont assembly.fasta read.fastq 2> out1.log | samtools view -b -F 0x4 - | samtools sort - -@ 8 -o aligned.bam
```
```
purge_haplotigs readhist aligned.bam
```
```
purge_haplotigs contigcov -i aligned.bam.genecov -l 0 -m 110 -h 200
```
```
purge_haplotigs purge -g assembly.fasta -c coverage_stats.csv -b aligned.bam -a 60
```

1.4 genome assembly completeness evaluation using BUSCO (v5.5.0)
-------
```
busco -i Gqin.genome.fasta -c 8 -o busco -m geno -l insecta_odb10 -offline
```

2.Genome annotation
-----
2.1 Repeat annotation
-------
De novo predictions for long terminal repeat element (LTR) using LTR_Finder (v1.0.7) and LTR_retriever 
```
ltr_finder Gqin.genome.fasta > Gqin.scn
```
```
LTR_retriever -genome Gqin.LTR.fa -infinder Gqin.scn
```
De novo predictions for non-LTR using Repeatmodeler (v2.0.1)
```
BuildDatabase -engine ncbi -name Gqin Gqin.genome.fasta
```
```
RepeatModeler -engine ncbi -database Gqin -pa 80 -ninja_dir pwd_to_NINJA-0.95-cluster_only/NINJA >& run.out
```

repeat sequence annotaion using RepeatMasker v4.0.7
```
cat Gqin.LTR.fa consensi.fa.classified > repeat_ref.fa
```
```
RepeatMasker -a -lib  repeat_ref.fa -pa 80 -gff -dir ./masker Gqin.genome.fasta
```

2.2 Protein coding gene annotation
-------
2.2.1 ab initio prediction using Braker2(v2.1.2)
-------
```
braker.pl --genome=Gqin.genome.masked.fasta —species=Gqin —prot_seq=ref_prot_filter.fasta --threads 50 --gff3
```

2.2.2 homology-based prediction
-------
Exonerate(v2.2.023) 
```
exonerate --model protein2genome --percent 60 ref_prot_filter.fasta Gqin.genome.masked.fasta --showtargetgff yes --showalignment no --score 100 > exonerate_result.txt
```
```
perl Convert_exonerate.pl exonerate_result.txt > Gqin.exonerate.gff
```
GenomeThreader(v1.7.1)
```
gth -species drosophila -genomic Gqin.genome.masked.fasta -protein ref_prot_filter.fasta -intermediate -gff3out -o gth_result.txt
```
```
perl Convert_gth.pl all.gth > Gqin.gth.gff
```
combine results
```
cat Gqin.exonerate.gff Gqin.gth.gff > Gqin.homolog.gff3
```

2.2.3 Transcriptome-based prediction 
-------
hisat (v2.2.1) & samtools (v1.16.1)
```
hisat2-build Gqin.genome.masked.fasta Gqin_genome
```
```
hisat2 --dta -p 20 -x Gqin_genome -1 E1.R1.fastq.gz -2 E1.R2.fastq.gz -S E1.sam
```
```
samtools sort -@ 20 -o E1.bam E1.sam
```
stringtie (v2.1.7)
```
stringtie -p 20 -o E1.gtf -l E1 E1.bam
```
```
stringtie --merge -p 50 -o stringtie_merge.gtf mergelist.txt
```

2.2.4 TransDecoder(v5.5.0) & hmmer (v3.3.2) & blast+ (v2.13.0)
-------
```
pwd_to_TransDecoder/util/gtf_genome_to_cdna_fasta.pl stringtie_merge.gtf Gqin.genome.masked.fasta > transcripts.fasta
```
```
pwd_to_TransDecoder/util/gtf_to_alignment_gff3.pl stringtie_merge.gtf > stringtie_merge.gff3
```
```
TransDecoder.LongOrfs -t transcripts.fasta
```
```
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 60 > blastp.outfmt6
```
```
hmmscan --cpu 40 --domtblout pfam.domtblout Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep
```
```
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
```
```
pwd_to_TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 stringtie_merge.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
```

2.2.5 EVidenceModeler (v1.1.1)
-------
```
pwd_to_EVM/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome Gqin.genome.masked.fasta --gene_predictions gff3/braker.gff3 --protein_alignments gff3/Gqin.homolog.gff3 --transcript_alignments
gff3/transcripts.fasta.transdecoder.genome.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
```
```
pwd_to_EVM/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome Gqin.genome.masked.fasta --weights /data/xiaoshan/Gqin/Annotation/EVM/weights.txt --gene_predictions gff3/braker.gff3 --protein_alignments gff3/Gqin.homolog.gff3 --transcript_alignments gff3/transcripts.fasta.transdecoder.genome.gff3 --output_file_name evm.out --partitions partitions_list.out > commands.list
```
```
pwd_to_EVM/EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
```
```
pwd_to_EVM/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
```
```
pwd_to_EVM/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out —genome Gqin.genome.masked.fasta
```
```
find . -regex ".*evm.out.gff3" -exec cat {} \; > Gqin.evm.gff3
```
```
gffread Gqin.evm.gff3 -g Gqin.genome.softmasked.fasta -x GynQin.OGS.cds.fa -y GynQin.OGS.pep.fa
```

3.Genome functional annotation
-----
3.1 blast to Nr database
-------
```
blastp --db /database/nr -query GynQin.OGS.pep.fa --outfmt 6 -o Gqin.diamond.nr.txt --quiet -e 1e-5 --threads 60 > nr.blastp.outfmt6
```
```
perl get_blast_annotation.pl nr.blastp.outfmt6 uniprot_sprot.fasta > Gqin.nr.Annotation.txt
```

3.2 blast to SwissProt database
-------
```
blastp -query GynQin.OGS.pep.fa -db uniprot_sprot.fasta -outfmt 6 -evalue 1e-5 -num_threads 60 > uniprot.blastp.outfmt6
```
```
perl get_blast_annotation.pl uniprot.blastp.outfmt6 uniprot_sprot.fasta > Gqin.uniprot.Annotation.txt
```

3.3 protein domain searching
-------
```
hmmscan --cpu 40 -E 1e-5 --domtblout pfam.domtblout Pfam-A.hmm GynQin.OGS.pep.fa
```
```
perl get_hmmer_result.pl pfam.domtblout
```

3.4 KEGG（Kyoto Encyclopedia of Genes and Genomes）annotation
-------
Submit annotated amino acid sequences to BlastKOALA（https://www.kegg.jp/blastkoala/)

3.5 GO（Gene Ontology）annotation
-------
Submit annotated amino acid sequences to PANNZER (http://ekhidna2.biocenter.helsinki.fi/sanspanz/)

4.Comparative genomics
-----
4.1 orthofinder
-------
extract the longest sequence for each genes
```
perl get_protein_and_cds.pl GCF_002156985.1_Harm_1.0_protein.faa cds_sequence.txt GCF_002156985.1_Harm_1.0_genomic.gff -prefix HelArm
```

run orthofinder(v2.5.4)
```
orthofinder -f protein
```

4.2 Phylogenetic tree construction
-------

aligned with MAFFT (v7.123b) and trim using trimAl (v1.4.rev22)
```
./mafft_trimal.sh
```
concatenated into a supergene sequence
```
perl combination_fasta.pl > supergene_prot.fa
```
phylogenetic tree construction using iqtree (v2.1.4-beta)
```
iqtree2 -s supergene_prot.fa -T AUTO -m MF
```
```
iqtree2 -s supergene_prot.fa -T AUTO -m  Q.insect+F+R6 -B 1000
```

4.3 estimate divergence times using r8s （v1.81）
-------
```
r8s -b -f r8s_in.txt > r8s_out.txt
```

4.4 Cafe (v5.0)
-------
```
cafe5 -i gene_families.txt -t tree.txt
```


















