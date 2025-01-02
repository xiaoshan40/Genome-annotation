#!/bin/bash
mkdir group_mafft
mkdir group_trimal
for file in /data/xiaoshan/sound/2024_09/protein/OrthoFinder/Results_Aug25/Single_Copy_Orthologue_Sequences/*.fa
do
        mafft --auto $file >  ${file##*/}.mafft.fasta
        mv  ${file##*/}.mafft.fasta group_mafft
        trimal -in group_mafft/${file##*/}.mafft.fasta -out ${file##*/}.trimal.fasta -automated1
        mv  ${file##*/}.trimal.fasta group_trimal
done
