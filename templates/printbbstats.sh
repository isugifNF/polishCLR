#! /usr/bin/env bash
# Desc: Prints out length distibutions, GC, etc of each assembly, could be added to pipeline at the end

module load bbtools

echo "Assmbly stats of $PWD according to bbtools stats.sh"

echo "=== original primary assembly"
stats.sh /project/ag100pest/Pgos/RawData/3-unzip/all_p_ctg.fasta 
echo "=== Arrow polish 1"
stats.sh PolishCLR_Results/02_ArrowPolish/consensus.fasta
echo "=== Arrow polish 2"
stats.sh PolishCLR_Results/02b_ArrowPolish/consensus_02b.fasta
echo "=== Freebayes polish"
stats.sh PolishCLR_Results/04_FreeBayesPolish/consensus_04.fasta
echo "=== Freebayes polish"
stats.sh PolishCLR_Results/06_FreeBayesPolish/final_polished_assembly.fasta
