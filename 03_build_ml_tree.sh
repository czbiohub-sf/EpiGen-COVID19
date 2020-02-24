#!/bin/bash
mkdir -p tree/
input=$(ls msa/msa_muscle_* | tail -n 1)
msa=$(basename $input)
msa=${msa/msa_muscle_/iqtree_}
treeprefix=${msa/.fasta/}
  
if [ ! -f tree/${treeprefix}* ]
then
  iqtree -nt AUTO -s $input -pre tree/$treeprefix -b 500 -bnni -czb
fi