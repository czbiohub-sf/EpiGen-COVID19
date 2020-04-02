#!/bin/bash

treedir=$1
msadir=${treedir/tree/msa}
mkdir -p $treedir
input=$(ls ${msadir}/msa_mafft_* | tail -n 1)
msa=$(basename $input)
msa=${msa/msa_mafft_/iqtree_}
treeprefix=${msa/.fasta/}
  
if [ ! -f ${treedir}/${treeprefix}* ]
then
  iqtree -nt AUTO -s $input -pre ${treedir}/$treeprefix -b 500 -bnni -czb
fi
