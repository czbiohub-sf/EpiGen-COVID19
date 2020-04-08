#!/bin/bash

treedir=$1
msadir=${treedir/tree/msa}
aln=$(ls ${msadir}/msa_mafft_* | tail -n 1)
treefile=${aln/msa_mafft_/iqtree_}
treefile=${treefile/.fasta/.treefile}
treefile=${treefile/msa/tree}
datefile=${treefile/.treefile/_dates.tsv}
outdir=$(basename $treefile)
outdir=${outdir/iqtree_/}
outdir=${outdir/.treefile/}
outdir=${treedir}/treetime_$outdir

echo -e "name\tdate" > $datefile
grep '>' $aln | sed "s/>//" | xargs -I{} grep {} data/sequences/gisaid_metadata.tsv | awk -v OFS='\t' '{print $3,$5}' >> $datefile

mkdir -p $outdir
treetime --aln $aln --tree $treefile --dates $datefile --outdir $outdir --reroot least-squares --clock-filter 3 --tip-slack 3 --confidence --clock-rate 0.0008 --clock-std-dev 0.0005
