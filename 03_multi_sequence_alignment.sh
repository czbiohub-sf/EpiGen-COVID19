#!/bin/bash
# conda activate covid19

msadir=$1

inputfile=$(ls $msadir/input_mafft_* | tail -n 1)
msafile=${inputfile/input_/msa_}
if [ ! -f $msafile ]
then
  if [[ $(ls $msadir/msa_mafft*) ]]
  then
    # If there is an existing MSA, then add sequences to existing MSA
    lastest=$(ls $msadir/msa_mafft* | tail -n 1)
    mafft --add $inputfile $latest > $msafile
  else
    mafft $inputfile > $msafile
  fi
fi
