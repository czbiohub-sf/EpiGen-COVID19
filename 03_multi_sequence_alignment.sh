#!/bin/bash
# conda activate covid19

msadir=$1

inputfile=$(ls $msadir/input_muscle_* | tail -n 1)
msafile=${inputfile/input_/msa_}
if [ ! -f $msafile ]
then
  if [[ $(ls $msadir/msa_muscle*) ]]
  then
    # If there is an existing MSA, then add sequences to existing MSA
    lastest=$(ls $msadir/msa_muscle* | tail -n 1)
    muscle -in $inputfile -out ${msafile}temp_msa.fasta
    muscle -profile -in1 $lastest -in2 ${msafile}temp_msa.fasta -out $msafile
    rm ${msafile}temp_msa.fasta
  else
    muscle -in $inputfile -out $msafile
  fi
fi