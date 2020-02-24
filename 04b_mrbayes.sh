#!/bin/bash
mkdir -p tree/
input=$(ls tree/mb_input*.nex | tail -n 1)

cd tree

if [ ! -f ${input}.mcmc ]
then
  mb $input
fi
