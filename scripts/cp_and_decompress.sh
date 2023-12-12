#!/bin/bash

INFILE=$1
OUTFILE=$2

if [[ $INFILE == *.gz ]]; then
  pigz -dc $INFILE > $OUTFILE
elif [[ $INFILE == *.zst ]]; then
  zstd -d $INFILE -o $OUTFILE
else
  cp $INFILE $OUTFILE
fi
