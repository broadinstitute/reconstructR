#!/bin/bash

INFILES=${@:1:$#-1}
OUTDIR=${@:$#}

for INFILE in $INFILES; do
    if [[ $INFILE == *.gz ]]; then
        OUTFNAME=$(basename $INFILE .gz)
        pigz -dc $INFILE > "$OUTDIR/$OUTFNAME"
    elif [[ $INFILE == *.zst ]]; then
        OUTFNAME=$(basename $INFILE .zst)
        zstd -d $INFILE -o "$OUTDIR/$OUTFNAME"
    else
        OUTFNAME=$(basename $INFILE)
        cp $INFILE "$OUTDIR/$OUTFNAME"
    fi
done
