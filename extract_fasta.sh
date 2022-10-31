#!/bin/bash

grep -v "^>" query0.fa | awk 'BEGIN { ORS=""; print ">query\n" } { print }' > query.fasta
rm query0.fa