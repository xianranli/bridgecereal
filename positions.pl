#!/usr/bin/perl5.30.0 -w
use strict;

open(POSITIONS,"positions.txt");
while(<POSITIONS>){
    chomp;
    
    my $fa_gz=$ARGV[0];
    my ($seqName,$begin,$end) = split(/\s/);
    open(SAMTOOLS,"samtools faidx $fa_gz $seqName:$begin-$end >> query0.fa |");
    while(my $line = <SAMTOOLS>){
    print $line;
    }
    close(SAMTOOLS);
}
close(POSITIONS);

system("bash /mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/script/extract_fasta.sh");

