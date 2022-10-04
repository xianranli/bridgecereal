#!/usr/bin/perl5.30.0 -w
use strict;
use lib '/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Programs/BioPerl-1.7.8/lib';
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;


my $pre_dir='/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/';
my $sp='IWGS';
my $gene='TraesCS7D02G524200';
my $target_ch='chr7D'; ### modify this to the corresponding chr
my $ref_gb='IWGSCNEW'; ## 'Ojap' IWGSCNEW;

my $sp_dir = $pre_dir.$sp.'C/';
my $sp_gb_dir = '/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/'; ##$sp_dir.'gbs/';
mkdir $sp_dir.'Candidate_genes/' unless -e $sp_dir.'Candidate_genes/';
my $gene_dir = $sp_dir.'Candidate_genes/'.$gene.'/'; # modifify this as the directory you can write into.
mkdir $gene_dir unless -e $gene_dir;

# my $gene_mrna_file = $gene_dir.$gene.'_ref';

my $Working_dir='/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/'; ## 05/17/2022
my $gene_mrna_file = $Working_dir.$gene.'_ref'; ## 05/17/2022

opendir (SpDir, $sp_gb_dir) || die;
my @genomes = ($ref_gb);
foreach (readdir SpDir) {
	push @genomes, $_ unless $_ =~ /\./ || $_ eq $ref_gb;
	}
my $md = '"m D"';

my $flag = $ARGV[0];

my $p_size=10000; ## for upstream segment based on gene direction
my $a_size=1000; ## for downstream segment

#my %g_p_size = ("Mo17" => 15000,  
#                 "CML247" => 32000, 
#                 "CML333" => 36000, 
#                 "Ki3" => 47000,
#                 "HP301"=> 40000,
#                 "Mo18W" => 49000, 
#                 "CML322" => 40000
#                 );

#my %g_a_size = ('IWGSCNEW' => 10000,
#                'Tx303' => 10000,
#                'CML277' => 4000,
#                'Ky21' => 4500,
#                );
#my $sp_inbreds = ''; #F7_DK105_HZS_EP1
if ($flag == 1) {
	my $gene_syn_file_raw = $gene_dir.$gene.'_Haplotype_syn';
	open (Syn, '>'.$gene_syn_file_raw) || die;
	print Syn "query\tqry_S\tqry_E\tGenome\tch\tsbj_St\tsbj_E\tSize\tSimilarity\n";
	foreach my $g (@genomes) {
		my $g_fas = $sp_gb_dir.$g.'/'.$g.'_'.$target_ch.'.fa';
		
		my $IWGSCNEW_mrna_out = $gene_dir.$gene.'-'.$g.'gb_out_m8';


		# system("blastall -p /usr/local/bin/blastn -i $gene_mrna_file -d $g_fas -o $IWGSCNEW_mrna_out -m 8 -e 1e-20 -F $md");
          system("/usr/local/bin/blastn -db $g_fas -query $gene_mrna_file -out $IWGSCNEW_mrna_out -outfmt 6 -evalue 20");
		
		Extract_syn($IWGSCNEW_mrna_out, $g, $target_ch, \*Syn);
		unlink $IWGSCNEW_mrna_out  unless  $g eq 'IWGSCNEW'
		}
	close Syn;	
}

if ($flag == 2) {
	my $gene_syn_file = $gene_dir.$gene.'_Haplotype_syn';
	my $gene_syn_fa = $gene_dir.$gene.'_Haplotype.fa';
	my $blast_self_o = $gene_dir.$gene.'_Haplotype-Self_out_m8';
	my $blast_self_o2 = $gene_dir.$gene.'_Haplotype-Self_out';
	my $blast_mrna_o = $gene_dir.$gene.'_ref_mRNA-Haplotype_out_m8';
	my $blast_mrna_o2 = $gene_dir.$gene.'_ref_mRNA-Haplotype_out';
	my $region_anno_file = $gene_dir.$gene.'_Haplotype_anno';
	my $seq_gap = $gene_dir.$gene.'_Haplotype_N_Gaps';
	my $var_file = $gene_dir.$gene.'_Variations';



   # my $grass_repeats_fa = $gene_dir.'grasrep.fa'; # Grass Repeats 

     my $grass_repeats_fa = $Working_dir.'grasrep.fa'; ## 05/17/2022

     my $blast_repeats_o = $gene_dir.$gene.'_repMask2'; # Repeats


	#
	open (Gap, '>'.$seq_gap) || die;
	print Gap "Genome\tGAP_Start\tGAP_End\n";
	open (Anno, '>'.$region_anno_file) || die;
	print Anno "Genome\tch\tType\tStart\tEnd\tStrand\trefID\n";
	my ($syn_block_hashref, $syn_strand_hashref) = Parse_Syn_file ($gene_syn_file, $ref_gb);

	my $gene_block_fa = Bio::SeqIO->new(-file => '>'.$gene_syn_fa, -format => 'fasta');

	foreach my $g (keys %$syn_block_hashref) {
#		next unless exists $g_p_size{$g}; ## mask this unless only a subset of genome is needed
		my $g_strand = $$syn_strand_hashref{$g};
		my $g_fas = $sp_gb_dir.$g.'/'.$g.'_'.$target_ch.'.fa';
		next unless -e $g_fas;
		my $g_fas_db = Bio::DB::Fasta->new($g_fas);
		my @chs = keys %{ $$syn_block_hashref{$g} };
		my $ch = $chs[0];
		my @bps;
		push @bps, split /\t/, $_ foreach (@{ $$syn_block_hashref{$g}{$ch} });
		@bps = sort { $a <=> $b} @bps;
		my ($s1, $e1);
#		$p_size = $g_p_size{$g} if $g_p_size && exists $g_p_size{$g};
#		$a_size = $g_a_size{$g} if $g_a_size && exists $g_a_size{$g};
		if ($g_strand > 0) {
			  $s1 = $bps[0]  - $p_size; $e1 = $bps[-1] + $a_size; #### modify the numbers to fit your target
#			  if ($sp_inbreds =~ $g ) { $s1 = $bps[0]  - $p_size; $e1 = $bps[-1] + $a_size + 3000;} ## unmask this line and adjust the number to fit your target
			} else {  
				$s1 = $bps[-1] + $p_size; $e1 = $bps[0]  - $a_size;
#			  if ($sp_inbreds =~ $g ) { $s1 = $bps[-1]  + $p_size; $e1 = $bps[0] - $a_size  - 0000;} ## unmask this line and adjust the number to fit your target
#				unless ($g eq 'IWGSCNEW' || $g eq 'W22' || $g eq 'F7') { $s1 = $bps[-1] + 00000; $e1 = $bps[0] - 60000; } ## unmask this line and adjust the number to fit your target
				}
	
	
		my $fas = $g_fas_db->seq($ch, $s1 => $e1);
		my $gap_hashref = Ns_position_size($fas);
		foreach my $gap (sort { $a <=> $b} keys %$gap_hashref) {
			my $size =  $$gap_hashref{$gap};
			next unless $size > 100;
			my $end = $gap + $size;
			print Gap $g."\t".$gap."\t".$end."\n";
			}
		my $fas_obj = Bio::Seq->new(-display_id => $g, -seq => $fas);
		 $gene_block_fa->write_seq($fas_obj);
			
		my $gb_gff_file = $sp_gb_dir.$g.'/'.$g.'_gb.gff3';
		if (-e $gb_gff_file) {
#			my $target_gff_file = $gene_dir.$g.'_tmp_gff3';
			my $awk_para = "'".'$1=='.$target_ch.'&&$4>='.$s1.'&&$5<='.$e1."'";
#			print "awk $awk_para $gb_gff_file\n" if $g eq 'HZS';
			my @CDSs = `awk $awk_para $gb_gff_file`;
			foreach my $info (@CDSs) {
				my @t = split /\s+/, $info;
				next unless $t[2] eq 'CDS' || $t[2] eq 'gene';
				my $arrow = $t[6] eq '-' ? 1 : 2;
				my ($gID, $x) = ('.', '.'); 
				if ($g eq $ref_gb && $t[2] eq 'gene') { ($gID, $x) = $t[-1] =~ /ID\=(\S+)\;(N)ame/}
				print Anno join "\t", ($g, $t[0], $t[2], $t[3] - $s1, $t[4] - $s1, $arrow, $gID);
				print Anno "\n";
				}
			}
		}
	
	#system("formatdb -i $gene_syn_fa -p f");
	#system("blastall -p /usr/local/bin/blastn -i $gene_syn_fa -d $gene_syn_fa -e 1e-10 -m 8 -F $md -o $blast_self_o");
	#system("blastall -p /usr/local/bin/blastn -i $gene_mrna_file -d $gene_syn_fa -e 1e-10 -m 8 -F $md -o $blast_mrna_o");
	#system("blastall -p /usr/local/bin/blastn -i $gene_mrna_file -d $gene_syn_fa -e 1e-10 -F $md -o $blast_mrna_o2");
	#system("blastall -p /usr/local/bin/blastn -i $gene_syn_fa -d $gene_syn_fa -e 1e-10 -F $md -o $blast_self_o2");

	system("/usr/local/bin/makeblastdb -in $gene_syn_fa -dbtype nucl");

	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_syn_fa -out $blast_self_o -outfmt 6 -evalue 10");
	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_mrna_file -out $blast_mrna_o -outfmt 6 -evalue 10");
	system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_mrna_file -out $blast_mrna_o2 -outfmt 0 -evalue 10");
	# system("/usr/local/bin/blastn -db $gene_syn_fa -query $gene_syn_fa -out $blast_self_o2 -outfmt 0 -evalue 10");
	
	system("/usr/local/bin/makeblastdb -in $grass_repeats_fa -dbtype nucl"); # Repeats
  system("/usr/local/bin/blastn -db $grass_repeats_fa -query $gene_syn_fa -out $blast_repeats_o -outfmt 6 -evalue 10"); #Repeats

	
	Parse_mRNA_haplotype_blast($blast_mrna_o2, $var_file);
	}
##################################################
sub Parse_mRNA_haplotype_blast {
	my ($f, $o) = @_;
	my $blast_io = Bio::SearchIO->new(-format => 'blast', -file => $f);
	my @types = ('snp', 'ins', 'del');
	my%hash;
	while (my $result = $blast_io->next_result) {
		next unless $result->query_name =~ /\_CDS/;
		while (my $hit = $result->next_hit) {
			my $sbjt_id = $hit->accession;
#			next unless $sbjt_id eq 'W22'; 
			while (my $hsp = $hit->next_hsp) {
				next unless $hsp->evalue < 1e-40;
				foreach my $range ( $hsp->seq_inds('sbjct', 'nomatch', 1)) {
					if ($range =~ /(\d+)\-(\d+)/) {
						$hash{$sbjt_id}{'ins'}{$1} = $2 - $1 + 1
						}
						else { $hash{$sbjt_id}{'snp'}{$range} = 1}
					}
				$hash{$sbjt_id}{'del'}{$_} ++ foreach ( $hsp->seq_inds('sbjct', 'gap', 0));
				}
			}
		}
	open (O, '>'.$o) || die;
	print O "Inbred\tType\tBp\tSize\n";
	foreach my $acc (keys %hash ) {
		foreach my $type (@types) {
			next unless exists $hash{$acc}{$type};
			foreach my $bp (keys %{ $hash{$acc}{$type} }) {
				print O $acc."\t".$type."\t".$bp."\t".$hash{$acc}{$type}{$bp}."\n";
				}
			}
		}	
	}

sub Extract_syn {
	my ($f, $g, $ch, $o) = @_;
	open (F, $f) || die;
	my %hash;
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next unless $t[0] =~ /\_mRNA/;
		next unless $t[1] eq $ch; ## if $t[1] =~ /^\d/ && 
		next unless $t[2] > 90;
#		next unless $t[9] < 1e-10
		my $info = $t[0]."\t".$t[6]."\t".$t[7]."\t".$g."\t".$t[1]."\t".$t[8]."\t".$t[9]."\t".$t[3]."\t".$t[2]."\n";
		print $o $info; 
	}
	close F;
#	foreach my $s_s (sort {$a <=> $b} keys %hash) {
#		foreach my $q_s (sort { $a <=> $b} keys %{ $hash{$s_s} }) {
#			print $o $hash{$s_s}{$q_s};
#			}
#		}
 }
sub Parse_Syn_file {
	my ($f, $ref_gb) = @_;
	my (%hash, %hash2);
	open (F, $f) || die;
	my $IWGSCNEW_strand;
	while (<F>) { 
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'query';
		
		my ($genome, $ch, $s_s, $s_e) = @t[3..$#t];
		  $IWGSCNEW_strand =  $s_e - $s_s if $genome eq $ref_gb;
		my $g_strand = $s_e - $s_s;
		if ($g_strand * $IWGSCNEW_strand > 0) {
			$hash2{$genome} =  1;
		} else {
			$hash2{$genome} =  -1;
		}
		 
		push @ { $hash{$genome}{$ch} } , $s_s."\t".$s_e;
		}
	close F;
	return (\%hash, \%hash2);	
	}
sub Ns_position_size {
	my ($fas) = @_;
	my %hash;
	my @t = split //, $fas;
	my $n_0 = 0;
	my $n_tag = 0;
	for (my $i = 0; $i <= $#t; $i ++) {
		if ($t[$i] eq 'N') {
			$n_tag = $i;
			$hash{$n_0} ++
			}
		my $diff = $i - $n_tag;
		$n_0 = $i if $diff > 1;
		}
	return \%hash;	
	}
#my $fas_obj = Bio::SeqIO->new(-file => $gene_syn_fa, -format => 'fasta');
#while (my $seq = $fas_obj->next_seq()) {
#	my $fas = $seq->seq();
#	my $g = $seq->display_id;
#	my $gap_hashref = Ns_position_size($fas);
#	foreach my $gap (sort { $a <=> $b} keys %$gap_hashref) {
#		my $size =  $$gap_hashref{$gap};
#		next unless $size > 100;
#		my $end = $gap + $size;
#		print Gap $g."\t".$gap."\t".$end."\n";
#		}
#	
#	}
