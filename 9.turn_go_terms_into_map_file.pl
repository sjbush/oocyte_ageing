use strict;
use warnings;

# REQUIREMENTS
my $go_cats = '/data/home/Stephen/oocyte_atlas/prerequisites/Ens112.Homo_sapiens.GO_terms.txt'; # from BioMart (Ens112/GRCh38): Gene stable ID, Gene name, GO term evidence code, GO term accession
if (!(-e($go_cats))) { print "ERROR: cannot find $go_cats\n"; exit 1; }

# OUTPUT
my $mappings1 = '/data/home/Stephen/oocyte_atlas/prerequisites/go_categories.by_ens_id.map';
my $mappings2 = '/data/home/Stephen/oocyte_atlas/prerequisites/go_categories.by_name.map';
open(OUT_MAP1,'>',$mappings1) or die $!; open(OUT_MAP2,'>',$mappings2) or die $!;

# CONVERT GO TERM LIST INTO A FORMAT SUITABLE FOR THE R PACKAGE TOPGO
my %go_cats_by_id = (); my %go_cats_by_name = ();
open(IN,$go_cats) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1];
	  next if (!(defined($line[2])));
	  my $evidence_code = $line[2]; my $go_term_accession = $line[3];
	  next if (($evidence_code eq 'NAS') or ($evidence_code eq 'ND')); # non-traceable author statement, or no biological data available
	  $go_cats_by_id{$gene_id}{$go_term_accession}++;
	  $go_cats_by_name{$gene_name}{$go_term_accession}++ unless ($gene_name eq '');	  
	}
close(IN) or die $!;
for(my $x=0;$x<=1;$x++)
	{ my %go_cats = ();
	  if 	($x == 0) { %go_cats = %go_cats_by_id;   }
	  elsif ($x == 1) { %go_cats = %go_cats_by_name; }
	  while((my $gene_id,my $irrel)=each(%go_cats))
		{ my $go_term_line = '';
		  while((my $go_term,my $irrel)=each(%{$go_cats{$gene_id}}))
			{ $go_term_line .= "$go_term, "; }
		  $go_term_line =~ s/\, $//;
		  if    ($x == 0) { print OUT_MAP1 "$gene_id\t$go_term_line\n"; }
		  elsif ($x == 1) { print OUT_MAP2 "$gene_id\t$go_term_line\n"; }
		}
	}
close(OUT_MAP1) or die $!; close(OUT_MAP2) or die $!;

exit 1;