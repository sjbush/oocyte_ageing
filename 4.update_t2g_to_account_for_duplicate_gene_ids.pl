use strict;
use warnings;

# REQUIREMENTS
my $species = 'Homo_sapiens';
my $in = "/data/home/Stephen/oocyte_atlas/indexes/$species/t2g.txt"; # from 3.make_kb_index.sh
if (!(-e($in))) { print "ERROR: cannot find $in\n"; exit 1; }

# OUTPUT
my $out = "/data/home/Stephen/oocyte_atlas/indexes/$species/t2g2.txt";
open(OUT,'>',$out) or die $!;

# STORE GENE NAMES FOR EACH TRANSCRIPT ID
# we are going to create a list of UNIQUE gene names. To do this, we need to find out whether a gene name is assigned to multiple gene IDs.
my %gene_ids_per_gene_name = ();
open(IN,$in) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $transcript_id = $line[0]; my $gene_id = $line[1]; my $gene_name = $line[2]; my $transcript_name = $line[3]; my $chr = $line[4]; my $start = $line[5]; my $end = $line[6]; my $strand = $line[7];
	  if ($gene_id =~ /^(.*?)\.\d+$/) { $gene_id = $1; }
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  my $loc = "$chr:$start-$end:$strand";
	  $gene_ids_per_gene_name{$gene_name}{$gene_id} = $loc;
	}
close(IN) or die $!;

my %transcript_to_gene_name_lookup = ();
open(IN,$in) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $transcript_id = $line[0]; my $gene_id = $line[1]; my $gene_name = $line[2]; my $transcript_name = $line[3]; my $chr = $line[4]; my $start = $line[5]; my $end = $line[6]; my $strand = $line[7];
	  if ($gene_id =~ /^(.*?)\.\d+$/) { $gene_id = $1; }
	  if ($gene_name eq '') { $gene_name = $gene_id; }	  
	  my $number_of_gene_ids_per_gene_name = scalar keys %{$gene_ids_per_gene_name{$gene_name}};
	  if ($number_of_gene_ids_per_gene_name > 1)
		{ my $new_gene_name = "$gene_name/$gene_id";
		  $gene_name = $new_gene_name;
		}
	  $transcript_to_gene_name_lookup{$transcript_id} = $gene_name;
	}
close(IN) or die $!;

open(IN,$in) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $transcript_id = $line[0];
	  $line[2] = $transcript_to_gene_name_lookup{$transcript_id};
	  my $new_line = join("\t",@line);
	  print OUT "$new_line\n";
	}
close(IN) or die $!;

close(OUT) or die $!;
exit 1;