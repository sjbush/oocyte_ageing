=head

AFTER USAGE:

library(tidyverse)
library(grid)
library(ComplexUpset)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/dataframe_for_making_upset_genes_upregulated_in_MII_oocytes.txt',header=T,sep='\t',check.names=FALSE)
rownames(df)<-df$Gene
df$Gene<-NULL
fig1 <- upset(df, colnames(df), name='Study', width_ratio=0.1, min_size=1, wrap=TRUE, set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))), base_annotations=list('No. of genes'=intersection_size(text=list(vjust=-0.1,hjust=-0.1,angle=45,size=3,col='black'))) ) + ggtitle('A. Genes differentially upregulated in older MII oocytes')

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/dataframe_for_making_upset_genes_downregulated_in_MII_oocytes.txt',header=T,sep='\t',check.names=FALSE)
rownames(df)<-df$Gene
df$Gene<-NULL
fig2 <- upset(df, colnames(df), name='Study', width_ratio=0.1, min_size=1, wrap=TRUE, set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90))), base_annotations=list('No. of genes'=intersection_size(text=list(vjust=-0.1,hjust=-0.1,angle=45,size=3,col='black'))) ) + ggtitle('B. Genes differentially downregulated in older MII oocytes')

png(file="C:/Users/User/Desktop/oocyte_atlas/Documents/figure_2.png",width=10,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=2,ncol=2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig1, vp = define_region(row=1,col=1:2))
print(fig2, vp = define_region(row=2,col=1:2))
dev.off()

=cut

use strict;
use warnings;

# REQUIREMENTS
my $in_file      = 'list_of_genes_up_or_downregulated_in_MII_oocytes.txt';   # manually created, sourced from the tables, text and figures of multiple publictions
my $gene_info_Hs = 'prerequisites/Ens111.Homo_sapiens.gene_annotations.txt'; # manually created from Ensembl BioMart: Gene stable ID, Gene stable ID version, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene name, Gene type
my $gene_info_Mm = 'prerequisites/Ens113.Mus_musculus.gene_annotations.txt'; # manually created from Ensembl BioMart: Gene stable ID, Gene stable ID version, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene name, Gene type
my $orthologues  = 'prerequisites/Ens113.mouse_to_human_orthologues.txt';    # manually created from Ensembl BioMart: Gene stable ID, Gene name, Human gene stable ID, Human gene name, Human homology type, %id. target Human gene identical to query gene, %id. query gene identical to target Human gene, Human Gene-order conservation score, Human Whole-genome alignment coverage, Human orthology confidence [0 low, 1 high]
my $fatal        = 0;
if (!(-e($in_file)))      { $fatal++; print "ERROR: cannot find $in_file\n";      }
if (!(-e($gene_info_Hs))) { $fatal++; print "ERROR: cannot find $gene_info_Hs\n"; }
if (!(-e($gene_info_Mm))) { $fatal++; print "ERROR: cannot find $gene_info_Mm\n"; }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file1 = 'list_of_genes_up_or_downregulated_in_MII_oocytes_reformatted.txt';
my $out_file2 = 'dataframe_for_making_upset_genes_upregulated_in_MII_oocytes.txt';
my $out_file3 = 'dataframe_for_making_upset_genes_downregulated_in_MII_oocytes.txt';
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!;
print OUT1 "Gene name\tEnsembl gene ID\tLocation\tGene type\tDescription\tSource for gene being significantly up-regulated in aging MII oocytes\tNo. of sources\tSource for gene being significantly down-regulated in aging MII oocytes\tNo. of sources\tIs this gene consistently up- or down-regulated? (i.e. no study draws an opposing conclusion)\n";

# STORE HUMAN GENE ANNOTATIONS
my @chrs = (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT/);
my %chrs = map {$_ => 1} @chrs;
my %human_annotations = (); my %human_name_to_id_lookup = ();
open(IN,$gene_info_Hs) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $desc = $line[2]; my $chr = $line[3]; my $gene_start = $line[4]; my $gene_end = $line[5]; my $strand = $line[6]; my $gene_name = $line[7]; my $gene_type = $line[8];
	  next if (!(exists($chrs{$chr}))); # CHECKPOINT: restrict analysis only to genes on the complete chromosomes (i.e. excluding patches)
	  my $loc = "$chr:$gene_start-$gene_end:$strand";
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  $human_annotations{$gene_id}{loc}  = $loc;
	  $human_annotations{$gene_id}{desc} = $desc;
	  $human_annotations{$gene_id}{name} = $gene_name;
	  $human_annotations{$gene_id}{type} = $gene_type;
	  $human_name_to_id_lookup{$gene_name}{$gene_id}++;
	}
close(IN) or die $!;

# STORE MOUSE GENE ANNOTATIONS
@chrs = (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT/);
%chrs = map {$_ => 1} @chrs;
my %mouse_annotations = (); my %mouse_name_to_id_lookup = ();
open(IN,$gene_info_Mm) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $desc = $line[2]; my $chr = $line[3]; my $gene_start = $line[4]; my $gene_end = $line[5]; my $strand = $line[6]; my $gene_name = $line[7]; my $gene_type = $line[8];
	  next if (!(exists($chrs{$chr}))); # CHECKPOINT: restrict analysis only to genes on the complete chromosomes (i.e. excluding patches)
	  my $loc = "$chr:$gene_start-$gene_end:$strand";
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  $mouse_annotations{$gene_id}{loc}  = $loc;
	  $mouse_annotations{$gene_id}{desc} = $desc;
	  $mouse_annotations{$gene_id}{name} = $gene_name;
	  $mouse_annotations{$gene_id}{type} = $gene_type;
	  $mouse_name_to_id_lookup{$gene_name}{$gene_id}++;
	}
close(IN) or die $!;

# STORE HIGH CONFIDENCE MOUSE-TO-HUMAN ONE-TO-ONE ORTHOLOGUES
my %mouse_gene_names_to_human_gene_ids = ();
open(IN,$orthologues) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if ( (!(defined($line[3]))) or (!(defined($line[4]))) );
	  my $mouse_gene_name = $line[1]; my $human_gene_id = $line[2]; my $orthology_type = $line[4]; my $ortho_confidence = $line[9];
	  if (($mouse_gene_name ne '') and ($orthology_type eq 'ortholog_one2one') and ($ortho_confidence == 1))
		{ $mouse_gene_names_to_human_gene_ids{$mouse_gene_name} = $human_gene_id; }
	}
close(IN) or die $!;

# STORE GENES UP- OR DOWNREGULATED IN MII OOCYTES, PER STUDY - BUT FIRST, LET'S MAKE SURE THAT FOR EACH STUDY NO GENE APPEARS ON BOTH THE "UPREGULATED" AND "DOWNREGULATED" GENE LISTS
my %sanity_test = (); my %genes_to_skip = ();
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $val = $line[0]; my $up_or_down = $line[1]; my $source = $line[2]; my $species = $line[3];
	  if ($val =~ /^(.*?)\:.+$/) { $val = $1; }
	  $sanity_test{$source}{$val}{$up_or_down}++;
	}
close(IN) or die $!;
while((my $source,my $irrel)=each(%sanity_test))
	{ while((my $gene_id,my $irrel)=each(%{$sanity_test{$source}}))
		{ my $num_categories = scalar keys %{$sanity_test{$source}{$gene_id}};
		  if ($num_categories != 1)
			{ $genes_to_skip{$source}{$gene_id}++;
			  print "ERROR: in $source, gene $gene_id appears in $num_categories categories of differentially expressed genes; henceforth, skipping this...\n";
			}
		}
	}

# STORE GENES UP- OR DOWNREGULATED IN MII OOCYTES, PER STUDY
my %genes = (); my %associations_per_gene = ();
my $skipped_genes = 0;
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $val = $line[0]; my $up_or_down = $line[1]; my $source = $line[2]; my $species = $line[3];
	  if ($val =~ /^(.*?)\:.+$/) { $val = $1; }
	  
	  next if (exists($genes_to_skip{$source}{$val})); # CHECKPOINT: skip this gene because it appears in both the up- and down-regulated gene lists for a given study, which of course makes no sense (e.g. CEP57 in Grondahl 2010)
	  
	  # if we aren't looking at human data, let's first see if we can change non-human gene names to those of a high-confidence one-to-one orthologue
	  my $human_ortho_id = '';
	  if ( ($species eq 'mouse') and (exists($mouse_gene_names_to_human_gene_ids{$val})) )
		{ $human_ortho_id = $mouse_gene_names_to_human_gene_ids{$val};
		  $val = $human_ortho_id;
		}
	  next if (($species eq 'mouse') and ($human_ortho_id eq '')); # CHECKPOINT: skip this gene because it is non-human, and because we are unable to identify a human orthologue for it
	  
	  # from here on, we are looking up gene names under the assumption they are human
	  if (exists($human_name_to_id_lookup{$val}))
		{ my $num_ids = scalar keys %{$human_name_to_id_lookup{$val}};
		  if ($num_ids > 1)
			{ print "WARNING: $num_ids IDs for $val\n";
			  $skipped_genes++;
			}
		  next if ($num_ids > 1);
		  while((my $id,my $irrel)=each(%{$human_name_to_id_lookup{$val}}))
			{ $genes{$id}{$up_or_down}{$source}++;
			  $associations_per_gene{$up_or_down}{$id}{"$source ($species)"}++;
			}
		}
	  elsif (exists($human_annotations{$val}))
		{ $genes{$val}{$up_or_down}{$source}++;
		  $associations_per_gene{$up_or_down}{$val}{"$source ($species)"}++;
		}
	  else
	   { print "WARNING: unable to find an ID for $val\n";
		 $skipped_genes++;
	   }
	}
close(IN) or die $!;
print "no. of genes skipped because an Ensembl ID could not be associated with them: $skipped_genes\n";

my @genes = ();
while((my $gene_id,my $irrel)=each(%genes))
	{ push(@genes,$gene_id); }
my @sorted_genes = sort {$a cmp $b} @genes;
my $num_up = 0; my $num_down = 0; my $both = 0; my $only_up_supported_by_multiple = 0; my $only_down_supported_by_multiple = 0;
foreach my $gene_id (@sorted_genes)
	{ my $up_sources = ''; my $down_sources = '';
	  my $num_up_sources = 0; my $num_down_sources = 0;
	  if (exists($genes{$gene_id}{'up'}))
		{ my @sources = ();
		  while((my $source,my $irrel)=each(%{$genes{$gene_id}{'up'}}))
			{ push(@sources,$source); }
		  my @sorted_sources = sort {$a cmp $b} @sources;
		  $num_up_sources = @sorted_sources;
		  $up_sources = join(", ",@sorted_sources);
		}
	  if (exists($genes{$gene_id}{'down'}))
		{ my @sources = ();
		  while((my $source,my $irrel)=each(%{$genes{$gene_id}{'down'}}))
			{ push(@sources,$source); }
		  my @sorted_sources = sort {$a cmp $b} @sources;
		  $num_down_sources = @sorted_sources;
		  $down_sources = join(", ",@sorted_sources);
		}
	  if ($up_sources   eq '') { $up_sources   = '.'; }
	  if ($down_sources eq '') { $down_sources = '.'; }
	  my $consistently_up_or_down = '';
	  if ($num_up_sources > 0)
		{ $num_up++;
		  if ($num_up_sources > 1)
			{ $only_up_supported_by_multiple++; }
		  if ($num_down_sources == 0)
			{ $consistently_up_or_down = 'yes'; }
		}
	  if ($num_down_sources > 0)
		{ $num_down++;
		  if ($num_down_sources > 1)
			{ $only_down_supported_by_multiple++; }
		  if ($num_up_sources == 0)
			{ $consistently_up_or_down = 'yes'; }
		}
	  if (($num_up_sources > 0) and ($num_down_sources > 0))
		{ $both++;
		  $consistently_up_or_down = 'no';
		}
	  elsif (($num_up_sources == 0) and ($num_down_sources == 0))
		{ print "ERROR: $gene_id has no evidence of differential gene expression in any study\n"; exit 1; }
	  if ($consistently_up_or_down eq '') { print "ERROR: unable to tell whether $gene_id is consistently up- or downregulated\n"; exit 1; }
	  print OUT1 "$human_annotations{$gene_id}{name}\t$gene_id\t$human_annotations{$gene_id}{loc}\t$human_annotations{$gene_id}{type}\t$human_annotations{$gene_id}{desc}\t$up_sources\t$num_up_sources\t$down_sources\t$num_down_sources\t$consistently_up_or_down\n";
	}
close(OUT1) or die $!;

# OUTPUT A DATA FRAME FOR MAKING AN UPSET PLOT USING THE R PACKAGE ComplexUpset, FOLLOWING THE VIGNETTE AT https://krassowski.github.io/complex-upset/articles/Examples_R.html
my @associations = ("Barone 2020 (human)","Barragán 2017 (human)","Grondahl 2010 (human)","Hamatini 2004 (mouse)","Llonch 2021 (human)","Mishina 2021 (mouse)","Ntostis 2021 (human)","Pan 2008 (mouse)","Reyes 2017 (human)","Steuerwald 2007 (human)","This study (human)","Yuan 2021 (human)","Zhang 2020 (human)");
my $associations_line = join("\t",@associations);
my $num_studies = @associations;
print OUT2 "Gene\t$associations_line\n"; print OUT3 "Gene\t$associations_line\n";
for(my $x=0;$x<=1;$x++)
	{ my $up_or_down = '';
	  if ($x == 0) { $up_or_down = 'up'; } elsif ($x == 1) { $up_or_down = 'down'; }
	  my @gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%{$associations_per_gene{$up_or_down}}))
		{ push(@gene_ids,$gene_id); }
	  my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
	  foreach my $gene_id (@sorted_gene_ids)
		{ my $out_line = '';
		  foreach my $association (@associations)
			{ my $true_or_false = 'FALSE'; # the default
			  if (exists($associations_per_gene{$up_or_down}{$gene_id}{$association}))
				{ $true_or_false = 'TRUE'; }
			  $out_line .= "$true_or_false\t";
			}
		  $out_line =~ s/\t$//;
		  if    ($x == 0) { print OUT2 "$gene_id\t$out_line\n"; }
		  elsif ($x == 1) { print OUT3 "$gene_id\t$out_line\n"; }
		}
	}
close(OUT2) or die $!; close(OUT3) or die $!;

my $num_genes = @sorted_genes;
print "there are $num_genes differentially expressed across $num_studies studies\n";
print "$num_up genes showed evidence of up-regulation with age\n";
print "$num_down genes showed evidence of down-regulation with age\n";
print "$both genes showed evidence both up- and down-regulation with age\n";
my $pc_only_up_supported_by_multiple   = sprintf("%.1f",(($only_up_supported_by_multiple/$num_genes)*100));
my $pc_only_down_supported_by_multiple = sprintf("%.1f",(($only_down_supported_by_multiple/$num_genes)*100));
print "$only_up_supported_by_multiple genes up-regulated with age ($pc_only_up_supported_by_multiple%) are supported by two or more studies\n";
print "$only_down_supported_by_multiple genes down-regulated with age ($pc_only_down_supported_by_multiple%) are supported by two or more studies\n";

exit 1;