=head

AFTER USAGE, TO VISUALISE THE RELATIONSHIP OF ANY GIVEN GENE'S TPM WITH AGE:

library(tidyverse)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/results/min_depth3000000.min_length70.max_pct_mt10.maturationBoth/MII_data_from_oocyte_atlas_median_downsampled.tpm.filtered.EXACT_AGES.txt',header=T,sep='\t')
df.sub<-subset(df,df$Gene.name == 'GLOD4')
df<-df.sub
cor.test(df$Age,df$TPM,method=c('spearman'))
ggplot(df, aes(x=Age, y=TPM, color=Study)) + geom_point() + ggtitle('GLOD4') + xlab('Age (years)') + ylab('Median TPM\n(across 25 replicates downsampled to 3m reads)')
mean<-mean(df.sub$TPM,na.rm=TRUE)
median<-median(df.sub$TPM,na.rm=TRUE)
range<-max(df.sub$TPM,na.rm=TRUE) - min(df.sub$TPM,na.rm=TRUE)

=cut

use strict;
use warnings;
use Acme::Tools qw(avg median);
use Capture::Tiny qw(capture);

# REQUIREMENTS (1)
my $gene_info = '/data/home/Stephen/oocyte_atlas/prerequisites/Ens111.Homo_sapiens.gene_annotations.txt'; # downloaded from Ensembl BioMart, in this order: Gene stable ID, Gene stable ID version, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene name, Gene type
my $R	   	  = '/data/home/Stephen/miniforge3/envs/R/bin/R';
if (!(-e($R))) 		   { print "ERROR: cannot find $R\n";         exit 1; }
if (!(-e($gene_info))) { print "ERROR: cannot find $gene_info\n"; exit 1; }

# PARAMETERS
my @max_pct_mt  = (qw/10 25 100/);
my $max_seeds   = 25; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my @min_depth   = (qw/3000000 6000000 10000000/); # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $min_length  = 70; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $vivo_vitro  = 'Both'; # 'Vivo'; # 'Vitro';
## the following are parameters for considering whether a gene is significantly differentially expressed
my $nonzero_threshold 	 = 80;
my $abs_rho_threshold 	 = 0.3;
my $min_mean_over_median = 0.8;
my $max_mean_over_median = 1.2;

# STORE GENE ANNOTATIONS
my %annotations = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $desc = $line[2]; my $gene_name = $line[7]; my $gene_type = $line[8];
	  if ($gene_name eq '') { $gene_name = $gene_id; }
	  $annotations{$gene_id}{desc} = $desc;
	  $annotations{$gene_id}{name} = $gene_name;
	  $annotations{$gene_id}{type} = $gene_type;
	}
close(IN) or die $!;

foreach my $max_pct_mt (@max_pct_mt)
	{ foreach my $min_depth (@min_depth)
		{ my $params = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt";

		  # REQUIREMENTS (2)
	  	  my $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/oocyte_atlas_median_downsampled.tpm.filtered.expression"; # from 7b.make_downsampled_tpm_matrix.pl
	  	  if (!(-e($in_file))) { print "ERROR: cannot find $in_file\n"; exit 1; }
		  
		  # OUTPUT
		  my $subdir    = "$params.maturation$vivo_vitro";
		  if (!(-d("/data/home/Stephen/oocyte_atlas/results/$subdir"))) { mkdir "/data/home/Stephen/oocyte_atlas/results/$subdir" or die $!; }
	 	  my $out_file1 = "/data/home/Stephen/oocyte_atlas/results/$subdir/correlations_of_MII_tpm_with_age.median_downsampled.tpm.filtered.EXACT_AGES.txt";
		  my $out_file2 = "/data/home/Stephen/oocyte_atlas/results/$subdir/MII_data_from_oocyte_atlas_median_downsampled.tpm.filtered.EXACT_AGES.txt";
		  my $out_file3 = "/data/home/Stephen/oocyte_atlas/results/$subdir/correlations_of_MII_tpm_with_age.median_downsampled.tpm.filtered.EXACT_AND_AVG_AGES.txt";
		  my $out_file4 = "/data/home/Stephen/oocyte_atlas/results/$subdir/MII_data_from_oocyte_atlas_median_downsampled.tpm.filtered.EXACT_AND_AVG_AGES.txt";
		  open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!;
		  print OUT1 "Gene ID\tGene name\tGene type\tDescription\tNo. of MII oocytes\tNo. of MII oocytes with non-zero expression\t% of MII oocytes with non-zero expression\tAge range (years)\tMin TPM - max TPM\tTPM range\tMean TPM\tMedian TPM\tMean TPM/median TPM\trho\tp-value\tBonferroni corrected p-value\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > $abs_rho_threshold, mean/median TPM > $min_mean_over_median and < $max_mean_over_median, and p-value < 0.05)\n";
		  print OUT2 "Gene ID\tGene name\tSample\tStudy\tAge\tTPM\n";
		  print OUT3 "Gene ID\tGene name\tGene type\tDescription\tNo. of MII oocytes\tNo. of MII oocytes with non-zero expression\t% of MII oocytes with non-zero expression\tAge range (years)\tMin TPM - max TPM\tTPM range\tMean TPM\tMedian TPM\tMean TPM/median TPM\trho\tp-value\tBonferroni corrected p-value\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > $abs_rho_threshold, mean/median TPM > $min_mean_over_median and < $max_mean_over_median, and p-value < 0.05)\n";
		  print OUT4 "Gene ID\tGene name\tSample\tStudy\tAge\tTPM\n";

		  # CORRELATE TPM WITH AGE, BOTH FOR EXACT AGES ONLY AND FOR ALL AGES (i.e. EXACT PLUS AVERAGED AGES)
		  for(my $z=0;$z<=1;$z++)
			{ my %data_per_header = (); my %out_lines = ();
			  my $eof;
			  open(IN,$in_file) or die $!;
			  while(<IN>) { $eof = $.; }
			  close(IN) or die $!;
			  open(IN,$in_file) or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  my $pc = sprintf("%.2f",(($./$eof)*100)); print "$params: $z of 1: $pc%\n";
				  if (($. == 1) or ($. == 2) or ($. == 3) or ($. == 4) or ($. == 5) or ($. == 6))
					{ for(my $x=5;$x<@line;$x++)
						{ my $val = $line[$x];
						  if 	($. == 1) { $data_per_header{$x}{sampl} = $val; }
						  elsif ($. == 2) { $data_per_header{$x}{study} = $val; }
						  elsif ($. == 3) { $data_per_header{$x}{stage} = $val; }
						  elsif ($. == 4) { $data_per_header{$x}{age}   = $val; }
						  elsif ($. == 5) { $data_per_header{$x}{exact} = $val; }
						  elsif ($. == 6) { $data_per_header{$x}{cat}   = $val; }
						}
					}
				  next if ($. <= 6);
				  my $gene_id = $line[0];
				  my $gene_name = $annotations{$gene_id}{name};
				  my $gene_type = $annotations{$gene_id}{type};
				  my $gene_desc = $annotations{$gene_id}{desc};
				  my @tpms_MII = (); my @ages_MII = (); my @samples_MII = (); my @studies_MII = ();
				  for(my $x=5;$x<@line;$x++)
					{ my $tpm   = $line[$x];
					  my $age   = $data_per_header{$x}{age};
					  my $sampl = $data_per_header{$x}{sampl};
					  my $study = $data_per_header{$x}{study};
					  my $stage = $data_per_header{$x}{stage};
					  my $exact = $data_per_header{$x}{exact};
					  if ($z == 0) # i.e. EXACT ages only
						{ if (($age =~ /^\d+$/) and ($exact eq 'yes'))
							{ if ($stage eq 'MII')
								{ push(@tpms_MII,$tpm); push(@ages_MII,$age); push(@samples_MII,$sampl); push(@studies_MII,$study); }
							}
						}
					  elsif ($z == 1) # i.e. ALL ages
						{ if ($age =~ /\d/)
							{ if ($stage eq 'MII')
								{ push(@tpms_MII,$tpm); push(@ages_MII,$age); push(@samples_MII,$sampl); push(@studies_MII,$study); }
							}
						}
					}
				  my $num_MII = @tpms_MII;
				  my $num_MII_non_zero = 0;
				  foreach my $tpm (@tpms_MII) { if ($tpm > 0) { $num_MII_non_zero++; } }
				  my @sorted_ages_MII = sort {$a <=> $b} @ages_MII;
				  my $min_age_MII = $sorted_ages_MII[0]; my $max_age_MII = $sorted_ages_MII[$#sorted_ages_MII]; my $age_range_MII = "$min_age_MII-$max_age_MII";
				  my @sorted_tpms_MII = sort {$a <=> $b} @tpms_MII;
				  my $mean_tpm_MII    = avg(@sorted_tpms_MII);
				  my $median_tpm_MII  = median(@sorted_tpms_MII);
				  my $min_tpm_MII     = $sorted_tpms_MII[0]; 				 $min_tpm_MII = sprintf("%.2f",$min_tpm_MII);
				  my $max_tpm_MII     = $sorted_tpms_MII[$#sorted_tpms_MII]; $max_tpm_MII = sprintf("%.2f",$max_tpm_MII);
				  my $min_max_MII     = "$min_tpm_MII-$max_tpm_MII";
				  my $tpm_range_MII   = $max_tpm_MII-$min_tpm_MII;
				  next if ($max_tpm_MII == 0); # CHECKPOINT: exclude genes which are not detectably expressed in MII oocytes
				  my $mean_over_median_tpm_MII = 'NA';
				  if ($median_tpm_MII > 0) { $mean_over_median_tpm_MII = $mean_tpm_MII/$median_tpm_MII; }
				  my $pct_MII_non_zero = sprintf("%.2f",(($num_MII_non_zero/$num_MII)*100));
				  
				  # DETERMINE SPEARMAN'S RANK CORRELATION BETWEEN AGE AND TPM FOR MII OOCYTES
				  my $ages_MII_list = join(",",@ages_MII); my $tpms_MII_list = join(",",@tpms_MII);
				  my $cmd = "$R -e \"cor.test(c($ages_MII_list),c($tpms_MII_list),method=c('spearman'))\"";
				  my ($stdout, $stderr, $exit) = capture { system($cmd) == 0 or die "system $cmd failed: $?\n" }; # see http://rgr-cyt.org/2017/04/capturing-system-command-output-in-perl
				  if ($exit != 1) { print "ERROR: failed at line $. of $in_file, with $gene_id: $cmd\n"; }
				  my @out = split(/\n/,$stdout);
				  my $p = $out[24]; my $rho = $out[28];
				  if ($p =~ /^.*?p\-value \= (.*?)$/) { $p = $1; }
				  if (($rho =~ /\d/) and (($rho < -1) or ($rho > 1))) { print "ERROR: failed at line $. of $in_file, with $gene_id: output of $cmd was considered to have rho of $rho\n"; exit 1; }
				  my $p_MII = $p; my $rho_MII = $rho;
				  if ( (!(defined($p_MII))) or (!(defined($rho_MII))) )
					{ print "ERROR: unable to determine rho and p-value at line $. of $in_file, with $gene_id\n"; exit 1; }

				  my $adj_p = 'NA';
				  my $is_signif = 'no';
				  if (($pct_MII_non_zero > $nonzero_threshold) and (abs($rho) > $abs_rho_threshold) and ($p_MII < 0.05) and ($mean_over_median_tpm_MII > $min_mean_over_median) and ($mean_over_median_tpm_MII < $max_mean_over_median))
					{ $is_signif = 'yes'; }
				  
				  my $out_line = "$gene_id\t$gene_name\t$gene_type\t$gene_desc\t$num_MII\t$num_MII_non_zero\t$pct_MII_non_zero\t$age_range_MII\t$min_max_MII\t$tpm_range_MII\t$mean_tpm_MII\t$median_tpm_MII\t$mean_over_median_tpm_MII\t$rho_MII\t$p_MII\t$adj_p\t$is_signif";
				  $out_lines{$gene_name} = $out_line;
				  
				  # EXPORT MII DATA IN A FORMAT SUITABLE FOR VISUALISATION WITH R
				  for(my $x=0;$x<@samples_MII;$x++)
					{ my $sample = $samples_MII[$x]; my $study = $studies_MII[$x]; my $age = $ages_MII[$x]; my $tpm = $tpms_MII[$x];
					  if ($z == 0) # i.e. EXACT ages only
						{ print OUT2 "$gene_id\t$gene_name\t$sample\t$study\t$age\t$tpm\n"; }
					  elsif ($z == 1) # i.e. ALL ages
						{ print OUT4 "$gene_id\t$gene_name\t$sample\t$study\t$age\t$tpm\n"; }
					}
				}
			  close(IN) or die $!;

			  my @genes = ();
			  while((my $gene_name,my $irrel)=each(%out_lines))
				{ push(@genes,$gene_name); }
			  my @sorted_genes = sort {"\L$a" cmp "\L$b"} @genes;
			  foreach my $gene (@sorted_genes)
				{ my $out_line = $out_lines{$gene};
				  if ($z == 0) # i.e. EXACT ages only
					{ print OUT1 "$out_line\n"; }
				  elsif ($z == 1) # i.e. ALL ages
					{ print OUT3 "$out_line\n"; }
				}
			}
		  close(OUT1) or die $!; close(OUT2) or die $!; close(OUT3) or die $!; close(OUT4) or die $!;
		}
	}

exit 1;