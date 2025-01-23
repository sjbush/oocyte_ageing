# PURPOSE: how many genes are consistently found to be correlated with age, irrespective of downsampling depth or mitochondrial content?

use strict;
use warnings;

# REQUIREMENTS
my $p_or_adj_p  = 'adj_p'; # 'p';
my @max_pct_mt  = (qw/10 25 100/);
my $max_seeds   = 25; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my @min_depth   = (qw/3000000 6000000 10000000/); # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $min_length  = 70; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $vivo_vitro  = 'Both'; # 'Vivo'; # 'Vitro';
my $header_line = '';
foreach my $max_pct_mt (@max_pct_mt)
	{ foreach my $min_depth (@min_depth)
		{ my $params  = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
		  my $in_file = '';
		  if ($p_or_adj_p eq 'p')
			{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary_without_bonferroni_correction.txt"; } # from 8b.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
		  elsif ($p_or_adj_p eq 'adj_p')
			{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary.txt"; } # from 8b.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
		  next if (!(-e($in_file)));
		  $header_line .= "$params\t";
		}
	}
$header_line =~ s/\t$//;

# OUTPUT
my $out_file0a = "/data/home/Stephen/oocyte_atlas/results/summary_of_meta_dataset_contents.txt";
my $out_file0b = "/data/home/Stephen/oocyte_atlas/results/all_meta_datasets_EXACT_AGES.txt";
my $out_file0c = "/data/home/Stephen/oocyte_atlas/results/all_meta_datasets_EXACT_AND_AVG_AGES.txt";
my $out_file1 = ''; my $out_file2 = ''; my $out_file3 = ''; my $out_file4 = ''; my $out_file5 = ''; my $out_file6 = '';
if ($p_or_adj_p eq 'p')
	{ $out_file1 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_age_without_bonferroni_correction.up.txt';
	  $out_file2 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_age_without_bonferroni_correction.down.txt';
	  $out_file3 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_or_avg_age_without_bonferroni_correction.up.txt';
	  $out_file4 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_or_avg_age_without_bonferroni_correction.down.txt';
	  $out_file5 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_AND_exact_or_avg_age_without_bonferroni_correction.up.txt';
	  $out_file6 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_AND_exact_or_avg_age_without_bonferroni_correction.down.txt';
	}
elsif ($p_or_adj_p eq 'adj_p')
	{ $out_file1 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_age.up.txt';
	  $out_file2 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_age.down.txt';
	  $out_file3 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_or_avg_age.up.txt';
	  $out_file4 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_or_avg_age.down.txt';
	  $out_file5 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_AND_exact_or_avg_age.up.txt';
	  $out_file6 = '/data/home/Stephen/oocyte_atlas/results/genes_whose_TPM_significantly_correlates_with_exact_AND_exact_or_avg_age.down.txt';
	}
open(OUT0A,'>',$out_file0a) or die $!; open(OUT0B,'>',$out_file0b) or die $!; open(OUT0C,'>',$out_file0c) or die $!;
print OUT0A "Minimum read depth (bp)\tMaximum mitochondrial content (as % of total TPM)\tAge of oocyte (i.e. is the meta-dataset restricted to oocytes where only exact ages are recorded, or does it also include those where an averaged age is available?)\tNo. of oocytes\tNo. of genes\t";
print OUT0A "No. of genes where expression level was considered significantly differentially expressed with age WITHOUT correcting for multiple hypothesis testing (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < 0.05)\t";
print OUT0A "% of genes where expression level was considered significantly differentially expressed with age WITHOUT correcting for multiple hypothesis testing\t";
print OUT0A "p-value threshold for Bonferroni correction\t";
print OUT0A "No. of genes where expression level was considered significantly differentially expressed with age AFTER correcting for multiple hypothesis testing (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < dataset-specific threshold, i.e. Bonferroni-corrected p < 0.05)\t";
print OUT0A "% of genes where expression level was considered significantly differentially expressed with age AFTER correcting for multiple hypothesis testing\n";
print OUT0B "Dataset number\tMaximum mitochondrial content (as % of total TPM)\tMinimum sequencing depth (no. of reads)\tGene ID\tGene name\tGene type\tDescription\tNo. of MII oocytes\tNo. of MII oocytes with non-zero expression\t% of MII oocytes with non-zero expression\tAge range (years)\tMin TPM - max TPM\tTPM range\tMean TPM\tMedian TPM\tMean TPM/median TPM\trho\tp-value\tp-value threshold for Bonferroni correction\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < 0.05)\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < dataset-specific threshold, i.e. Bonferroni-corrected p < 0.05)\n";
print OUT0C "Dataset number\tMaximum mitochondrial content (as % of total TPM)\tMinimum sequencing depth (no. of reads)\tGene ID\tGene name\tGene type\tDescription\tNo. of MII oocytes\tNo. of MII oocytes with non-zero expression\t% of MII oocytes with non-zero expression\tAge range (years)\tMin TPM - max TPM\tTPM range\tMean TPM\tMedian TPM\tMean TPM/median TPM\trho\tp-value\tp-value threshold for Bonferroni correction\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < 0.05)\tIs expression level considered significantly differentially expressed with age? (criteria: absolute rho > 0.3, mean/median TPM > 0.8 and < 1.2, and p-value < dataset-specific threshold, i.e. Bonferroni-corrected p < 0.05)\n";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!; open(OUT5,'>',$out_file5) or die $!; open(OUT6,'>',$out_file6) or die $!;
print OUT1 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";
print OUT2 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";
print OUT3 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";
print OUT4 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";
print OUT5 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";
print OUT6 "Gene name\tGene ID\tDescription\tNo. of datasets in which gene significantly correlates with age\t$header_line\n";

# CREATE FILES COMBINING DATA FROM ALL META-DATASETS, AND A SUMMARY FILE DETAILING THE CONTENTS OF EACH
my $dataset_num = 0; my %data = (); my %datasets = ();
foreach my $max_pct_mt (@max_pct_mt)
	{ foreach my $min_depth (@min_depth)
		{ $dataset_num++;
		  my $params = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
		  $datasets{$max_pct_mt}{$min_depth}{dataset_num} = $dataset_num;
		  my @age_status = ("EXACT_AGES","EXACT_AND_AVG_AGES");
		  foreach my $age_status (@age_status)
			{ my $in_file_p = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.median_downsampled.tpm.filtered.$age_status.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
			  my $in_file_q = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age_with_adj_p.median_downsampled.tpm.filtered.$age_status.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
			  if ( (!(-e($in_file_p))) or (!(-e($in_file_q))) ) { print "ERROR: unable to find either $in_file_p or $in_file_q; cannot make a combined file of all meta-datasets\n"; exit 1; }
			  next if ( (!(-e($in_file_p))) or (!(-e($in_file_q))) );
			  my $no_of_genes = 0; my $no_of_oocytes = 0; my $no_signif_p = 0; my $no_signif_q = 0;
			  open(IN,$in_file_p) or die $!;
			  while(<IN>)
				{ next if ($. == 1);
				  my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  $no_of_genes++;
				  my $gene_id = $line[0]; my $gene_name = $line[1]; my $gene_type = $line[2]; my $desc = $line[3]; $no_of_oocytes = $line[4]; my $no_of_oocytes_with_non_zero_exp = $line[5]; my $pc_of_oocytes_with_non_zero_exp = $line[6]; my $age_range = $line[7];	my $min_to_max_tpm = $line[8]; my $tpm_range = $line[9]; my $mean_tpm = $line[10]; my $median_tpm = $line[11]; my $mean_over_median_tpm = $line[12]; my $rho = $line[13]; my $p = $line[14]; my $is_signif = $line[16];
				  if ($is_signif eq 'yes') { $no_signif_p++; }
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_id} = $gene_id;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_type} = $gene_type;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{desc} = $desc;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes} = $no_of_oocytes;				  
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes_with_non_zero_exp} = $no_of_oocytes_with_non_zero_exp;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{pc_of_oocytes_with_non_zero_exp} = $pc_of_oocytes_with_non_zero_exp;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{age_range} = $age_range;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{min_to_max_tpm} = $min_to_max_tpm;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{tpm_range} = $tpm_range;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_tpm} = $mean_tpm;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{median_tpm} = $median_tpm;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_over_median_tpm} = $mean_over_median_tpm;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{rho} = $rho;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{p} = $p;
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_p} = $is_signif;
				}
			  close(IN) or die $!;
			  my $bonferroni_threshold = 'NA';
			  open(IN,$in_file_q) or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  if ($line[16] =~ /^.*? p\-value \< (.*?)\)$/) { $bonferroni_threshold = $1; }
				  next if ($. == 1);
				  my $gene_name = $line[1]; my $is_signif = $line[16];
				  if ($is_signif eq 'yes') { $no_signif_q++; }
				  $data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_q} = $is_signif;
				}
			  close(IN) or die $!;
			  $datasets{$max_pct_mt}{$min_depth}{bonferroni_threshold} = $bonferroni_threshold;
			  my $age = ''; if ($age_status eq 'EXACT_AGES') { $age = 'exact'; } elsif ($age_status eq 'EXACT_AND_AVG_AGES') { $age = 'exact or average'; }
			  my $pc_signif_p = sprintf("%.2f",(($no_signif_p/$no_of_genes)*100));
			  my $pc_signif_q = sprintf("%.2f",(($no_signif_q/$no_of_genes)*100));
			  print OUT0A "$min_depth\t$max_pct_mt\t$age\t$no_of_oocytes\t$no_of_genes\t$no_signif_p\t$pc_signif_p\t$bonferroni_threshold\t$no_signif_q\t$pc_signif_q\n";
			}
		}
	}
close(OUT0A) or die $!;
my @genes = ();
while((my $gene_name,my $irrel)=each(%data))
	{ push(@genes,$gene_name); }
my @sorted_genes = sort {$a cmp $b} @genes;
foreach my $gene_name (@sorted_genes)
	{ foreach my $max_pct_mt (@max_pct_mt)
		{ foreach my $min_depth (@min_depth)
			{ my $bonferroni_threshold = $datasets{$max_pct_mt}{$min_depth}{bonferroni_threshold};
			  my $dataset_num 		   = $datasets{$max_pct_mt}{$min_depth}{dataset_num};
			  my @age_status = ("EXACT_AGES","EXACT_AND_AVG_AGES");
			  foreach my $age_status (@age_status)
				{ if ( ($age_status eq 'EXACT_AGES') and (exists($data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status})) )
					{ print OUT0B "$dataset_num\t$max_pct_mt\t$min_depth\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_id}\t$gene_name\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_type}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{desc}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes_with_non_zero_exp}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{pc_of_oocytes_with_non_zero_exp}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{age_range}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{min_to_max_tpm}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{tpm_range}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_tpm}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{median_tpm}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_over_median_tpm}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{rho}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{p}\t";
					  print OUT0B "$bonferroni_threshold\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_p}\t";
					  print OUT0B "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_q}\n";
					}
				  elsif ( ($age_status eq 'EXACT_AND_AVG_AGES') and (exists($data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status})) )
					{ print OUT0C "$dataset_num\t$max_pct_mt\t$min_depth\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_id}\t$gene_name\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{gene_type}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{desc}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{no_of_oocytes_with_non_zero_exp}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{pc_of_oocytes_with_non_zero_exp}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{age_range}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{min_to_max_tpm}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{tpm_range}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_tpm}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{median_tpm}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{mean_over_median_tpm}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{rho}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{p}\t";
					  print OUT0C "$bonferroni_threshold\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_p}\t";
					  print OUT0C "$data{$gene_name}{$max_pct_mt}{$min_depth}{$age_status}{is_signif_q}\n";
					}
				}
			}
		}
	}
close(OUT0B) or die $!; close(OUT0C) or die $!;

my %annotations = (); my %de_genes_up = (); my %de_genes_down = ();
foreach my $max_pct_mt (@max_pct_mt)
	{ foreach my $min_depth (@min_depth)
		{ my $params  = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
		  my $in_file = '';
		  if ($p_or_adj_p eq 'p')
			{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary_without_bonferroni_correction.txt"; } # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
		  elsif ($p_or_adj_p eq 'adj_p')
			{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary.txt"; } # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
		  next if (!(-e($in_file)));
		  open(IN,$in_file) or die $!;
		  while(<IN>)
			{ next if ($. == 1);
			  my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  my $gene_name = $line[0]; my $gene_id = $line[1]; my $desc = $line[3]; my $rho_exact = $line[4]; my $rho_exact_or_avg = $line[6]; my $correl_with_exact_and_avg_age = $line[8];
			  $annotations{$gene_name}{gene_id} = $gene_id;
			  $annotations{$gene_name}{desc}    = $desc;
			  if ($rho_exact !~ /^NA$/) # if there is a value for $rho_exact, this means gene expression is significantly correlated with exact ages
				{ if 	($rho_exact > 0) { $de_genes_up{$gene_name}{exact}{$params}++;   }
				  elsif ($rho_exact < 0) { $de_genes_down{$gene_name}{exact}{$params}++; }
				}
			  if ($rho_exact_or_avg !~ /^NA$/) # if there is a value for $rho_exact_or_avg, this means gene expression is significantly correlated with either exact *or* average ages
				{ if 	($rho_exact_or_avg > 0) { $de_genes_up{$gene_name}{exact_or_avg}{$params}++;   }
				  elsif ($rho_exact_or_avg < 0) { $de_genes_down{$gene_name}{exact_or_avg}{$params}++; }
				}
			  if ($correl_with_exact_and_avg_age eq 'yes') # if gene expression is significantly correlated with both exact AND 'exact or average' age
				{ if (($rho_exact > 0) and ($rho_exact_or_avg > 0))
					{ $de_genes_up{$gene_name}{both}{$params}++; }
				  elsif (($rho_exact < 0) and ($rho_exact_or_avg < 0))
					{ $de_genes_down{$gene_name}{both}{$params}++; }
				}
			}
		  close(IN) or die $!;
		}
	}

for(my $x=0;$x<=1;$x++)
	{ my %de_genes = ();
	  if 	($x == 0) { %de_genes = %de_genes_up;   }
	  elsif ($x == 1) { %de_genes = %de_genes_down; }
	  my @gene_names = ();
	  while((my $gene_name,my $irrel)=each(%de_genes))
		{ push(@gene_names,$gene_name); }
	  my @sorted_gene_names = sort {$a cmp $b} @gene_names;
	  foreach my $gene_name (@sorted_gene_names)
		{ my $out_line_exact = ''; my $out_line_exact_or_avg = ''; my $out_line_both = '';
		  my $num_of_datasets_exact = 0; my $num_of_datasets_exact_or_avg = 0; my $num_of_datasets_both = 0;
		  foreach my $max_pct_mt (@max_pct_mt)
			{ foreach my $min_depth (@min_depth)
				{ my $params  = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
				  my $in_file = '';
			      if ($p_or_adj_p eq 'p')
					{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary_without_bonferroni_correction.txt"; } # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
				  elsif ($p_or_adj_p eq 'adj_p')
					{ $in_file = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.summary.txt"; } # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
			      next if (!(-e($in_file)));
				  
				  if (exists($de_genes{$gene_name}{exact}{$params}))
					{ $out_line_exact .= "yes\t";
					  $num_of_datasets_exact++;
					}
				  else
					{ $out_line_exact .= "no\t"; }
				  
				  if (exists($de_genes{$gene_name}{exact_or_avg}{$params}))
				  	{ $out_line_exact_or_avg .= "yes\t";
					  $num_of_datasets_exact_or_avg++;
					}
				  else
					{ $out_line_exact_or_avg .= "no\t"; }
				  
				  if (exists($de_genes{$gene_name}{both}{$params}))
				  	{ $out_line_both .= "yes\t";
					  $num_of_datasets_both++;
					}
				  else
					{ $out_line_both .= "no\t"; }
				}
			}
		  $out_line_exact =~ s/\t$//; $out_line_exact_or_avg =~ s/\t$//; $out_line_both =~ s/\t$//;
		  if ($x == 0)
			{ print OUT1 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_exact\t$out_line_exact\n" unless ($num_of_datasets_exact == 0);
			  print OUT3 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_exact_or_avg\t$out_line_exact_or_avg\n" unless ($num_of_datasets_exact_or_avg == 0);
			  print OUT5 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_both\t$out_line_both\n" unless ($num_of_datasets_both == 0);
			}
		  elsif ($x == 1)
			{ print OUT2 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_exact\t$out_line_exact\n" unless ($num_of_datasets_exact == 0);
			  print OUT4 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_exact_or_avg\t$out_line_exact_or_avg\n" unless ($num_of_datasets_exact_or_avg == 0);
			  print OUT6 "$gene_name\t$annotations{$gene_name}{gene_id}\t$annotations{$gene_name}{desc}\t$num_of_datasets_both\t$out_line_both\n" unless ($num_of_datasets_both == 0);
			}
		}
	}
close(OUT1) or die $!; close(OUT2) or die $!; close(OUT3) or die $!; close(OUT4) or die $!; close(OUT5) or die $!; close(OUT6) or die $!;
exit 1;