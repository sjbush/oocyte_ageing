use strict;
use warnings;
use LWP::Simple;

# PARAMETERS
my $total_no_of_genes_with_this_go_term = 0;
my $p_or_adj_p = 'p'; # 'adj_p';
my @up_or_down = (qw/up down/); # does expression go up or down with age?
my @max_pct_mt = (qw/10 25 100/);
my $max_seeds  = 25; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my @min_depth  = (qw/3000000 6000000 10000000/); # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $min_length = 70; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $vivo_vitro = 'Both'; # 'Vivo'; # 'Vitro';

# OVERALL OUTPUT FILES
my $out_file1_all = ''; my $out_file2_all = '';
if ($p_or_adj_p eq 'p')
	{ $out_file1_all = "/data/home/Stephen/oocyte_atlas/results/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.EXACT_AGES.txt";
	  $out_file2_all = "/data/home/Stephen/oocyte_atlas/results/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.EXACT_AND_AVG_AGES.txt";
	}
elsif ($p_or_adj_p eq 'adj_p')
	{ $out_file1_all = "/data/home/Stephen/oocyte_atlas/results/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.EXACT_AGES.txt";
	  $out_file2_all = "/data/home/Stephen/oocyte_atlas/results/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.EXACT_AND_AVG_AGES.txt";
	}
open(OUT1_ALL,'>',$out_file1_all) or die $!; open(OUT2_ALL,'>',$out_file2_all) or die $!;
print OUT1_ALL "Genes up- or down-regulated with age\tMinimum read depth (bp)\tMaximum mitochondrial content (as % of total TPM)\tGO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n"; # alternative column 3: "GO ID (depth in GO tree, i.e. maximum distance from a parent term)"
print OUT2_ALL "Genes up- or down-regulated with age\tMinimum read depth (bp)\tMaximum mitochondrial content (as % of total TPM)\tGO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n"; # alternative column 3: "GO ID (depth in GO tree, i.e. maximum distance from a parent term)"

foreach my $up_or_down (@up_or_down)
	{ foreach my $max_pct_mt (@max_pct_mt)
		{ foreach my $min_depth (@min_depth)
			{ my $params = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
			  
			  # REQUIREMENTS
			  my $in_dir1 = ''; my $in_dir2 = '';
			  if ($p_or_adj_p eq 'p')
				{ $in_dir1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AGES"; # from 10.run_topgo_on_genes_whose_TPM_signif_correls_with_age.pl
				  $in_dir2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AND_AVG_AGES"; # from 10.run_topgo_on_genes_whose_TPM_signif_correls_with_age.pl
				}
			  elsif ($p_or_adj_p eq 'adj_p')
				{ $in_dir1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AGES"; # from 10.run_topgo_on_genes_whose_TPM_signif_correls_with_age.pl
				  $in_dir2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AND_AVG_AGES"; # from 10.run_topgo_on_genes_whose_TPM_signif_correls_with_age.pl
				}
			  if (!(-d($in_dir1))) { print "ERROR: cannot find $in_dir1\n"; exit 1; }
			  if (!(-d($in_dir2))) { print "ERROR: cannot find $in_dir2\n"; exit 1; }
			  
			  # OUTPUT FILES PER SAMPLE
			  my $out_file1 = ''; my $out_file2 = '';
			  if ($p_or_adj_p eq 'p')
				{ $out_file1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AGES.txt";
				  $out_file2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AND_AVG_AGES.txt";
				}
			  elsif ($p_or_adj_p eq 'adj_p')
				{ $out_file1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AGES.txt";
				  $out_file2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AND_AVG_AGES.txt";
				}
			  open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
			  print OUT1 "GO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n"; # alternative column 3: "GO ID (depth in GO tree, i.e. maximum distance from a parent term)"
			  print OUT2 "GO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n"; # alternative column 3: "GO ID (depth in GO tree, i.e. maximum distance from a parent term)"
			  
			  for(my $y=0;$y<=1;$y++)
				{ my $in_dir;
				  if 	($y == 0) { $in_dir = $in_dir1; }
				  elsif ($y == 1) { $in_dir = $in_dir2; }
				  for(my $x=0;$x<=2;$x++)
					{ my $in_file; my $go_category;
					  if 	($x == 0) { $in_file = "$in_dir/bp.txt"; $go_category = 'biological process'; }
					  elsif ($x == 1) { $in_file = "$in_dir/mf.txt"; $go_category = 'molecular function'; }
					  elsif ($x == 2) { $in_file = "$in_dir/cc.txt"; $go_category = 'cellular component'; }
					  next if (!(-e($in_file)));
					  # next if ($go_category ne 'biological process');
					  open(GO,$in_file) or die $!;
					  while(<GO>)
						{ next if ($. == 1);
						  my $go_line = $_; chomp($go_line);
						  my @go_line = split(/\t/,$go_line);
						  my $go_id   = $go_line[0];
						  my $go_term = $go_line[1];
						  my $total_no_with_this_term = $go_line[2];
						  my $obs_with_this_term 	  = $go_line[3];
						  my $exp_with_this_term 	  = $go_line[4];
						  my $go_p 					  = $go_line[5];
						  my $pc_with_this_term 	  = sprintf("%.2f",(($obs_with_this_term/$total_no_with_this_term)*100));
						  next if ($exp_with_this_term == 0);
						  my $fold_increase = sprintf("%.2f",($obs_with_this_term/$exp_with_this_term));
						  # next if ($fold_increase < 2); # FILTER: if the observed number of genes with this GO term in the cluster exceeds the expected number by fewer than 2-fold
						  next if ($total_no_with_this_term < $total_no_of_genes_with_this_go_term); # FILTER: GO terms represented by too few a set of genes
						  next if (($go_p =~ /^\d+/) && ($go_p > 0.05)); # FILTER: GO terms insignificantly represented in this cluster				  
							  
=cut
						  # what is the 'depth' of the GO ID, i.e. furthest distance from a parent, in the GO tree?
						  my $go_digits;
						  if ($go_id =~ /^GO\:(\d+)$/) { $go_digits = $1; }
						  my $url = "http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/go.cgi?search=GO%3A"."$go_digits";
						  my $content = get($url);
						  my $data_available = 0;
						  my @distances_to_parent = ();
						  if (defined($content))
							{ my @content = split(/\n/,$content);
							  foreach my $line (@content)
								{ if ($line =~ /\<tr\>\<td align\=right\>\<font color\=\#FF0000\>/)
									{ my @result = split(/\:/,$line);
									  my $partition = $result[0];
									  if ($partition =~ /^.*?(\d+)$/)
										{ my $distance_to_parent = $1;
										  push(@distances_to_parent,$distance_to_parent);
										  $data_available++;
										}
									}
								}
							}
						  my $furthest_distance_to_parent = '.';
						  if ($data_available > 0)
							{ my @sorted_distances_to_parent = sort {$b <=> $a} @distances_to_parent;
							  $furthest_distance_to_parent = $sorted_distances_to_parent[0];
							}
						  my $go_id_with_parent_depth = "$go_id ($furthest_distance_to_parent)";
=cut
						  if ($y == 0)
							{ print OUT1 "$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
							  print OUT1_ALL "$up_or_down\t$min_depth\t$max_pct_mt\t$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
							}
						  elsif ($y == 1)
							{ print OUT2 "$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
							  print OUT2_ALL "$up_or_down\t$min_depth\t$max_pct_mt\t$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
							}
						}
					  close(GO) or die $!;
					}
				}
			  close(OUT1) or die $!; close(OUT2) or die $!;
			}
		}
	}
close(OUT1_ALL) or die $!; close(OUT2_ALL) or die $!;
exit 1;