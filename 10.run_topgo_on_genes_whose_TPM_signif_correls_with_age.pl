=head

BEFORE USAGE:

conda activate R

=cut

use strict;
use warnings;

# PARAMETERS
my $p_or_adj_p = 'p'; # 'adj_p';
my @up_or_down = (qw/up down/); # does expression go up or down with age?
my @max_pct_mt = (qw/10 25 100/);
my $max_seeds  = 25; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my @min_depth  = (qw/3000000 6000000 10000000/); # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $min_length = 70; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl and 7.make_downsampled_tpm_matrix.pl
my $vivo_vitro = 'Both'; # 'Vivo'; # 'Vitro';

# REQUIREMENTS (1)
my $mappings = '/data/home/Stephen/oocyte_atlas/prerequisites/go_categories.by_ens_id.map'; # from 9.turn_go_terms_into_map_file.pl
if (!(-e($mappings))) { print "ERROR: cannot find $mappings\n"; exit 1; }

# OUTPUT (1)
my $r_script = "/data/home/Stephen/oocyte_atlas/run_topgo_analysis.$p_or_adj_p.R";
open(OUT_R,'>',$r_script) or die $!;
print OUT_R "library(topGO)\n";

foreach my $up_or_down (@up_or_down)
	{ foreach my $max_pct_mt (@max_pct_mt)
		{ foreach my $min_depth (@min_depth)
			{ my $params = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt.maturation$vivo_vitro";
			  
			  # REQUIREMENTS (2)
			  my $in_file1 = ''; my $in_file2 = '';
			  if ($p_or_adj_p eq 'p')
				{ $in_file1 = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.median_downsampled.tpm.filtered.EXACT_AGES.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
				  $in_file2 = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age.median_downsampled.tpm.filtered.EXACT_AND_AVG_AGES.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
				}
			  elsif ($p_or_adj_p eq 'adj_p')
				{ $in_file1 = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age_with_adj_p.median_downsampled.tpm.filtered.EXACT_AGES.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
				  $in_file2 = "/data/home/Stephen/oocyte_atlas/results/$params/correlations_of_MII_tpm_with_age_with_adj_p.median_downsampled.tpm.filtered.EXACT_AND_AVG_AGES.txt"; # from 8c.compare_spearmans_correls_for_tpm_vs_exact_age_and_avg_age.pl
				}
			  my $fatal    = 0;
			  if (!(-e($in_file1))) { $fatal++; print "ERROR: cannot find $in_file1\n"; }
			  if (!(-e($in_file2))) { $fatal++; print "ERROR: cannot find $in_file2\n"; }
			  if ($fatal > 0) { exit 1; }
			  
			  # OUTPUT (2)
			  my $out_file1 = ''; my $out_file2 = ''; my $out_dir_for_go_script1 = ''; my $out_dir_for_go_script2 = '';
			  if ($p_or_adj_p eq 'p')
				{ $out_file1 			  = "/data/home/Stephen/oocyte_atlas/results/$params/genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AGES.txt";
				  $out_file2 			  = "/data/home/Stephen/oocyte_atlas/results/$params/genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AND_AVG_AGES.txt";
				  $out_dir_for_go_script1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AGES";
				  $out_dir_for_go_script2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age_without_bonferroni_correction.$up_or_down.EXACT_AND_AVG_AGES";
				}
			  elsif ($p_or_adj_p eq 'adj_p')
				{ $out_file1 			  = "/data/home/Stephen/oocyte_atlas/results/$params/genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AGES.txt";
				  $out_file2 			  = "/data/home/Stephen/oocyte_atlas/results/$params/genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AND_AVG_AGES.txt";
				  $out_dir_for_go_script1 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AGES";
				  $out_dir_for_go_script2 = "/data/home/Stephen/oocyte_atlas/results/$params/go_terms_enriched_in_genes_with_exp_signif_correl_with_age.$up_or_down.EXACT_AND_AVG_AGES";
				}
			  if (!(-d($out_dir_for_go_script1))) { mkdir "$out_dir_for_go_script1" or die $!; }
			  if (!(-d($out_dir_for_go_script2))) { mkdir "$out_dir_for_go_script2" or die $!; }
			  open(OUT1,'>',$out_file1) or die $!;
			  open(OUT2,'>',$out_file2) or die $!;
			  
			  for(my $x=0;$x<=1;$x++)
				{ my $in_file; my $out_file; my $out_dir_for_go_script;
				  if 	($x == 0) { $in_file = $in_file1; $out_file = $out_file1; $out_dir_for_go_script = $out_dir_for_go_script1; }
				  elsif ($x == 1) { $in_file = $in_file2; $out_file = $out_file2; $out_dir_for_go_script = $out_dir_for_go_script2; }
				  
				  # STORE THE LIST OF GENES WHICH HAVE EXPRESSION CORRELATED WITH AGE IN MII OOCYTES
				  my $num_genes = 0;
				  open(IN,$in_file) or die $!;
				  while(<IN>)
					{ next if ($. == 1);
					  my $line = $_; chomp($line); $line =~ s/[\r\n]//;
					  my @line = split(/\t/,$line);
					  my $gene_id = $line[0]; my $rho = $line[13]; my $is_signif = $line[16];
					  if ($is_signif eq 'yes')
						{ if ( (($up_or_down eq 'up') and ($rho > 0)) or (($up_or_down eq 'down') and ($rho < 0)) )
							{ if 	($x == 0) { print OUT1 "$gene_id\n"; }
							  elsif ($x == 1) { print OUT2 "$gene_id\n"; }
							  $num_genes++;
							}
						}
					}
				  next if ($num_genes == 0);

				  # CREATE R SCRIPT OF TOPGO CODE SO AS TO RUN A GO TERM ENRICHMENT ANALYSIS ON THE GENE SET
				  if ( (!(-e("$out_dir_for_go_script/mf.txt"))) or (!(-e("$out_dir_for_go_script/bp.txt"))) or (!(-e("$out_dir_for_go_script/cc.txt"))) )
					{ print OUT_R "table<-read.table('$out_file',sep='\\t',header=F)\n";
					  print OUT_R "genes<-table\$V1\n";
					  print OUT_R "geneID2GO <- readMappings(file = file('$mappings'))\n";
					  print OUT_R "geneNames <- names(geneID2GO)\n";
					  print OUT_R "geneList <- factor(as.integer(geneNames %in% genes))\n";
					  print OUT_R "names(geneList) <- geneNames\n";
					  if (!(-e("$out_dir_for_go_script/mf.txt")))
						{ print OUT_R "GOdata <- new('topGOdata', ontology = 'MF', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
						  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
						  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
						  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
						  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/mf.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
						}
					  if (!(-e("$out_dir_for_go_script/bp.txt")))
						{ print OUT_R "GOdata <- new('topGOdata', ontology = 'BP', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
						  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
						  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
						  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
						  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/bp.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
						}
					  if (!(-e("$out_dir_for_go_script/cc.txt")))
						{ print OUT_R "GOdata <- new('topGOdata', ontology = 'CC', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
						  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
						  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
						  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
						  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/cc.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
						}
					}
				}
			  close(OUT1) or die $!; close(OUT2) or die $!;
			}
		}
	}
close(OUT_R) or die $!;
exit 1;