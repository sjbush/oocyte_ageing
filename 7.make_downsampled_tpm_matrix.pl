use strict;
use warnings;
use Acme::Tools qw(median);

# PARAMETERS
my $max_pct_mt = 100; # 10; # 25; # 100;
my $max_seeds  = 25; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl
my $min_depth  = 6000000; # 10000000; # 3000000; # 6000000; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl
my $min_length = 70; # should be identical to that used in 6.run_kb_downsampled_oocyte.pl
my $params 	   = "min_depth$min_depth.min_length$min_length.max_pct_mt$max_pct_mt";
my @studies    = ("Reyes 2017","Ferrero 2019","Leng 2019","Li 2020","Yu 2020","Zhang 2020","Lee 2021","Llonch 2021","Ntostis 2021","Yuan 2021","Takeuchi 2022","Li 2024","This study");
my %studies    = map {$_ => 1} @studies;

# REQUIREMENTS
my $species   = 'Homo_sapiens';
my $root      = '/data/home/Stephen/oocyte_atlas';
my $in_dir    = "$root/kb_downsampled.min_depth$min_depth.min_length$min_length/$species"; # from 6b.run_kb_downsampled_oocyte.pl
my $depths    = "$root/read_length_and_depth_per_sample.$species.txt"; # from 5b.determine_length_and_depth_per_sample.pl
my $metadata  = "$root/prerequisites/metadata.$species.txt"; # manually created
my $gene_info = "$root/prerequisites/Ens111.$species.gene_annotations.txt"; # manually created
my $fatal     = 0;
if (!(-d($in_dir)))    { $fatal++; print "ERROR: cannot find $in_dir\n";    }
if (!(-e($depths)))    { $fatal++; print "ERROR: cannot find $depths\n";    }
if (!(-e($metadata)))  { $fatal++; print "ERROR: cannot find $metadata\n";  }
if (!(-e($gene_info))) { $fatal++; print "ERROR: cannot find $gene_info\n"; }
exit 1 if ($fatal > 0);

# OUTPUT
if (!(-d("$root/results"))) 		{ mkdir "$root/results" 		or die $!; }
if (!(-d("$root/results/$params"))) { mkdir "$root/results/$params" or die $!; }
my $out_file0 = "$root/results/$params/samples_excluded_from_downsampled_oocyte_atlas.txt";
my $out_file1 = "$root/results/$params/pct_mitochondrial_expression_per_downsampled_sample.txt";
my $out_file2 = "$root/results/$params/oocyte_atlas_downsampled.tpm.expression";
my $out_file3 = "$root/results/$params/oocyte_atlas_downsampled.tpm.filtered.expression";
my $out_file4 = "$root/results/$params/oocyte_atlas_median_downsampled.tpm.filtered.expression";
open(OUT0,'>',$out_file0) or die $!; open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!;
print OUT0 "Sample ID\tStudy\tReason(s) for exclusion\n";
print OUT1 "Sample ID\tStudy\tMaturation stage\tAge (exact or average)\tIs this an exact age?\tAge category\tDownsampled replicate\tSummed TPM for mitochondrial genes\t% mitochondrial expression per million transcripts\n";

# STORE GENE ANNOTATIONS
my %annotations = ();
open(IN,$gene_info) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $desc = $line[2]; my $chr = $line[3]; my $gene_start = $line[4]; my $gene_end = $line[5]; my $strand = $line[6]; my $gene_name = $line[7]; my $gene_type = $line[8];
	  my $loc = "$chr:$gene_start-$gene_end:$strand";
	  $annotations{$gene_id}{chr}  = $chr;
	  $annotations{$gene_id}{loc}  = $loc;
	  $annotations{$gene_id}{desc} = $desc;
	  $annotations{$gene_id}{name} = $gene_name;
	  $annotations{$gene_id}{type} = $gene_type;
	}
close(IN) or die $!;

# STORE METADATA PER SAMPLE
my %metadata = (); my %exclusions = (); my %studies_per_sample_id = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $stage = $line[1]; my $sample_id = $line[4]; my $layout = $line[7]; my $exact_age_available = $line[8]; my $exact_age = $line[9]; my $exact_OR_avg_age = $line[10]; my $age_category = $line[11];
	  my $error = 0;
	  if (($exact_age_available eq 'yes') and ($exact_age !~ /^\d+$/))     { $error++; print "ERROR: exact age is said to be available for $sample_id but an age is given of $exact_age\n";   }
	  if (($exact_age_available eq 'no')  and ($exact_age !~ /^unknown$/)) { $error++; print "ERROR: exact age is said to be unavailable for $sample_id but an age is given of $exact_age\n"; }
	  if (($exact_age_available eq 'no')  and ($exact_age !~ /^unknown$/) and ($exact_OR_avg_age !~ /^\d+$/) and ($exact_OR_avg_age !~ /^\d+\.\d+$/))
		{ $error++; print "ERROR: exact age is said to be unavailable for $sample_id but an 'exact or average' age is given of $exact_OR_avg_age\n"; }
	  exit 1 if ($error > 0); # CHECKPOINT: fail if there are metadata-parsing incongruities
	  $studies_per_sample_id{$sample_id} = $study;
	  if ($layout ne 'paired') { $exclusions{$sample_id}{'not paired end'}++; }
	  next if ($layout ne 'paired'); # CHECKPOINT: exclude all datasets which are not paired-end
	  if (!(exists($studies{$study}))) { $exclusions{$sample_id}{"sample from an excluded study, $study"}++; }
	  next if (!(exists($studies{$study}))); # CHECKPOINT: exclude all datasets which are not from an approved study
	  $metadata{$sample_id}{exact_age_available} = $exact_age_available;
	  $metadata{$sample_id}{age_category}	  	 = $age_category;
	  $metadata{$sample_id}{exact_OR_avg_age} 	 = $exact_OR_avg_age;
	  $metadata{$sample_id}{exact_age}        	 = $exact_age;
	  $metadata{$sample_id}{stage}            	 = $stage;
	  $metadata{$sample_id}{study}            	 = $study;
	}
close(IN) or die $!;

# WHAT ARE THE READ LENGTHS AND DEPTHS PER SAMPLE? (NOTE THAT SAMPLES HAD TO MEET MINIMUM REQUIREMENTS OF LENGTH AND DEPTH IN ORDER TO QUALIFY FOR DOWNSAMPLING IN THE FIRST PLACE)
open(IN,$depths) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $sample_id = $line[0]; my $length = $line[3]; my $depth = $line[4];
	  $metadata{$sample_id}{length} = $length;
	  $metadata{$sample_id}{depth}  = $depth;
	}
close(IN) or die $!;

# STORE TPMs PER GENE, PER SAMPLE, FOR APPROVED STUDIES ONLY
opendir(DIR,$in_dir) or die $!;
my @sample_ids_and_seeds = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids_and_seeds = sort {$a cmp $b} @sample_ids_and_seeds;
my %tpms = (); my %mt_tpms = ();
my $samples_seen = 0; my $samples_total = @sorted_sample_ids_and_seeds; $samples_total = $samples_total-2;
my %samples_we_downsampled = ();
foreach my $sample_id_and_seed (@sorted_sample_ids_and_seeds)
	{ next if (($sample_id_and_seed eq '.') or ($sample_id_and_seed eq '..'));
	  my $sample_id; my $seed;
	  if ($sample_id_and_seed =~ /^(.*?)\.(\d+)$/) { $sample_id = $1; $seed = $2; }
	  if ( (!(defined($sample_id))) or (!(defined($seed))) ) { print "ERROR: cannot interpret the filename in $in_dir\n"; exit 1; }
	  next if ( (!(defined($sample_id))) or (!(defined($seed))) ); # CHECKPOINT: we are unable to interpret the filename in $in_dir
	  my $study = $studies_per_sample_id{$sample_id};
	  next if (!(exists($studies{$study})));
	  if (!(exists($metadata{$sample_id}))) { print "ERROR: unable to find metadata for $sample_id\n"; exit 1; }
	  $samples_seen++; print "$samples_seen of $samples_total: $sample_id.$seed\n";
	  my $in_file = "$in_dir/$sample_id.$seed/gene_tpm_matrix.txt";
	  if (!(-e($in_file))) { print "ERROR: cannot find $in_file\n"; }
	  next if (!(-e($in_file)));
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/\"//g;
		  my @line = split(/\t/,$line);
		  my $gene_id = $line[0]; my $tpm = $line[1];
		  if ($gene_id =~ /^(.*?)\.\d+$/) { $gene_id = $1; }
		  $tpms{$gene_id}{$sample_id}{$seed} = $tpm;
		  if ($annotations{$gene_id}{chr} eq 'MT')
			{ $mt_tpms{$sample_id}{$seed} += $tpm; }
		}
	  close(IN) or die $!;
	  $samples_we_downsampled{$sample_id}++;
	}
	
# WHAT SAMPLES - WHICH ARE FROM OTHERWISE APPROVED STUDIES - HAVE WE NOT ACTUALLY DOWNSAMPLED? WE SHOULD CHECK THAT THE REASON WE HAVEN'T DOWNSAMPLED THEM IS BECAUSE THEY FAILED MINIMUM THRESHOLDS FOR DOING SO.
while((my $sample_id,my $irrel)=each(%studies_per_sample_id))
	{ my $study = $studies_per_sample_id{$sample_id};
	  next if (!(exists($studies{$study}))); # CHECKPOINT: we only care about the reason we haven't downsampled samples from approved studies
	  if (!(exists($samples_we_downsampled{$sample_id})))
		{ my $length = $metadata{$sample_id}{length};
		  my $depth  = $metadata{$sample_id}{depth};
		  if ($length < $min_length) { $exclusions{$sample_id}{"could not be downsampled because read length was too short ($length)"}++; }
		  if ($depth  < $min_depth)  { $exclusions{$sample_id}{"could not be downsampled because read depth was too low ($depth)"}++;     }
		  if (($length >= $min_length) and ($depth >= $min_depth))
			{ print "ERROR: we should have downsampled $sample_id as it has long enough reads ($length bp) at high enough depth ($depth)\n"; exit 1; }
		}
	}

# IDENTITY THE MITOCHONDRIAL CONTENT FOR EACH SAMPLE, EXCLUDE THOSE SAMPLES WHERE IT IS TOO HIGH (>25%), AND THEN CREATE HEADER LINES FOR EACH OUTPUT FILE WE'RE GOING TO LATER POPULATE
my %mt_too_high = (); my %sample_ids_seen = ();
my $sample_list = ''; my $study_list = ''; my $stage_list = ''; my $age_list = ''; my $exact_list = ''; my $cat_list = '';
my $sample_list_no_rep = ''; my $study_list_no_rep = ''; my $stage_list_no_rep = ''; my $age_list_no_rep = ''; my $exact_list_no_rep = ''; my $cat_list_no_rep = '';
foreach my $sample_id_and_seed (@sorted_sample_ids_and_seeds)
	{ next if (($sample_id_and_seed eq '.') or ($sample_id_and_seed eq '..'));
	  my $sample_id; my $seed;
	  if ($sample_id_and_seed =~ /^(.*?)\.(\d+)$/) { $sample_id = $1; $seed = $2; }
	  next if ( (!(defined($sample_id))) or (!(defined($seed))) );
	  next if (!(exists($metadata{$sample_id})));
	  next if (!(exists($mt_tpms{$sample_id}{$seed})));
	  my $study = $studies_per_sample_id{$sample_id};
	  next if (!(exists($studies{$study})));
	  my $total_mt = sprintf("%.2f",$mt_tpms{$sample_id}{$seed});
	  my $pct_mt   = sprintf("%.2f",(($total_mt/1000000)*100)); # why this denominator? because TPMs, by definition, sum to 1000000
	  print OUT1 "$sample_id\t$metadata{$sample_id}{study}\t$metadata{$sample_id}{stage}\t$metadata{$sample_id}{exact_OR_avg_age}\t$metadata{$sample_id}{exact_age_available}\t$metadata{$sample_id}{age_category}\t$seed\t$total_mt\t$pct_mt\n";
	  if ($pct_mt > $max_pct_mt)
		{ $mt_too_high{$sample_id}{$seed}++;
		  print "excluding $sample_id replicate $seed (from $metadata{$sample_id}{study}: $metadata{$sample_id}{stage}, $metadata{$sample_id}{exact_OR_avg_age} years) as its % MT is $pct_mt\n";
		}
	  next if ($pct_mt > $max_pct_mt); 
	  $sample_list .= "$metadata{$sample_id}{stage}, $metadata{$sample_id}{exact_OR_avg_age} years ($sample_id), replicate $seed\t";
	  $study_list  .= "$metadata{$sample_id}{study}\t";
	  $stage_list  .= "$metadata{$sample_id}{stage}\t";
	  $age_list    .= "$metadata{$sample_id}{exact_OR_avg_age}\t";
	  $exact_list  .= "$metadata{$sample_id}{exact_age_available}\t";
	  $cat_list    .= "$metadata{$sample_id}{age_category}\t";
	  
	  if (!(exists($sample_ids_seen{$sample_id})))
		{ $sample_list_no_rep .= "$metadata{$sample_id}{stage}, $metadata{$sample_id}{exact_OR_avg_age} years ($sample_id)\t";
		  $study_list_no_rep  .= "$metadata{$sample_id}{study}\t";
		  $stage_list_no_rep  .= "$metadata{$sample_id}{stage}\t";
		  $age_list_no_rep    .= "$metadata{$sample_id}{exact_OR_avg_age}\t";
		  $exact_list_no_rep  .= "$metadata{$sample_id}{exact_age_available}\t";
		  $cat_list_no_rep    .= "$metadata{$sample_id}{age_category}\t";
		  $sample_ids_seen{$sample_id}++;
		}
	}
close(OUT1) or die $!;
$sample_list =~ s/\t$//g; $study_list =~ s/\t$//g; $stage_list =~ s/\t$//g; $age_list =~ s/\t$//g; $exact_list =~ s/\t$//g; $cat_list =~ s/\t$//g;
$sample_list_no_rep =~ s/\t$//g; $study_list_no_rep =~ s/\t$//g; $stage_list_no_rep =~ s/\t$//g; $age_list_no_rep =~ s/\t$//g; $cat_list_no_rep =~ s/\t$//g;
while((my $sample_id,my $irrel)=each(%mt_too_high))
	{ my $num_seeds = scalar keys %{$mt_too_high{$sample_id}};
	  if ($num_seeds == $max_seeds)
		{ $exclusions{$sample_id}{"% mitochondrial reads too high for each of $num_seeds downsampled replicates"}++; }
	}
my $num_included = scalar keys %sample_ids_seen;

print OUT2 "Gene ID\tGene name\tGene type\tLocation\tDescription\t$sample_list\n";
print OUT2 "Study\t\t\t\t\t$study_list\n";
print OUT2 "Maturation stage\t\t\t\t\t$stage_list\n";
print OUT2 "Age\t\t\t\t\t$age_list\n";
print OUT2 "Is age exact?\t\t\t\t\t$exact_list\n";
print OUT2 "Age category\t\t\t\t\t$cat_list\n";
print OUT3 "Gene ID\tGene name\tGene type\tLocation\tDescription\t$sample_list\n";
print OUT3 "Study\t\t\t\t\t$study_list\n";
print OUT3 "Maturation stage\t\t\t\t\t$stage_list\n";
print OUT3 "Age\t\t\t\t\t$age_list\n";
print OUT3 "Is age exact?\t\t\t\t\t$exact_list\n";
print OUT3 "Age category\t\t\t\t\t$cat_list\n";
print OUT4 "Gene ID\tGene name\tGene type\tLocation\tDescription\t$sample_list_no_rep\n";
print OUT4 "Study\t\t\t\t\t$study_list_no_rep\n";
print OUT4 "Maturation stage\t\t\t\t\t$stage_list_no_rep\n";
print OUT4 "Age\t\t\t\t\t$age_list_no_rep\n";
print OUT4 "Is age exact?\t\t\t\t\t$exact_list_no_rep\n";
print OUT4 "Age category\t\t\t\t\t$cat_list_no_rep\n";

my @gene_ids = ();
while((my $gene_id,my $irrel)=each(%tpms))
	{ push(@gene_ids,$gene_id); }
my @sorted_gene_ids = sort {$a cmp $b} @gene_ids;
foreach my $gene_id (@sorted_gene_ids)
	{ my $out_line = '';
	  my %usable_samples = (); my %tpms_per_sample = ();
	  foreach my $sample_id_and_seed (@sorted_sample_ids_and_seeds)
		{ next if (($sample_id_and_seed eq '.') or ($sample_id_and_seed eq '..'));
		  my $sample_id; my $seed;
		  if ($sample_id_and_seed =~ /^(.*?)\.(\d+)$/) { $sample_id = $1; $seed = $2; }
		  next if ( (!(defined($sample_id))) or (!(defined($seed))) );
		  next if (!(exists($metadata{$sample_id})));
		  next if (!(exists($mt_tpms{$sample_id}{$seed})));
		  next if (exists($mt_too_high{$sample_id}{$seed})); # CHECKPOINT: exclude samples were the mitochondrial content is too high
		  my $study = $studies_per_sample_id{$sample_id};
		  next if (!(exists($studies{$study})));
		  my $tpm = 0;
		  if (exists($tpms{$gene_id}{$sample_id}{$seed}))
			{ $tpm = $tpms{$gene_id}{$sample_id}{$seed};
			  push(@{$tpms_per_sample{$sample_id}},$tpm);
			}
		  $out_line .= "$tpm\t";
		  if ($tpm >= 10)
			{ $usable_samples{$sample_id}++; }
		}
	  $out_line =~ s/\t$//g;
	  
	  my $out_line_median = '';
	  my @sample_ids = ();
	  while((my $sample_id,my $irrel)=each(%tpms_per_sample))
		{ push(@sample_ids,$sample_id); }
	  my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
	  foreach my $sample_id (@sorted_sample_ids)
		{ my $median_tpm = median(@{$tpms_per_sample{$sample_id}});
		  $out_line_median .= "$median_tpm\t";
		}	  
	  $out_line_median =~ s/\t$//g;
	  
	  # print $out_line for all genes, and all samples
	  print OUT2 "$gene_id\t$annotations{$gene_id}{name}\t$annotations{$gene_id}{type}\t$annotations{$gene_id}{loc}\t$annotations{$gene_id}{desc}\t$out_line\n";
	  
	  # print $out_line only if the gene is expressed at >= 10 TPM in at least three samples
	  # the purpose of this is facilitate use with BioLayout, as described at https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000859 and https://academic.oup.com/nargab/article/4/1/lqac017/6543592
	  my $num_usable_samples = scalar keys %usable_samples;
	  if ($num_usable_samples >= 3)
		{ print OUT3 "$gene_id\t$annotations{$gene_id}{name}\t$annotations{$gene_id}{type}\t$annotations{$gene_id}{loc}\t$annotations{$gene_id}{desc}\t$out_line\n";
		  print OUT4 "$gene_id\t$annotations{$gene_id}{name}\t$annotations{$gene_id}{type}\t$annotations{$gene_id}{loc}\t$annotations{$gene_id}{desc}\t$out_line_median\n";
		}
	}
close(OUT2) or die $!; close(OUT3) or die $!; close(OUT4) or die $!;

# IF WE EXCLUDED SOME SAMPLES FROM THE FINAL DATASET, WHAT WERE THE REASONS WHY?
my @sample_ids = ();
while((my $sample_id,my $irrel)=each(%exclusions))
	{ push(@sample_ids,$sample_id); }
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
my $num_excluded = 0;
foreach my $sample_id (@sorted_sample_ids)
	{ my $reasons = '';
	  while((my $reason,my $irrel)=each(%{$exclusions{$sample_id}}))
		{ $reasons .= "$reason | "; }
	  $reasons =~ s/\| $//;
	  print OUT0 "$sample_id\t$studies_per_sample_id{$sample_id}\t$reasons\n";
	  $num_excluded++;
	}
close(OUT0) or die $!;

my $total_processed = $num_excluded+$num_included;
print "samples excluded: $num_excluded. samples included: $num_included. total processed: $total_processed\n";

exit 1;