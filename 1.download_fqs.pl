use strict;
use warnings;

# REQUIREMENTS
my $root      = '/data/home/Stephen/oocyte_atlas';
my $species   = 'Homo_sapiens';
my $metadata  = "$root/prerequisites/metadata.$species.txt"; # manually created
my $ssh_file  = ''; # '/home/s/sbush/.aspera/connect/etc/asperaweb_id_dsa.openssh'; # see https://www.biostars.org/p/93482/
my $fatal     = 0;
#if (!(-e($ssh_file)))  { $fatal++; print "ERROR: cannot find $ssh_file\n";  }
if (!(-e($metadata)))  { $fatal++; print "ERROR: cannot find $metadata\n";  }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;
my $ascp_or_wget = 'wget';

# OUTPUT
my $out_dir = "$root/fq";
if (!(-d($out_dir))) 			{ mkdir  $out_dir 			or die $!; }
if (!(-d("$out_dir/$species"))) { mkdir "$out_dir/$species" or die $!; }
my $out_sh  = "$root/run_fq_downloader.$species.sh";
open(SH,'>',$out_sh) or die $!;
print SH "#!/bin/bash\n";
#print SH "module add aspera/3.9.8\n";

# WHAT SAMPLE IDs ARE WE GOING TO RUN?
my %samples = (); my %library_layouts = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $sample_id = $line[4]; my $run_ids = $line[5]; my $library_layout = $line[7];
	  my @run_ids = split(/\,/,$run_ids);
	  $samples{$sample_id}{layout} = $library_layout;
	  $library_layouts{$sample_id}{$library_layout}++;
	  foreach my $run_id (@run_ids)
		{ my $first_3; my $first_6; my $digits;
		  if ($run_id =~ /^(.{3}).*?$/) { $first_3 = $1; }
		  if ($run_id =~ /^(.{6}).*?$/) { $first_6 = $1; }
		  if ($run_id =~ /^.*?(\d+)$/)  { $digits  = $1; }
		  if ((!(defined($first_3))) or (!(defined($first_6))) or (!(defined($digits)))) { print "ERROR: cannot interpret the run IDs ('$run_ids') for $sample_id\n"; }
		  next if ((!(defined($first_3))) or (!(defined($first_6))) or (!(defined($digits)))); # CHECKPOINT: skip downloading if we are unable to interpret the run accession
		  my $number_of_digits = length($digits);
		  my $ena_url_1 = ''; my $ena_url_2 = '';
		  my $wget_url_1 = ''; my $wget_url_2 = '';
		  if ($number_of_digits == 6)
			{ if ($library_layout eq 'paired')
				{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
				  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
				  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
				  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
					}
			  elsif ($library_layout eq 'single')
				{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id.fastq.gz";
				  $ena_url_2  = 'NA';
				  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id.fastq.gz";
				  $wget_url_2 = 'NA';
				}
			}
		  elsif ($number_of_digits == 7)
			{ if ($digits =~ /^.+?(\d{1})$/)
				{ my $last_digit = $1;
				  if ($library_layout eq 'paired')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
					  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
					  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
					}
				  elsif ($library_layout eq 'single')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id.fastq.gz";
					  $ena_url_2  = 'NA';
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id.fastq.gz";
					  $wget_url_2 = 'NA';
					}
				}
			}
		  elsif ($number_of_digits == 8)
			{ if ($digits =~ /^.+?(\d{2})$/)
				{ my $last_two_digits = $1;
				  if ($library_layout eq 'paired')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
					  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
					  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
					}
				  elsif ($library_layout eq 'single')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id.fastq.gz";
					  $ena_url_2  = 'NA';
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id.fastq.gz";
					  $wget_url_2 = 'NA';
					}
				}
			}
		  elsif ($number_of_digits == 9)
			{ if ($digits =~ /^.+?(\d{3})$/)
				{ my $last_three_digits = $1;
				  if ($library_layout eq 'paired')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
					  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
					  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
					}
				  elsif ($library_layout eq 'single')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id.fastq.gz";
					  $ena_url_2  = 'NA';
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id.fastq.gz";
					  $wget_url_2 = 'NA';
					}
				}
			}
		  push(@{$samples{$sample_id}{fq_urls}},[$ena_url_1,$ena_url_2,$wget_url_1,$wget_url_2]);
		}
	}
close(IN) or die $!;

# FOR EACH SAMPLE...
my @samples = ();
while((my $sample_id,my $irrel)=each(%samples))
	{ push(@samples,$sample_id); }
my @sorted_samples = sort {$a cmp $b} @samples;
my $files_seen = 0; my $files_total = @sorted_samples;
foreach my $sample_id (@sorted_samples)
	{ $files_seen++;
	  print "$files_seen of $files_total\n";
	  my $num_library_layouts = scalar keys %{$library_layouts{$sample_id}};
	  next if ($num_library_layouts != 1); # CHECKPOINT: skip if the sample ID is associated with both paired-end and single-end fqs; these have already been pre-processed in some way
	  my $library_layout = '';
	  while((my $this_layout,my $irrel)=each(%{$library_layouts{$sample_id}}))
		{ $library_layout = $this_layout; }
	  
	  # SKIP IF WE'VE OBTAINED THESE FQs ALREADY
	  next if ( ($library_layout eq 'single') and (-e("$out_dir/$species/$sample_id/$sample_id.fq.gz")) );
	  next if ( ($library_layout eq 'paired') and (-e("$out_dir/$species/$sample_id/$sample_id.1.fq.gz")) and (-e("$out_dir/$species/$sample_id/$sample_id.2.fq.gz")) );

	  # IF WE HAVEN'T OBTAINED THESE FQs ALREADY BUT WE HAVE CREATED THE OUTPUT DIR (PERHAPS IN A PREVIOUS FAILED RUN?), THEN PURGE ITS CONTENTS - WE DON'T TRUST THEM
	  if (-e("$out_dir/$species/$sample_id"))
		{ print SH "rm -r $out_dir/$species/$sample_id\n"; }
	  print SH "mkdir $out_dir/$species/$sample_id\n";
	  print SH "cd $out_dir/$species/$sample_id\n";
	  
	  my $fq1 = "$out_dir/$species/$sample_id/$sample_id.1.fq.gz";
	  my $fq2 = "$out_dir/$species/$sample_id/$sample_id.2.fq.gz";
	  my $fq  = "$out_dir/$species/$sample_id/$sample_id.fq.gz";
	  if (exists($samples{$sample_id}{fq_urls}))
		{ # [FQ ROUTE] DOWNLOAD ALL FQs FOR EACH RUN OF THIS SAMPLE. THE VAST MAJORITY OF SAMPLES IN $metadata WILL TAKE THIS ROUTE.
		  my $file_line1 = ''; my $file_line2 = '';
		  my @fqs_to_delete = ();
		  my @arr = @{$samples{$sample_id}{fq_urls}};
		  for(my $x=0;$x<@arr;$x++)
			{ my $url1 = $arr[$x][0]; my $url2 = $arr[$x][1]; my $wget1 = $arr[$x][2]; my $wget2 = $arr[$x][3];
			  my $fq_1 = $url1; my $fq_2 = $url2;
			  if ($fq_1 =~ /^.+\/(.*?)$/) { $fq_1 = $1; }
			  if ($fq_2 =~ /^.+\/(.*?)$/) { $fq_2 = $1; }
			  if ($ascp_or_wget eq 'ascp')
				{ print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url1 $fq_1\n"; }
			  elsif ($ascp_or_wget eq 'wget')
				{ print SH "wget $wget1 -O $fq_1\n"; }
			  if (($url2 eq 'NA') and ($library_layout eq 'single'))
				{ $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  push(@fqs_to_delete,$fq_1);
				}
			  elsif (($url2 ne 'NA') and ($library_layout eq 'paired'))
				{ if ($ascp_or_wget eq 'ascp')
					{ print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url2 $fq_2\n"; }
				  elsif ($ascp_or_wget eq 'wget')
					{ print SH "wget $wget2 -O $fq_2\n"; }
				  $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  $file_line2 .= "$out_dir/$species/$sample_id/$fq_2 ";
				  push(@fqs_to_delete,$fq_1,$fq_2);
				}
			}
		
		  # [FQ ROUTE] IF THERE HAVE BEEN MULTIPLE SETS OF RUN FQS PER SAMPLE, THEN MERGE THEM ALL TOGETHER INTO ONE FINAL SAMPLE.FQ
		  if ($#arr > 0)
			{ if ($library_layout eq 'paired')
				{ print SH "cat $file_line1 > $fq1\n";
				  print SH "cat $file_line2 > $fq2\n";
				}
			  elsif ($library_layout eq 'single')
				{ print SH "cat $file_line1 > $fq\n";
				}
			  foreach my $fq (@fqs_to_delete)
				{ print SH "rm $fq\n"; }
			}
		  elsif ($#arr == 0) # ELSE JUST RENAME THE SINGLE RUN.FQ TO SAMPLE.FQ
			{ if ($library_layout eq 'paired')
				{ print SH "mv $file_line1 $fq1\n";
				  print SH "mv $file_line2 $fq2\n";
				}
			  elsif ($library_layout eq 'single')
				{ print SH "mv $file_line1 $fq\n";
				}
			}
		}
	}
close(SH) or die $!;
exit 1;