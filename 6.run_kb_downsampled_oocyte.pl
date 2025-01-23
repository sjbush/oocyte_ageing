=head

# NOTE:

The appropriate level to downsample each oocyte library to should be informed by the content of read_length_and_depth_per_sample.Homo_sapiens.txt, produced by 4b.determine_length_and_depth_per_sample.pl. Here we require that each library have a minimum depth of 10m reads, a minimum read length of 75bp, and be paired-end.

# BEFORE USAGE, INSTALL KB INTO CONDA:

conda create -n kb python=3.10
conda activate kb
conda install bioconda::kb-python

# TO SUBMIT THIS TO PBS, USE THE COMMANDS:
# (n.b. why 'source'? see https://community.anaconda.cloud/t/unable-to-activate-environment-prompted-to-run-conda-init-before-conda-activate-but-it-doesnt-work/68677/4)

#!/bin/bash
#PBS -l nodes=3:ppn=4
#PBS -l walltime=72:00:00
#PBS -N kb_oo

source activate base
conda activate kb
cd /data/home/Stephen/oocyte_atlas
/data/home/Stephen/human_adult_SSCs/run_kb_downsampled.Homo_sapiens.oocyte.sh

=cut

use strict;
use warnings;

# PARAMETERS
my $num_procs 	  = 15;
my $max_seeds 	  = 25;
my $downsample_to = 6000000; # 3000000; # 10000000; # also equals the minimum read depth required
my $min_length 	  = 70; # the minimum read length required

# REQUIREMENTS
my $root       = '/data/home/Stephen/oocyte_atlas';
my $species    = 'Homo_sapiens';
my $fq_dir     = "$root/fq/$species"; # from 1.download_fqs.pl and 2.validate_downloaded_fqs.pl
my $t2g        = "$root/indexes/$species/t2g2.txt"; # from 4.update_t2g_to_replace_null_names.pl
my $index      = "$root/indexes/$species/index.idx"; # from 3.make_kb_index.sh
my $depths     = "$root/read_length_and_depth_per_sample.$species.txt"; # from 5b.determine_length_and_depth_per_sample.pl
my $get_matrix = "$root/prerequisites/extract_tpm_per_gene.R"; # manually created
my $metadata   = "$root/prerequisites/metadata.$species.txt"; # manually created
my $se_len     = "$root/determine_read_length_from_fastp.pl"; # manually created
my $rscript    = '/data/home/Stephen/miniforge3/envs/R/bin/Rscript';
my $seqtk      = '/data/home/Stephen/programs/seqtk/seqtk';
my $fastp	   = '/data/home/Stephen/programs/fastp';
my $fatal      = 0;
if (!(-e($get_matrix))) { $fatal++; print "ERROR: cannot find $get_matrix\n"; }
if (!(-e($metadata)))   { $fatal++; print "ERROR: cannot find $metadata\n";   }
if (!(-e($depths)))     { $fatal++; print "ERROR: cannot find $depths\n";     }
if (!(-e($se_len)))     { $fatal++; print "ERROR: cannot find $se_len\n";     }
if (!(-d($fq_dir)))     { $fatal++; print "ERROR: cannot find $fq_dir\n";     }
if (!(-e($rscript))) 	{ $fatal++; print "ERROR: cannot find $rscript\n"; 	  }
if (!(-e($seqtk))) 	    { $fatal++; print "ERROR: cannot find $seqtk\n"; 	  }
if (!(-e($fastp))) 	    { $fatal++; print "ERROR: cannot find $fastp\n"; 	  }
if (!(-e($index))) 	    { $fatal++; print "ERROR: cannot find $index\n"; 	  }
if (!(-e($t2g))) 	    { $fatal++; print "ERROR: cannot find $t2g\n"; 	      }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_dir = "$root/kb_downsampled.min_depth$downsample_to.min_length$min_length";
if (!(-d($out_dir))) 			{ mkdir $out_dir 			or die $!; }
if (!(-d("$out_dir/$species"))) { mkdir "$out_dir/$species" or die $!; }
my $out_sh = "$root/run_kb_downsampled.$species.oocyte.min_depth$downsample_to.min_length$min_length.sh";
open(SH,'>',$out_sh) or die $!;
print SH "#!/bin/bash\n";
print SH "mkdir $out_dir/$species\n" unless (-d("$out_dir/$species"));
print SH "cd $out_dir/$species\n";

# WHAT SAMPLE IDs ARE WE GOING TO RUN?
my %samples = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $sample_id = $line[4]; my $techno = $line[6]; my $layout = $line[7];
	  $samples{$sample_id}{study}  = $study;
	  $samples{$sample_id}{layout} = $layout;
	  $samples{$sample_id}{techno} = $techno;
	}
close(IN) or die $!;

# WHAT ARE THE READ LENGTHS AND DEPTHS PER SAMPLE?
open(IN,$depths) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $sample_id = $line[0]; my $length = $line[3]; my $depth = $line[4];
	  $samples{$sample_id}{length} = $length;
	  $samples{$sample_id}{depth}  = $depth;
	}
close(IN) or die $!;

# FOR EACH SAMPLE...
my @samples = ();
while((my $sample_id,my $irrel)=each(%samples))
	{ push(@samples,$sample_id); }
my @sorted_samples = sort {$a cmp $b} @samples;
my $num = 0; my $total = @sorted_samples;
for(my $x=0;$x<@sorted_samples;$x++)
	{ $num++;
	  print "$num of $total\n";
	  my $sample_id = $sorted_samples[$x];
	  if (!(exists($samples{$sample_id}{length}))) # CHECKPOINT: ensure that each sample has an entry in $depths. If one doesn't, it's because fastp could not be run on it - and that might be due to file-malformatting errors from a corrupted or truncated download, e.g. lines of sequence where the corresponding quality string is of a different length to it
		{ print "ERROR: cannot find data on read lengths and depths for $sample_id - it's not listed in $depths\n"; exit 1; }		  
	  my $study  = $samples{$sample_id}{study};
	  my $techno = $samples{$sample_id}{techno};
	  my $layout = $samples{$sample_id}{layout};
	  my $length = $samples{$sample_id}{length};
	  my $depth  = $samples{$sample_id}{depth};
	  next if ($layout ne 'paired'); # CHECKPOINT: exclude libraries which are not paired-end
	  next if ($length < $min_length); # CHECKPOINT: exclude libraries which do not meet a minimum read length
	  next if ($depth < $downsample_to); # CHECKPOINT: exclude libraries which do not meet a minimum read depth
	  
	  # CHECK WHETHER WE HAVE THE RAW DATA
	  my $fq1 = ''; my $fq2 = ''; my $fq = ''; my $files_available = 0;
	  if ($layout eq 'paired')
		{ $fq1 = "$fq_dir/$sample_id/$sample_id.1.fq.gz";
		  $fq2 = "$fq_dir/$sample_id/$sample_id.2.fq.gz";
		  if ( ( (-e($fq1)) and (-e($fq2)) ) and ( (-s($fq1)) and (-s($fq2)) ) ) { $files_available++; }
		}
	  elsif ($layout eq 'single')
		{ $fq = "$fq_dir/$sample_id/$sample_id.fq.gz";
		  if ((-e($fq)) and (-s($fq))) { $files_available++; }
		}
	  if ($files_available == 0)
		{ print "WARNING: skipping $sample_id as we do not have the fqs\n";
		  if    ( ($layout eq 'paired') and ( (!(-e($fq1))) or (!(-e($fq2))) ) ) { print "*** at least one of the fqs for $sample_id does not exist ***\n"; }
		  elsif ( ($layout eq 'paired') and ( (!(-s($fq1))) or (!(-s($fq2))) ) ) { print "*** at least one of the fqs for $sample_id is incomplete ***\n";  }
		  if    ( ($layout eq 'single') and   (!(-e($fq)))					   ) { print "*** fq for $sample_id does not exist ***\n"; }
		  elsif ( ($layout eq 'single') and   (!(-s($fq)))					   ) { print "*** fq for $sample_id is incomplete ***\n";  }
		}
	  next if ($files_available == 0); # CHECKPOINT: skip if the expected fastq files are unavailable or incomplete
	  
	  # CHECK WHETHER THE FASTQS HAVE BEEN VALIDATED BY FQTOOLS
	  my $passes_validation = 0;
	  my $validation = "$fq_dir/$sample_id/$sample_id.fqtools_validate"; # made by 2.validate_downloaded_fqs.pl
	  if (-e($validation))
		{ open(IN,$validation) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  if ($line =~ /^OK$/) { $passes_validation++; }
			}
		  close(IN) or die $!;
		}
	  if ($passes_validation == 0) { print "WARNING: skipping $sample_id as the fastqs have not been validated\n"; }
	  next if ($passes_validation == 0); # CHECKPOINT: skip if the fastq files have not been validated (using fqtools, as in 1post.validate_downloaded_fqs.pl)
	  
	  # DOWNSAMPLE THE FQs, x TIMES, THEN RUN KB AND EXPORT GENE TPM MATRIX
	  # see advice on commands here: https://github.com/pachterlab/kb_python/issues/167
	  for(my $seed_num=1;$seed_num<=$max_seeds;$seed_num++)
		{ # SKIP IF WE'VE *COMPLETELY* PROCESSED THIS SAMPLE ALREADY. IF WE HAVE, WE'LL HAVE PRODUCED AN OUTPUT DIRECTORY CONTAINING A TPM MATRIX: quant_unfiltered/matrix.abundance.gene.tpm.mtx
		  if ( (-e("$out_dir/$species/$sample_id.$seed_num/quant_unfiltered/matrix.abundance.gene.tpm.mtx")) and (-e("$out_dir/$species/$sample_id.$seed_num/gene_tpm_matrix.txt")) )
			{ print "skipping $sample_id.$seed_num as we've completely processed this sample already\n"; }
		  next if ( (-e("$out_dir/$species/$sample_id.$seed_num/quant_unfiltered/matrix.abundance.gene.tpm.mtx")) and (-e("$out_dir/$species/$sample_id.$seed_num/gene_tpm_matrix.txt")) );
 
	      # HAVE WE *INCOMPLETELY* PROCESS THIS SAMPLE?
		  # if we've already produced an output directory and yet haven't produced the cell x gene TPM matrix (see above checkpoint), then something has gone wrong - kb has not run to completion. In that case, we shall delete this directory and start again.
		  if ( (-d("$out_dir/$species/$sample_id.$seed_num")) and ( (!(-e("$out_dir/$species/$sample_id.$seed_num/quant_unfiltered/matrix.abundance.gene.tpm.mtx"))) ) )
			{ print SH "rm -r $out_dir/$species/$sample_id.$seed_num\n"; }
	  
		  if ( (!(-e("$out_dir/$species/$sample_id.$seed_num/quant_unfiltered/matrix.abundance.gene.tpm.mtx"))) or (!(-e("$out_dir/$species/$sample_id.$seed_num/gene_tpm_matrix.txt"))) )
			{ if (!(-e("$out_dir/$species/$sample_id.$seed_num/quant_unfiltered/matrix.abundance.gene.tpm.mtx")))
				{ if ($layout eq 'paired')
					{ my $seed 	   = int(rand($downsample_to));
					  my $seed_fq1 = "$sample_id.downsampled.1.fq";
					  my $seed_fq2 = "$sample_id.downsampled.2.fq";
					  
					  # downsample fqs
					  print SH "$seqtk sample -s $seed $fq1 $downsample_to > $seed_fq1\n";
					  print SH "$seqtk sample -s $seed $fq2 $downsample_to > $seed_fq2\n";
					  
					  # quantify expression
					  print SH "kb count -i $index -g $t2g -x $techno -o $sample_id.$seed_num --parity $layout --tcc -t $num_procs $seed_fq1 $seed_fq2\n";
					  
					  # delete redundant downsampled fqs
					  print SH "rm $seed_fq1 $seed_fq2\n";
					}
				  elsif ($layout eq 'single')
					{ # downsample fqs
					  my $seed = int(rand($downsample_to));
					  my $seed_fq = "$sample_id.downsampled.fq";
					  print SH "$seqtk sample -s $seed $fq $downsample_to > $seed_fq\n"; 
					  
					  # note that if --parity is "single", then you need to provide --fragment-l and --fragment-s parameters otherwise the pipeline will run with an assumption that all transcripts have the same length (which is obviously untrue): WARNING [main] `--fragment-l` and `--fragment-s` not provided. Assuming all transcripts have the exact same length.
					  my $clean_fq   = "$fq_dir/$sample_id.$seed_num/$sample_id.cleaned.fq.gz";
					  my $fastp_json = "$fq_dir/$sample_id.$seed_num/$sample_id.json";
					  my $fastp_html = "$fq_dir/$sample_id.$seed_num/$sample_id.html";
				  
					  # that in mind, we need to do something a little different with single-end data
					  # for paired-end samples, Kallisto estimates the fragment length from the reads so does not need to be user-specified
					  # for single-end samples, fragment length cannot be empirically derived from read mapping and is assumed to follow a truncated Gaussian distribution with user-specified mean and standard deviation
					  # as such, for the single-end libraries, we considered the mean fragment length to be 1.2 × the median read length and the standard deviation to be 0.1 × the mean fragment length
					  # but what is the median read length?
					  # to obtain this, we use the (v. fast) pre-processing tool fastp with default parameters (this performs some basic read cleaning). We then parse fastp's summary JSON to extract the median read length
					  # we could, if we wished, also make use of the cleaned reads in downstream processing but we choose not to - in previous work (https://www.biorxiv.org/content/10.1101/2021.11.07.467633v1), we found there was little substantive difference to expression estimates when doing so
					  
					  print SH "$fastp -i $seed_fq -o $clean_fq -j $fastp_json -h $fastp_html --thread $num_procs\n";
					  print SH "rm $fastp_html $clean_fq\n";
					  print SH "median_read_length=\$(perl $se_len $fastp_json)\n";
					  print SH "rm $fastp_json\n";
					  print SH "num1=1.2\n";
					  print SH "fragment_length=\$(echo \"\$num1 * \$median_read_length\" | bc -l | xargs printf %.0f) \n"; # bash does not support floating-point multiplication so we need to use bc for this. Discussed at https://stackoverflow.com/questions/11039876/multiplication-on-command-line-terminal
					  print SH "num2=0.1\n";
					  print SH "fragment_sd=\$(echo \"\$num2 * \$fragment_length\" | bc -l | xargs printf %.0f)\n";
					  
					  # quantify expression
					  print SH "kb count -i $index -g $t2g -x $techno -o $sample_id.$seed_num --parity $layout --fragment-l \$fragment_length --fragment-s \$fragment_sd --tcc -t $num_procs $seed_fq\n";
					
					  # delete redundant downsampled fq
					  print SH "rm $seed_fq\n";
					}
				  print SH "rm $out_dir/$species/$sample_id.$seed_num/output.bus $out_dir/$species/$sample_id.$seed_num/output.unfiltered.bus\n";
				}
			}
		  if (!(-e("$out_dir/$species/$sample_id.$seed_num/gene_tpm_matrix.txt")))
			{ print SH "cd $out_dir/$species/$sample_id.$seed_num\n";
			  print SH "$rscript $get_matrix\n";
			  print SH "cd $out_dir/$species\n";
			}
		}
	}
close(SH) or die $!;
exit 1;