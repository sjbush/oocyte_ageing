=head

# AFTER USAGE, SUBMIT TO PBS USING THE FOLLOWING SCRIPT:

#!/bin/bash
#PBS -l nodes=3:ppn=4
#PBS -l walltime=24:00:00
#PBS -N fastp

cd /data/home/Stephen/oocyte_atlas
/data/home/Stephen/oocyte_atlas/run_fastp.Homo_sapiens.sh

=cut

use strict;
use warnings;

# REQUIREMENTS
my $root     = '/data/home/Stephen/oocyte_atlas';
my $species  = 'Homo_sapiens';
my $fq_dir   = "$root/fq/$species"; # from 1.download_fqs.pl and 2.validate_downloaded_fqs.pl
my $fastp    = '/data/home/Stephen/programs/fastp';
my $metadata = "$root/prerequisites/metadata.$species.txt"; # manually created
my $fatal    = 0;
if (!(-e($fastp)))    { $fatal++; print "ERROR: cannot find $fastp\n";    }
if (!(-d($fq_dir)))   { $fatal++; print "ERROR: cannot find $fq_dir\n";   }
if (!(-e($metadata))) { $fatal++; print "ERROR: cannot find $metadata\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 15;

# OUTPUT
my $out_dir = "$root/fastp";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_sh = "$root/run_fastp.$species.sh";
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
	  my $sample_id = $line[4]; my $library_layout = $line[7];
	  $samples{$sample_id} = $library_layout;
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
	  
	  # CHECK WHETHER WE HAVE THE RAW DATA
	  my $layout = $samples{$sample_id};
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
	  if ($files_available == 0) { print "WARNING: skipping $sample_id as we do not have the fqs\n"; }
	  next if ($files_available == 0);

	  # CHECK WHETHER THE FASTQS HAVE BEEN VALIDATED BY FQTOOLS
	  my $passes_validation = 0;
	  my $validation = "$fq_dir/$sample_id/$sample_id.fqtools_validate"; # made by 1post.validate_downloaded_fqs.pl
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
	
	  # RUN FASTP TO CLEAN THE FASTQs
	  if (!(-d("$out_dir/$species/$sample_id"))) { print SH "mkdir $out_dir/$species/$sample_id\n"; }
	  if ($layout eq 'paired')
		{ my $clean_fq1  = "$out_dir/$species/$sample_id/$sample_id.cleaned.1.fq.gz";
		  my $clean_fq2  = "$out_dir/$species/$sample_id/$sample_id.cleaned.2.fq.gz";
		  my $fastp_json = "$out_dir/$species/$sample_id/$sample_id.json";
		  my $fastp_html = "$out_dir/$species/$sample_id/$sample_id.html";
		  next if ( (-e($clean_fq1)) and (-e($clean_fq2)) and (-e($fastp_json)) ); # CHECKPOINT: skip is we have run fastp before
		  print SH "$fastp -i $fq1 -I $fq2 -o $clean_fq1 -O $clean_fq2 -j $fastp_json -h $fastp_html --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 10 --thread $num_procs\n";
		  print SH "rm $fastp_html\n";
		}
	  elsif ($layout eq 'single')
		{ my $clean_fq   = "$out_dir/$species/$sample_id/$sample_id.cleaned.fq.gz";
		  my $fastp_json = "$out_dir/$species/$sample_id/$sample_id.json";
		  my $fastp_html = "$out_dir/$species/$sample_id/$sample_id.html";
		  next if ( (-e($clean_fq)) and (-e($fastp_json)) ); # CHECKPOINT: skip is we have run fastp before
		  print SH "$fastp -i $fq -o $clean_fq -j $fastp_json -h $fastp_html --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 10 --thread $num_procs\n";
		  print SH "rm $fastp_html\n";
		}
	}
close(SH) or die $!;
exit 1;