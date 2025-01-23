=head

BEFORE USAGE:

# install fqtools, using instructions from https://anaconda.org/bioconda/fqtools

conda create -n fqtools
conda activate fqtools
conda install bioconda::fqtools
conda install bioconda/label/cf201901::fqtools

AFTER USAGE:

# submit to PBS like so:

#!/bin/bash
#PBS -l nodes=3:ppn=2
#PBS -l walltime=02:00:00
#PBS -N validator

source activate base
conda activate fqtools
cd /data/home/Stephen/oocyte_atlas
/data/home/Stephen/oocyte_atlas/run_fq_validator.Homo_sapiens.sh

=cut

use strict;
use warnings;

# REQUIREMENTS
my $root     = '/data/home/Stephen/oocyte_atlas';
my $species  = 'Homo_sapiens';
my $metadata = "$root/prerequisites/metadata.$species.txt"; # manually created
my $fq_dir   = "$root/fq/$species"; # from 1.download_fqs.pl
my $fatal    = 0;
if (!(-d($fq_dir))) { $fatal++; print "ERROR: cannot find $fq_dir\n";  }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_sh  = "$root/run_fq_validator.$species.sh";
open(SH,'>',$out_sh) or die $!;
print SH "#!/bin/bash\n";

# WHAT SAMPLE IDs ARE WE GOING TO VALIDATE?
my %samples = (); my %library_layouts = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $sample_id = $line[4]; my $library_layout = $line[7];
	  $library_layouts{$sample_id}{$library_layout}++;
	}
close(IN) or die $!;

# RUN FQTOOLS VALIDATE FOR EVERY SAMPLE IN $fq_dir
opendir(DIR,$fq_dir) or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my $warning = 0;
foreach my $sample_id (@sample_ids)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  
	  if (!(exists($library_layouts{$sample_id}))) { $warning++; print "WARNING: skipping $sample_id as it is not listed in $metadata\n"; }
	  next if (!(exists($library_layouts{$sample_id}))); # CHECKPOINT: skip if the sample ID is not record in $metadata
	  
	  my $num_library_layouts = scalar keys %{$library_layouts{$sample_id}};
	  if ($num_library_layouts != 1) { $warning++; print "WARNING: skipping $sample_id as $metadata states it has $num_library_layouts possible library layouts - this cannot be right\n"; }
	  next if ($num_library_layouts != 1); # CHECKPOINT: skip if the sample ID is associated with both paired-end and single-end fqs; these have already been pre-processed in some way
	  my $library_layout = '';
	  while((my $this_layout,my $irrel)=each(%{$library_layouts{$sample_id}}))
		{ $library_layout = $this_layout; }
		
	  # CHECK WHETHER WE HAVE THE RAW DATA; WARN IF WE DON'T AND ATTEMPT TO VALIDATE IF WE DO
	  my $fq1 = ''; my $fq2 = ''; my $fq = '';
	  my $files_available = 0; my %expected_files = ();
	  if ($library_layout eq 'paired')
		{ $fq1 = "$fq_dir/$sample_id/$sample_id.1.fq.gz";
		  $fq2 = "$fq_dir/$sample_id/$sample_id.2.fq.gz";
		  if ( ( (-e($fq1)) and (-e($fq2)) ) and ( (-s($fq1)) and (-s($fq2)) ) )
			{ $files_available++;
			  $expected_files{"$sample_id.1.fq.gz"}++; $expected_files{"$sample_id.2.fq.gz"}++;
			}
		}
	  elsif ($library_layout eq 'single')
		{ $fq = "$fq_dir/$sample_id/$sample_id.fq.gz";
		  if ((-e($fq)) and (-s($fq)))
			{ $files_available++;
			  $expected_files{"$sample_id.fq.gz"}++;
			}
		}
	  
	  # ALSO CHECK WHETHER THERE IS ANYTHING ELSE IN $fq_dir/$sample_id OTHER THAN THE FASTQ FILES WE'RE EXPECTING
	  # AN EXCEPTION IS MADE FOR THE FQTOOLS VALIDATE OUTPUT FILE
	  opendir(DIR,"$fq_dir/$sample_id") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  foreach my $file (@files)
		{ next if (($file eq '.') or ($file eq '..'));
		  next if ($file =~ /^$sample_id\.fqtools\_validate$/);
		  if (!(exists($expected_files{$file})))
			{ $warning++; print "WARNING: unexpected file ($file) seen in $fq_dir/$sample_id\n"; }
		}
	  
	  # IF WE HAVE THE FQ FILES, ATTEMPT TO VALIDATE
	  if ($files_available == 0)
		{ $warning++; print "WARNING: skipping $sample_id as we do not have the fqs\n";
		  if    ( ($library_layout eq 'paired') and ( (!(-e($fq1))) or (!(-e($fq2))) ) ) { print "*** at least one of the fqs for $sample_id does not exist ***\n"; }
		  elsif ( ($library_layout eq 'paired') and ( (!(-s($fq1))) or (!(-s($fq2))) ) ) { print "*** at least one of the fqs for $sample_id is incomplete ***\n";  }
		  if    ( ($library_layout eq 'single') and   (!(-e($fq)))					   ) { print "*** fq for $sample_id does not exist ***\n"; }
		  elsif ( ($library_layout eq 'single') and   (!(-s($fq)))					   ) { print "*** fq for $sample_id is incomplete ***\n";  }
		}
	  else
		{ my $out_file = "$fq_dir/$sample_id/$sample_id.fqtools_validate";
		  if ( (-e($out_file)) and (!(-s($out_file))) )
			{ print "WARNING: $out_file exists but is empty - this is a sign that validation was previously attempted, but failed\n"; }
		  my $try_again = 0;
#		  if ( (-e($out_file)) and (!(-s($out_file))) ) { print SH "rm $out_file\n"; $try_again++; }
		  next if (((-e($out_file))) and ($try_again == 0)); # CHECKPOINT: skip this sample if we have successfully validated it before
		  if    ($library_layout eq 'paired') { print SH "fqtools validate $fq1 $fq2 > $out_file\n"; }
		  elsif ($library_layout eq 'single') { print SH "fqtools validate $fq  > $out_file\n"; 	 }
		}
	}
close(SH) or die $!;

if ($warning == 0) { print "completed without warning\n"; } else { print "completed BUT WITH WARNINGS\n"; }

exit 1;