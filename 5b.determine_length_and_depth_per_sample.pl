=head

# SUPPLEMENTARY FIGURE 1 of the manuscript: the characteristics of each sample (read depth, read length, Q20 rate, and mitochondrial content)

library(grid)
library(ggplot2)
theme_set(theme_bw())

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/read_length_and_depth_per_sample.Homo_sapiens.txt',header=T,sep='\t') # from 5b.determine_length_and_depth_per_sample.pl
df.sub<-subset(df,(df$Maturation.stage == 'MII' & df$Publication != 'Li 2019' & df$Publication != 'Qi 2020'))
df<-df.sub
fig1 <- ggplot(df, aes(x = Sample.ID, y = Total.reads.before.filtering, colour = factor(Publication))) + geom_point() + xlab('Sample') + ylab('Number of reads') + labs(colour='Study') + theme(axis.text.x=element_blank()) + scale_y_log10(limits=c(100000,160000000),breaks=c(100000,1000000,5000000,10000000,50000000,100000000),labels=c("100,000","1 million","5 million","10 million","50 million","100 million")) + ggtitle('A. Read depth per sample')

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/read_length_and_depth_per_sample.Homo_sapiens.txt',header=T,sep='\t') # from 5b.determine_length_and_depth_per_sample.pl
df.sub<-subset(df,(df$Maturation.stage == 'MII' & df$Publication != 'Li 2019' & df$Publication != 'Qi 2020'))
df<-df.sub
fig2 <- ggplot(df, aes(x = Sample.ID, y = Read.length.before.filtering, colour = factor(Publication))) + geom_point() + xlab('Sample') + ylab('Read length (bp)') + labs(colour='Study') + theme(axis.text.x=element_blank()) + ggtitle('B. Read length per sample')

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/read_length_and_depth_per_sample.Homo_sapiens.txt',header=T,sep='\t') # from 5b.determine_length_and_depth_per_sample.pl
df.sub<-subset(df,(df$Maturation.stage == 'MII' & df$Publication != 'Li 2019' & df$Publication != 'Qi 2020'))
df<-df.sub
fig3 <- ggplot(df, aes(x = Publication, y = Q20.rate.before.filtering)) + geom_boxplot() + xlab('Study') + ylab('Q20 rate') + ggtitle('C. Q20 rate per sample')

df<-read.table('C:/Users/User/Desktop/oocyte_atlas/results/pct_mitochondrial_expression_per_sample.txt',header=T,sep='\t') # from 7a.make_tpm_matrix.pl
df.sub<-subset(df,(df$Maturation.stage == 'MII' & df$Study != 'Li 2019' & df$Study != 'Qi 2020'))
fig4 <- ggplot(df, aes(x = Sample.ID, y = X..mitochondrial.expression.per.million.transcripts, colour = factor(Study))) + geom_point() + xlab('Sample') + ylab('% mitochondrial reads') + labs(colour='Study') + theme(axis.text.x=element_blank()) + ggtitle('D. Mitochondrial content per sample')

png(file="C:/Users/User/Desktop/oocyte_atlas/Documents/figure_S1.png",width=10,height=16,units='in',res=300) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow=4,ncol=2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(fig1, vp = define_region(row=1,col=1:2))
print(fig2, vp = define_region(row=2,col=1:2))
print(fig3, vp = define_region(row=3,col=1:2))
print(fig4, vp = define_region(row=4,col=1:2))
dev.off()

=cut

use strict;
use warnings;
use JSON qw(decode_json);

# REQUIREMENTS
my $species  = 'Homo_sapiens';
my $root     = '/data/home/Stephen/oocyte_atlas'; # 'C:/Users/Administrator/Desktop/oocyte_atlas';
my $in_dir   = "$root/fastp/$species"; # from 5a.run_fastp_to_determine_length_and_depth.pl
my $metadata = "$root/prerequisites/metadata.$species.txt"; # manually created
my $fatal    = 0;
if (!(-d($in_dir)))   { $fatal++; print "ERROR: cannot find $in_dir\n";   }
if (!(-e($metadata))) { $fatal++; print "ERROR: cannot find $metadata\n"; }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file = "$root/read_length_and_depth_per_sample.$species.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Sample ID\tPublication\tMaturation stage\tRead length before filtering\tTotal reads before filtering\tQ20 rate before filtering\tRead length after filtering\tTotal reads after filtering\tQ20 rate after filtering\n";

# WHAT SAMPLE IDs ARE WE GOING TO RUN?
my %samples = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $stage = $line[1]; my $sample_id = $line[4];
	  $samples{$sample_id}{study} = $study;
	  $samples{$sample_id}{stage} = $stage;
	}
close(IN) or die $!;

opendir(DIR,$in_dir) or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
my $samples_seen = 0; my $samples_total = @sorted_sample_ids; $samples_total = $samples_total-2;
foreach my $sample_id (@sorted_sample_ids)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  my $file = "$in_dir/$sample_id/$sample_id.json";
	  if (!(-e($file))) { print "ERROR: unable to find $file\n"; }
	  next if (!(-e($file)));
	  $samples_seen++; print "$samples_seen of $samples_total...\n";
	  open(JSON,'<',$file) or die $!;
	  local $/;
	  my $fastp = decode_json(<JSON>);
	  close(JSON);
	  my $raw_avg_read_length   = $fastp->{'summary'}{'before_filtering'}{'read1_mean_length'};
	  my $raw_depth 		    = $fastp->{'summary'}{'before_filtering'}{'total_reads'};
	  my $raw_q20 			    = $fastp->{'summary'}{'before_filtering'}{'q20_rate'};
	  my $clean_avg_read_length = $fastp->{'summary'}{'after_filtering'}{'read1_mean_length'};
	  my $clean_depth 		    = $fastp->{'summary'}{'after_filtering'}{'total_reads'};
	  my $clean_q20 			= $fastp->{'summary'}{'after_filtering'}{'q20_rate'};
	  my $study			  	    = $samples{$sample_id}{study};
	  my $stage			  	    = $samples{$sample_id}{stage};
	  print OUT "$sample_id\t$study\t$stage\t$raw_avg_read_length\t$raw_depth\t$raw_q20\t$clean_avg_read_length\t$clean_depth\t$clean_q20\n";
	}
close(OUT) or die $!;
exit 1;