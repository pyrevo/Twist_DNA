#!/usr/bin/perl -w
# Extracts information från Picard HsMetrics and samtools statistics and combines to a qc report.
use strict;
use warnings;
use Time::Piece;

#my $HsMinfile1=shift @ARGV;

#my $name=$HsMinfile1;
#my $iterationnumber=$HsMinfile1;
#$iterationnumber =~ s/^[^_]*_//;
#$iterationnumber =~ s/.bam_HsM.txt//;
#$name =~ s/-ready_$iterationnumber.bam_HsM.txt//;
#my $Bamtoolsinfile1=$name."_samtools_stats_".$iterationnumber.".txt";
#my $outfile=$name."_qc_rapport.txt";
my $statfile=shift @ARGV;
my $outfile=shift @ARGV;


open OUT,">$outfile";
# print OUT "*********************************************\n";
# print OUT "QC rapport NGS\n";
# print OUT "DNAnr: $name\n";
# print OUT "Metod: TWIST Exome\n";
#
# my $date = localtime-> ymd;
# print OUT "Datum: $date\n";
# print OUT "vcf: $name-gatk-haplotype_cg_$iterationnumber.vcf\n";
# print OUT "*********************************************\n";
# print OUT "Bioinformatisk analys\n\n";

my @data1;
my $AA_count=0;
my $RA_count;
my $hetero;
my $ratio;
my $total_count;
open(IN1, '<', $statfile)
#open IN1,$statfile;
while (<IN1>) {
    chomp;
   	s/\,//g;
    @data1=split(/\s+/, $_);

    if ($_ =~ m/undef/) {
    	$RA_count = 1;
    	$AA_count = 0;
    	$total_count = 3;
    	next;
    }
    if ($_ =~ m/'count'/) {
    	$total_count=$data1[3];
    }
    if ($_ =~ m/het_RA_count/) {
        $RA_count=$data1[3];
    }
    if ($_ =~ m/het_AA_count/) {
        $AA_count=$data1[3];
   }
}
$hetero=$RA_count + $AA_count;
$ratio=$hetero/$total_count;
my $ratio_rounded = sprintf("%.2f", $ratio);
if ($ratio < 0.22) {
   	print OUT "Kön baserat på analys: M ($ratio_rounded)\n";
}
elsif ($ratio > 0.45 ) {
   	print OUT "Kön baserat på analys: K ($ratio_rounded)\n";
}
else {
	print OUT "Kön baserat på analys: Kan ej avgöras ($ratio_rounded)\n";
}


# my @data2;
# my @data2_percentage;
#
# open IN2,$Bamtoolsinfile1;
# while (<IN2>) {
#     chomp;
#
#     if ($_ =~ m/SN\tsequences:/) {
#         @data2=split(/\s+/, $_);
# 	@data2_percentage=split(/\s+/, $_);
#         my $total_reads=$data2[2];
#         print OUT "Antal reads: $total_reads\n";
#     }
#     if ($_ =~ m/SN\treads\ mapped:/) {
#     	@data2=split(/\s+/, $_);
#    	my $mapped_reads=$data2[3];
#     	my $pct_mapped_reads=eval sprintf('%.2f',$data2[3]/$data2_percentage[2]*100);
# 	$pct_mapped_reads="(".$pct_mapped_reads."%)";
#     	print OUT "Antal mappade reads:\t$mapped_reads\t$pct_mapped_reads\n";
#     }
# }
#
# my @data3;
# open IN3,$HsMinfile1;
# while (<IN3>) {
#     chomp;
#     if ($_ =~ m/^s/) { next; }
#     if ($_ =~ m/^#/) { next; }
#     if ($_ =~ m/^B/) { next; }
#     if ($_ =~ m/^$/) { next; }
#     @data3=split(/\s+/, $_);
#     my $andel_reads=$data3[17];
#     my $total_mean=$data3[21];
#     my $total_20x=eval sprintf('%.6f',$data3[29]);
#     my $total_10x=eval sprintf('%.6f',$data3[28]);
#     print OUT "Andel reads på target: $andel_reads\n";
#     print OUT "Täckningsdjup medel (x): $total_mean\n";
#     print OUT "Täckning (% ≥10x): $total_10x\n";
#     print OUT "Täckning (% ≥20x): $total_20x\n";
#     last;
# }
# print OUT "\n";
# print OUT "*********************************************\n";
# my $whoami = `bash -i -c whoami`;
# chomp $whoami;
# print OUT "Analys av $whoami";
# close IN2;
# close IN3;
close IN1;
close OUT;
