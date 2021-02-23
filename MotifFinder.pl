#!/usr/bin/perl
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Tree::Suffix;
use Enumeration qw(seed);

my ($gap,$dir,$output,$threads,$window,$ref,$extension,$mcs,$background);

GetOptions(
	'd=s'  => \$dir,
	'm=i'  => \$gap,
	't=i'  => \$threads,
	'o=s'  => \$output,
	'w=i'  => \$window,
	'r=s'  => \$ref,
	'e=i'  => \$extension,
	'c=i'  => \$mcs,
	'b=s'  => \$background
);

unless(defined $dir){print <DATA>;exit}

$gap //= 11;
$threads //= 11;
$output //= 'motif.out2';
$window //= 10;
$ref //= 'Scer';
$extension //= 4;
$mcs //= 5;
$background //= 'promoter';


$dir =~s/\/$//;

print STDERR "\n";
seed($dir,$gap,$window,$threads,$ref,$output,$extension,$mcs,$background);





__DATA__

Useage:
	
	perl MotifFinder.pl [options]

[Options]

	-d	dir for .aln files [Required]
	-m	Gaps for seed search (XYZ.UVW) [11]
	-t 	threads	[12]
	-o	output [motif.out]
	-w	distant of motif movement [10], 10bp to the left and to the right
	-r	reference speceis [Scer]
	-e	extension of both ends : 4bp
	-c	cutoff for MCS: 5
	-b  	backgroud: ('promoter','coding','terminator') ['promoter']
	
	



