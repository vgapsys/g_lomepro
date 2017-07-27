use strict;
use warnings;
use Getopt::Std;

sub do_well
{
	my $inp = $_[0];
	my $out = $_[1];
	my $first = $_[2];
	my $last = $_[3];

	open(INP,$inp);
	open(OUT,">$out");

	my $num = 0;
	my $mid = 0;
	
	while(my $line = <INP>)
	{
		if($line =~ /FRAME (\d{1,})/g)
		{
			my $num = $1;
			if($num>=$first && $num<=$last)
			{
				$mid = 1;
				print OUT $line;
			}
			elsif($num>$last)
			{
				last;
			}
		}
		elsif($mid==1)
		{
			print OUT $line;
		}

	}
	
	close(INP);
	close(OUT);
}

if(scalar@ARGV==0)
{
        print "\n\nScript cuts out a part of mov_mat.dat file.\nCalling: perl trjconv_mat.pl \n\nOptions:\n-h this help message\n-i input.mat\n-o output.mat\n-b first_frame \n-e last_frame\n\n";
        exit;
}

my %args = ();
getopts('hi:o:b:e:',\%args);

if($args{h})
{	
        print "\n\nScript cuts out a part of mov_mat.dat file.\nCalling: perl trjconv_mat.pl \n\nOptions:\n-h this help message\n-i input.mat\n-o output.mat\n-b first_frame \n-e last_frame\n\n";
	exit;
}

my $inp = "";
my $out = "out.dat";
my $first = -42;
my $last = 9999999;

if($args{i})
{ $inp = $args{i}; }
if($args{o})
{ $out = $args{o}; }
if(exists $args{b})
{ $first = $args{b}; }
if(exists $args{e})
{ $last = $args{e}; }

do_well($inp,$out,$first,$last);


