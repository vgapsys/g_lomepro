use strict;
use warnings;
use POSIX;
use Getopt::Std;

sub create
{
	my $out = $_[0];
	my $num = $_[1];

	my $denom = floor($num);
	my $b = 1.0;
	my $g = 1.0;
	my $r = 1.0;
	my $delta = 1/$denom;

	open(OUT,">$out");
	for(my $i=0; $i<$num; $i++)
	{
			$r-=$delta;
			$g-=$delta;

		print OUT "$r	$g	$b\n";
	}
	close(OUT);
}

if(scalar@ARGV==0)
{
        print "\n\nScript creates 'wb' colormap.\nCalling: perl create_wb_colormap.pl \n\nOptions:\n-h this help message\n-o cmap.dat\n-n number_of_rows\n\n";
        exit;
}

my %args = ();
getopts('ho:n:',\%args);

if($args{h})
{
        print "\n\nScript creates 'wb' colormap.\nCalling: perl create_wb_colormap.pl \n\nOptions:\n-h this help message\n-o cmap.dat\n-n number_of_rows\n\n";
        exit;
}

my $out = "cmap.dat";
my $num = 64;

if($args{o})
{ $out = $args{o}; }
if($args{n})
{ $num = $args{n}; }

create($out,$num);


