use strict;
use warnings;
use Getopt::Std;

sub read_frame
{
	my $file = $_[0];
	my $first = $_[1];
	my $reverse = $_[2];

	my $line = "";
	if($first) #read line with FRAME
	{
		$line = <$file>;	
	}
	
	$line = "";
	# read further
	my @res = ();
	while($line = <$file>)
	{
		if($line =~ /FRAME/g)
		{
			last;
		}
		$line =~ s/(.*)\n/$1/;
		my @arr = split(/[\s\t]{1,}/,$line);
		foreach my $el(@arr)
		{
			if($reverse)
			{
				unshift(@res,$el);
			}
			else
			{
				push(@res,$el);
			}
		}
	}
	return(\@res);
}

sub do_well
{
	my $inp = $_[0];
	my $top = $_[1];
	my $out = $_[2];
	my $down = $_[3];

	open(MAT1,$top);
	open(PDB,$inp);
	open(OUT,">$out");

	my $mat1_fr;
	my $first = 1;
	my $foo = 0.0;
	my $pdbform = "";

	while(my $line = <PDB>)
	{
		if($line =~ /ATOM/g)
		{
			#check if there is only one frame in the pdb, i.e. pdb of an averaged property
			if($first==1)
			{
				$first=0;
                                $mat1_fr = read_frame(\*MAT1,$first,$down);
			}

			$line =~ s/(.*)\n/$1/;
			my @arr = split(/[\s\t]{1,}/,$line);
			$foo = shift(@$mat1_fr);
	                $pdbform = sprintf "%-6s%5u  %-4.4s%3.3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
          	                   $arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8],$arr[9],$foo;
	                print OUT $pdbform;

		}
		else
		{
			if($line =~ /MODEL/g)
			{
				$mat1_fr = read_frame(\*MAT1,$first,$down);
				$first = 0;
			}
			print OUT $line;
		}
	}

	close(OUT);
	close(PDB);
	close(MAT1);
}

if(scalar@ARGV==0)
{
        print "\n\nScript inserts values in the b-factor fields of one leaflet data matrix.\nCalling: perl mat2pdb_single_leaflet.pl \n\nOptions:\n-h this help message\n-i input.pdb\n-l leaflet.mat\n-d put this flag if using lower leaflet\n-o output.pdb\n\n";
        exit;
}

my %args = ();
getopts('hi:l:o:d',\%args);

if($args{h})
{	
        print "\n\nScript inserts values in the b-factor fields of one leaflet data matrix.\nCalling: perl mat2pdb_single_leaflet.pl \n\nOptions:\n-h this help message\n-i input.pdb\n-l leaflet.mat\n-d put this flag if using lower leaflet\n-o output.pdb\n\n";
	exit;
}

my $inp = "";
my $top = "";
my $down = 0;
my $out = "out.pdb";

if($args{i})
{ $inp = $args{i}; }
if($args{l})
{ $top = $args{l}; }
if($args{o})
{ $out = $args{o}; }
if($args{d})
{ $down = 1; }

do_well($inp,$top,$out,$down);

