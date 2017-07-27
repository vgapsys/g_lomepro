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
	my $bot = $_[2];
	my $out = $_[3];

	open(MAT1,$top);
	open(MAT2,$bot);
	open(PDB,$inp);
	open(OUT,">$out");

	my $mat1_fr;
	my $mat2_fr;
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
                                $mat1_fr = read_frame(\*MAT1,$first,0);
                                $mat2_fr = read_frame(\*MAT2,$first,1);
			}

			#top leaflet
			$line =~ s/(.*)\n/$1/;
			my @arr = split(/[\s\t]{1,}/,$line);
			$foo = shift(@$mat1_fr);
	                $pdbform = sprintf "%-6s%5u  %-4.4s%3.3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
          	                   $arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8],$arr[9],$foo;
	                print OUT $pdbform;

			#bottom leaflet
			$line = <PDB>;
                        $line =~ s/(.*)\n/$1/;
                        @arr = split(/[\s\t]{1,}/,$line);
                        $foo = shift(@$mat2_fr);
                        $pdbform = sprintf "%-6s%5u  %-4.4s%3.3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                                   $arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5],$arr[6],$arr[7],$arr[8],$arr[9],$foo;
                        print OUT $pdbform;
		}
		else
		{
			if($line =~ /MODEL/g)
			{
				$mat1_fr = read_frame(\*MAT1,$first,0);
				$mat2_fr = read_frame(\*MAT2,$first,1);
				$first = 0;
			}
			print OUT $line;
		}
	}

	close(OUT);
	close(PDB);
	close(MAT1);
	close(MAT2);
}

if(scalar@ARGV==0)
{
        print "\n\nScript inserts values in the b-factor fields of both bilayers from the top and bottom leaflet data matrices.\nCalling: perl mat2pdb.pl \n\nOptions:\n-h this help message\n-i input.pdb\n-u upper_leaflet.mat\n-d down_leaflet.mat \n-o output.pdb\n\n";
        exit;
}

my %args = ();
getopts('hi:u:d:o:',\%args);

if($args{h})
{	
        print "\n\nScript inserts values in the b-factor fields of both bilayers from the top and bottom leaflet data matrices.\nCalling: perl mat2pdb.pl \n\nOptions:\n-h this help message\n-i input.pdb\n-u upper_leaflet.mat\n-d down_leaflet.mat \n-o output.pdb\n\n";
	exit;
}

my $inp = "";
my $top = "";
my $bot = "";
my $out = "out.pdb";

if($args{i})
{ $inp = $args{i}; }
if($args{u})
{ $top = $args{u}; }
if($args{d})
{ $bot = $args{d}; }
if($args{o})
{ $out = $args{o}; }

do_well($inp,$top,$bot,$out);


