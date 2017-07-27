use strict;
use warnings;
use Getopt::Std;

sub read_file
{
	my $in_name = $_[0];
	my $out_name = $_[1];

	open(INP,$in_name);	

	my $switch = 0;
	my $first = 0;
	my $out_file = "";
	while(my $line = <INP>)
	{
		if($switch > 0)
		{
			if($line =~ /^ATOM.*/)
			{
				print OUT $line;
			}
			$first = 1;
		}
		if($line =~ /^MODEL.*/)
		{
			if($first == 1)
			{ close(OUT); }
			$out_file = $out_name."_".$switch.".pdb";
			open(OUT,">$out_file");
			$switch++;
		}	
	}

	close(OUT);
	close(INP);
}

if(scalar@ARGV==0)
{
        print "\n\nScript splits pdb into separate frames but retains B-factor values.\nCalling: perl split_pdb_Bfact_CMD.pl \n\nOptions:\n-h this help message\n-i input\n-o output_name_root\n\n";
        exit;
}

my %args = ();
getopts('hi:o:',\%args);

if($args{h})
{	
        print "\n\nScript splits pdb into separate frames but retains B-factor values.\nCalling: perl split_pdb_Bfact_CMD.pl \n\nOptions:\n-h this help message\n-i input\n-o output_name_root\n\n";
        exit;
}

my $in_name = "mov.pdb";
my $out_name = "frame_.pdb";

if($args{i})
{ $in_name = $args{i}; }
if($args{o})
{ $out_name = $args{o}; }

read_file($in_name,$out_name);



