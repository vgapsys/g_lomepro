use strict;
use warnings;
use Getopt::Std;

sub read_file
{
	my $path = $_[0];
	my $num_init = $_[1];
	my $num_last = $_[2];
	my $out_name = $_[3];
	my $name = $_[4];

	my $png_string = "";
	for(my $i=$num_init; $i<=$num_last; $i++)
	{
		my $frame = "$path/$name$i.png";
                if($i==$num_init)
                {
                        $png_string = "$frame";
                }
                else
                {
                        $png_string .= ",$frame";
                }
	}

	my $run_string = "mencoder 'mf://$png_string' -mf type=png:fps=10 -ovc copy -oac copy -o $out_name";

	system("$run_string");
}

if(scalar@ARGV==0)
{
        print "\n\nScript runs mencoder to merge PNGs into an avi movie in a correct order of frames.\nCalling: perl make_a_movie.pl \n\nOptions:\n-h this help message\n-p /path/to/the/png/folder/\n-b first_frame_number\n-e last_frame_number\n-n root_of_the_png_file_name\n-o output_file\n\n";
        exit;
}

my %args = ();
getopts('hp:b:e:o:n:',\%args);

if($args{h})
{	
        print "\n\nScript runs mencoder to merge PNGs into an avi movie in a correct order of frames.\nCalling: perl make_a_movie.pl \n\nOptions:\n-h this help message\n-p /path/to/the/png/folder/\n-b first_frame_number\n-e last_frame_number\n-n root_of_the_png_file_name\n-o output_file\n\n";
	exit;
}

my $path = ".";
my $num_init = 0;
my $num_last = 1;
my $out_name = "out.avi";
my $name = "frame_";

if($args{p})
{ $path = $args{p}; }
if(exists $args{b})
{ $num_init = $args{b}; }
if(exists $args{e})
{ $num_last = $args{e}; }
if($args{o})
{ $out_name = $args{o}; }
if($args{n})
{ $name = $args{n}; }

read_file($path,$num_init,$num_last,$out_name,$name);



