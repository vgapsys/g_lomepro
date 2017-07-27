use strict;
use warnings;
use Getopt::Std;

sub octave_png
{
	my $dat = $_[0];
	my $png = $_[1];
	my $script = "script.m";
	my $res = $_[2];
	my $colormap = $_[3];
	my $min = $_[4];
	my $max = $_[5];
	my $lcmap = $_[6];

	#create a script
        unless(open(SCRIPT,">$script"))
        {       print "Cannot open $script\n" and exit; }
	
		print SCRIPT "a = load(\"$dat\");";
		print SCRIPT "figure(\"visible\",\"off\");";
		if($min!=-42 && $max!=-42)
		{
			print SCRIPT "imagesc(a,limits=[$min,$max]);";
		}
		elsif($min!=-42 && $max==-42)
		{
                        print SCRIPT "imagesc(a,limits=[$min,max(max(a))]);";
		}
		elsif($min==-42 && $max!=-42)
		{
                        print SCRIPT "imagesc(a,limits=[min(min(a)),$max]);";
		}
		else
		{
                        print SCRIPT "imagesc(a,limits=[min(min(a)),max(max(a))]);";
		}
                print SCRIPT "set(gca,'ydir','normal');";
		if($lcmap)
		{
			print SCRIPT "colormap(load('$lcmap'));";
		}
		else
		{
               		print SCRIPT "colormap($colormap);";
		}

                print SCRIPT "colorbar;";
		print SCRIPT "print('$png','-dpng','-r$res');";
	close(SCRIPT);

	#execute script
	system("octave --silent $script");

	#remove script
	system("rm $script");
}


sub matlab_png
{
        my $dat = $_[0];
        my $png = $_[1];
        my $script = "script.m";
        my $res = $_[2];
	my $colormap = $_[3];
	my $min = $_[4];
        my $max = $_[5];
	my $lcmap = $_[6];

        #create a script
#	system("unset DISPLAY");
        unless(open(SCRIPT,">$script"))
        {       print "Cannot open $script\n" and exit; }

                print SCRIPT "a = load(\'$dat\');";
                print SCRIPT "figure(\'visible\',\'off\');";
                print SCRIPT "imagesc(a);";
                if($min!=-42 && $max!=-42)
                {
			print SCRIPT "set(gca,'CLim',[$min,$max]);";
                }
                elsif($min!=-42 && $max==-42)
                {
			print SCRIPT "set(gca,'CLim',[$min,max(max(a))]);";
                }
                elsif($min==-42 && $max!=-42)
                {
			print SCRIPT "set(gca,'CLim',[min(min(a)),$max]);";
                }
                else
                {
			print SCRIPT "set(gca,'CLim',[min(min(a)),max(max(a))]);";
                }
		print SCRIPT "set(gca,'ydir','normal');";
                if($lcmap)
                {
                        print SCRIPT "colormap(load('$lcmap'));";
                }
                else
                {
                        print SCRIPT "colormap($colormap);";
                }
		print SCRIPT "colorbar;";
                print SCRIPT "print('$png','-dpng','-r$res');";
        close(SCRIPT);
	
        #execute script
        system("matlab -nodisplay -nosplash -nodesktop -r \"script;quit;\" > /dev/null 2>&1");
	system("stty echo");

        #remove script
        system("rm $script");
}


sub read_file
{
	my $input = $_[0];
	my $png = $_[1];
	my $flag = lc($_[2]);
	my $res = $_[3];
	my $colormap = $_[4];
	my $min = $_[5];
	my $max = $_[6];
	my $lcmap = $_[7];

	if($flag =~ /octave/)
	{
		octave_png($input,$png,$res,$colormap,$min,$max,$lcmap);
	}
	else #matlab
	{
		matlab_png($input,$png,$res,$colormap,$min,$max,$lcmap);
	}
}


sub do_well
{
#do_well($inp,$path,$xsoft,$res,$cmap,$min,$max,$lcmap);

	my $input = $_[0];
	my $path = $_[1];
	my $flag = lc($_[2]);
	my $res = $_[3];
	my $colormap = $_[4];
	my $min = $_[5];
	my $max = $_[6];
	my $lcmap = $_[7];

	unless(open(FILE,$input))
	{	print "Cannot open $input\n" and exit; }


	my $line = "";
	my $file = "";
	my $out_name = "";
	my $png_name = ""; #may not be used, depends on the flag
	my $frame_count = 0;
	while($line = <FILE>)
	{
		if($line =~ /^FRAME.*/)
		{
			#close previous file, doesn't hold for the first file
			if($frame_count > 0)
			{
				close(OUT);
				if($flag =~ /octave/) #octave
				{
					octave_png($out_name,$png_name,$res,$colormap,$min,$max,$lcmap);
				}
				elsif($flag =~ /matlab/) #matlab
				{
					matlab_png($out_name,$png_name,$res,$colormap,$min,$max,$lcmap);
				}
			}

			print "PROCESSING: $line";

                        chomp($line);
                        my @arr = split(/[\s\t]{1,}/,$line);
                        $file = "frame_$frame_count.dat";
                        $out_name = "$path/$file";
                        $png_name = "$path/frame_$frame_count.png";

		        unless(open(OUT,">$out_name"))
		        {       print "Cannot open $file\n" and exit; }

			$frame_count++;
		}
		else
		{
			print OUT $line;
		}
	}

	close(OUT);
	# png for the last frame
        if($flag =~ /octave/) #octave
        {
        	octave_png($out_name,$png_name,$res,$colormap,$min,$max,$lcmap);
        }
        elsif($flag =~ /matlab/) #matlab
        {
                matlab_png($out_name,$png_name,$res,$colormap,$min,$max,$lcmap);
        }

	close(FILE);
}


if(scalar@ARGV==0)
{	
	print "\n\nScript extracts and splits frames from a matrix file with many frames. Can convert dat files to png images: for that add a flag matlab or octave.\nCalling: perl split_mov_mat_CMD.pl\n\nOptions:\n-h this help message\n-i input_matrix \n-p /path/to/output/ \n-x choose \"matlab\" or \"octave\" (default matlab)\n-r resolution_dpi (default 100)\n-c colormap, e.g. jet, gray... [optional]\n-l load_external_colormap (provide filename)\n-a min_value [optional]\n-z max_value [optional]\n\nREMARK: every frame is printed into a new file, therefore many files may be produced\n\n";
	exit;
}

my %args = ();
getopts('hi:p:x:r:c:l:a:z:',\%args);

if($args{h})
{
        print "\n\nScript extracts and splits frames from a matrix file with many frames. Can convert dat files to png images: for that add a flag matlab or octave.\nCalling: perl split_mov_mat_CMD.pl\n\nOptions:\n-h this help message\n-i input_matrix \n-p /path/to/output/ \n-x choose \"matlab\" or \"octave\" (if none is chosen, image is not produced)\n-r resolution_dpi (default 100)\n-c colormap, e.g. jet, gray... [optional]\n-l load_external_colormap (provide filename)\n-a min_value [optional]\n-z max_value [optional]\n\nREMARK: every frame is printed into a new file, therefore many files may be produced\n\n";
        exit;
}

my $inp = "mat.dat";
my $path = "./";
my $xsoft = "none";
my $res = "100";
my $cmap = "jet";
my $lcmap = "";
my $min = -42;
my $max = -42;

if($args{i})
{ $inp = $args{i}; }
if($args{p})
{ $path = $args{p}; }
if($args{x})
{ $xsoft = $args{x}; }
if($args{r})
{ $res = $args{r}; }
if($args{c})
{ $cmap = $args{c}; }
if($args{l})
{ $lcmap = $args{l}; }
if(exists $args{a})
{ $min = $args{a}; }
if(exists $args{z})
{ $max = $args{z}; }

do_well($inp,$path,$xsoft,$res,$cmap,$min,$max,$lcmap);







