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

if(scalar@ARGV==0)
{	
	print "\n\nScript converts .dat to .png.\nCalling: perl dat2png.pl \n\nOptions:\n-h this help message\n-i input_matrix \n-o out.png \n-x choose \"matlab\" or \"octave\" (default matlab)\n-r resolution_dpi (default 100)\n-c colormap, e.g. jet, gray... [optional]\n-l load_external_colormap (provide filename)\n-a min_value [optional]\n-z max_value [optional]\n\n";
	exit;
}

my %args = ();
getopts('hi:o:x:r:c:l:a:z:',\%args);

if($args{h})
{
        print "\n\nScript converts .dat to .png.\nCalling: perl dat2png.pl \n\nOptions:\n-h this help message\n-i input_matrix \n-o out.png \n-x choose \"matlab\" or \"octave\" (default matlab)\n-r resolution_dpi (default 100)\n-c colormap, e.g. jet, gray... [optional]\n-l load_external_colormap (provide filename)\n-a min_value [optional]\n-z max_value [optional]\n\n";
        exit;
}

my $inp = "mat.dat";
my $out = "out.png";
my $xsoft = "matlab";
my $res = "100";
my $cmap = "jet";
my $lcmap = "";
my $min = -42;
my $max = -42;

if($args{i})
{ $inp = $args{i}; }
if($args{o})
{ $out = $args{o}; }
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

read_file($inp,$out,$xsoft,$res,$cmap,$min,$max,$lcmap);



