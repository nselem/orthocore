#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long 'HelpMessage';
no warnings 'experimental::smartmatch';

#####################################################################################
=head1 NAME

scan - Look for core in a reduced set!

=head1 SYNOPSIS

  --set_name    Outputfolder	
  --verbose,-v  Verbose mode
  --help,-h     This help
  --my_blast    Your allvsall blast in proper format
  --e_core      E value for consider orthologous
  --rast_ids    Required. List of genomes to be procesed 
  --mode, -m    Accepted values:g,f 
		mode -g is For gorwing it will grow the genome set until no core is found
		mode -r will remove genomes one by one and try to find Core.
 
 Please consult https://github.com/nselem/orthocore/wiki 
 
=head1 VERSION

0.01

=cut

###################################################################################
GetOptions(
        'verbose' => \my $verbose,
        'set_name=s' => \my $setname,
        'my_blast=s' => \my $myblast,
        'e_core=f'=>\(my $e_core=0.001) ,
	'initial=f'=>\(my $initial=2) ,
        'add=f'=>\(my $add=1) ,
        'list=s'=>\(my $lista=""), ##Wich genomes would you process in case you might, otherwise left empty for whole DB search
        'rast_ids=s' => \my $rast_ids,
	'mode=s'=> \(my $mode="g"),
        'help'     =>   sub { HelpMessage(0) },
        ) or HelpMessage(1);

die "$0 requires the rast_ids file (--rast_ids\nfor help type:\nscan.pl -h" unless $rast_ids;  ## A genome names list is mandatory
die "$0 requires the set_name  (--set_name\nfor help type:\nscan.pl -h" unless $setname;  ## A genome names list is mandatory


open(FILE, "$rast_ids") or die "Couldnt open $rast_ids\n $!";
my  @lines = <FILE>;
close(FILE);
my $num = @lines;

if ($mode eq "g"){
	growing($num,$rast_ids,$setname,$initial,$add);
	}
elsif ($mode eq "r"){
	removing($num,$rast_ids,$setname);
	}
else {
	print "$mode is not an accepted algorithm for mode\n";
	}
#################################### Subs ##############################################
sub removing{
	my $num=shift;
	my $rast_ids=shift;
	my $setname=shift;
	for (my $i=1;$i<=$num;$i++){
		#my $i=10;
		my $newIds="s".$i."rast_ids";
		my $newSet="s".$i."_".$setname;
		######### Print a file without some lines
		system ("perl -ne \'\$\. != $i && print\' $rast_ids > $newIds");
		print ("CoreCluster.pl -rast_ids $newIds -v -set_name $newSet -my_blast $myblast\n");
		system ("CoreCluster.pl -rast_ids $newIds -v -set_name $newSet -my_blast $myblast");
		print "##################################\n\n";
	}
}


sub growing{
	my $num=shift;
        my $rast_ids=shift;
        my $setname=shift;
        my $initial=shift;
        my $add=shift;

	for (my $i=$initial;$i<=$num;$i=$i+$add){
		 my $newIds="g".$i."rast_ids";
                 my $newSet="g".$i."_".$setname;
                 ######### Print a file without some lines
                 system (" perl -ne \'print if \$\. <= $i\'  $rast_ids > $newIds");
                 print ("CoreCluster.pl -rast_ids $newIds -v -set_name $newSet -my_blast $myblast\n");
                 system ("CoreCluster.pl -rast_ids $newIds -v -set_name $newSet -my_blast $myblast");
		 my $coresize=`grep 'Aminoacid array size' $newSet/$newSet\_Report`;
		 $coresize=~s/Aminoacid array size \= //;
		 $coresize=int ($coresize);
		 print "Core Size of this set $coresize\n";
		 if ($coresize < 1){exit;}
                 print "##################################\n\n";
		
	}
}
