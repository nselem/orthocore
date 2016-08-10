#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long 'HelpMessage';
no warnings 'experimental::smartmatch';
######################################################################
###This is the main script to run the ortho-group tools
######################################################################

=head1 NAME

Corecluster - get license texts at the command line!

=head1 SYNOPSIS

  --list,-l     Holder name (required)
  --num,-n       License year (defaults to current year)
  --help,-h       Print this help

=head1 VERSION

0.01

=cut

#print "This script will help you to find a core between genomes\n";
#print "How many genomes would you like to process?\n";

##################################################################
######## Subs on this file #######################################
sub cleanFiles;
sub specialCluster;
sub printVariables;
sub getDrawInputs;

##################################################################
############# Input variables

GetOptions(
        'verbose' => \my $verbose,
        'set_name=s' => \my $setname,
        'my_blast=s' => \my $myblast,
        'e_value=f'=> \(my $e_value=0.000001), 		# E value. Minimal for a gene to be considered a hit.
        'bitscore=i'=>\(my $bitscore=0),  		## Revisar el archivo .BLAST.pre para tener idea de este parámetro.
        'cluster_radio=i'=>\(my $cluster_radio=10), 	#number of genes in the neighborhood to be analized
        'e_cluster=f'=>\(my $e_cluster=0.001), #Query search e-value for homologies from reference cluster, values above this will be colored
        'e_core=f'=>\(my $e_core=0.001) ,  
        'list=s'=>\(my $lista=""), ##Wich genomes would you process in case you might, otherwise left empty for whole DB search
        'rescale=i'=>\(my $rescale = 85000) ,
        'num=i'=>\my $num ,  #the number of genomes to be analized in case you used the option $LIST, comment if $LIST is empty
        'rast_ids=s' => \my $rast_ids,
	'help'     =>   sub { HelpMessage(0) },
        ) or HelpMessage(1);


die "$0 requires the rast_ids file (--rast_ids\nfor help type:\ncore.pl -h" unless $rast_ids;  ## A genome names list is mandatory
die "$0 requires the set_name  (--set_name\nfor help type:\ncore.pl -h" unless $setname;  ## A genome names list is mandatory

my $list_all=get_lista($lista,$verbose,$rast_ids);
$num=get_number($lista,$list_all,$rast_ids);
$lista=$list_all;

my $dir=&Cwd::cwd();            ##The path of your directory
#my $name=pop @{[split m|/|, $dir]};             ##The path of your directory
my $name=$setname;             ##The path of your directory

my $blast;
if ($myblast){
	$blast=$myblast;
}
else{
	## Here run header and blast Pending
	my $blast="$name.blast";
}


printVariables($verbose);

#####################################################################
########## Main ######################################################

my @LISTA=split(",",$list_all);
my $outname=$setname;
#$outname=~s/\.query//;
if(!-e $outname) {system("mkdir $outname");}
if ($verbose ){print "Your courrent directory: $name\n";}

my $report="";
if (-e "$outname/$outname\_Report"){`rm $outname/$outname\_Report`;}
$report=$report."Set name $setname\t
                My blast $myblast\t
		e_core $e_core\t
		list $lista\t
		number $num\t
		name folder $name\t
		dir $dir\t
		blast $blast\t";

		#print "
                #e_value $e_value\t
                #bitscore $bitscore\t
		#cluster radio $cluster_radio\t
		#rescale $rescale\t
		#";

	print "Searching genetic core on selected clusters\n";
	print"2_OrthoGroups.pl -e_core $e_core -list $lista -num $num -rast_ids $rast_ids -outname $outname -name $name -blast $myblast\n ";
	#print "Enter to continue\n";
	#my $pause=<STDIN>;
	if($blast){
		system("2_OrthoGroups.pl -e_core $e_core -list $lista -num $num -rast_ids $rast_ids -outname $outname  -name $name -blast $myblast");
	}
	else{
		system("2_OrthoGroups.pl -e_core $e_core -list $lista -num $num -rast_ids $rast_ids -outname $outname  -name $name ");
	}
	print "Core finished!\n\n";
	my $boolCore= `wc -l <$outname/Core`;
	chomp $boolCore;
	#$boolCore=~s/[^0-9]//g;
	$boolCore=int($boolCore);
	print "Elements on core: $boolCore!\n";
#____________________________________________________________________________________________________________
if ($boolCore>=1){
	print "There is a core with at least two genes on this cluster\n";
	$report=$report."\nThere is a core composed by $boolCore orhtolog on this cluster\n";
	$report=$report. "Enzyme functions on reference organisms are given by:\n";
	## Obteniendo el cluster del organismo de referenecia mas parecido al query
	# Abrimos los input files de ese organismo y tomamos el de mejor score	
	#my $specialCluster=specialCluster($special_org);
	#print "Best cluster $specialCluster\n";
       	#my $functions=`cut -f1,2 $outname/FUNCTION/$specialCluster.core.function `;
#       	print "cut -f1,2 $name/FUNCTION/$specialCluster.core.function ";
#	print "Function $functions#\n";
	#$report=$report."\n".$functions;
	print "Aligning...\n";
	print "system (perl multiAlign_gb.pl $num $lista $outname)\n";
	#print "Enter to continue\n";
	#my $paus=<STDIN>;
	system ("multiAlign_gb.pl $num $lista $outname");
	print "Sequences were aligned\n\n";

	print "Creating aminoacid core cluster matrix..\n";
	system("ChangeName.pl $outname");
	system("EliminadorLineas.pl $outname");

	system("Concatenador.pl $outname");
	system("Rename_Ids_Star_Tree.pl $rast_ids $outname");
	my $line =`perl -ne \'print if \$\. == 2\' $outname/RightNames.txt `;
	#print "`perl -ne 'print if \$\. == 2' RightNames.txt `";
	#print "Line $line\n";
 	my $len = map $_, $line =~ /(.)/gs;
	$len--;
	$report=$report."\nAminoacid array size = $len \n\n";
	print "Formating matrix..\n";
	system ("converter.pl $outname/RightNames.txt");

	print "Constructing a tree with quicktree with a 100 times bootstrap\n";
	system "quicktree -i a -o t -b 100 $outname/RightNames.stockholm > $outname/BGC_TREE.tre";
# 	$orderFile="$outname/$outname\_BGC_TREE.order";
#	print "I will draw SVG clusters with concatenated tree order\n";
#	$INPUTS=getDrawInputs($orderFile);
	}
	else{  ### If there is no core, then sort according to principal hits
		$report=$report. "The only gen on common on every cluster is the main hit\n";
#		if (-e $orderFile){
#			print "I will draw SVG clusters with the single hits order\n";
#			$report=$report. "I will draw with the single hits order\n";
#			$INPUTS=getDrawInputs($orderFile);
 #       		}
		my $line =`perl -ne \'print if \$\. == 2\' $outname/PrincipalHits `;
 		my $len = map $_, $line =~ /(.)/gs;
		$len--;
		$report=$report."\nAminoacid array size = $len \n\n";
        	}
#_____________________________________________________________________________________________

#print "Now SVG file will be generated with inputs: $INPUTS\n\n";
#	print "3_Draw.pl $rescale $INPUTS $outname";
#print "pausa to continue\n";
#my $stdin=<STDIN>;
#	system("3_Draw.pl $rescale $INPUTS $outname");

#print "SVG  file generated\n\n";
#`mv $outname/Contextos.svg $outname/$outname\.svg`;

open (REPORTE, ">$outname/$outname\_Report") or die "Couldn't open reportfile $!";
print REPORTE $report;
close REPORTE;

print "Cleaning temporary files\n";
cleanFiles($outname);

print "Done\n";
print "Have a nice day\n\n";
exit;
######################################################################
######################################################################
###   Sub  Rutinas (llamadas a los distintos pasos del script
#######################################################################
#######################################################################
sub specialCluster{
	my $special_org=shift;
	my @CLUSTERS=qx/ls $outname\/$special_org\_*.input/;
	my $specialCluster="";
	my $score=0;
	foreach my $cluster (@CLUSTERS){
		chomp $cluster;
		#print "I will open #$cluster#\n";
		open (FILE, $cluster) or die "Couldn't open $outname/$cluster\n"; 
		my $firstLine = <FILE>; 
		chomp $firstLine;
		close FILE;
		#print "Primera linea $firstLine\n";
		my @sp=split(/\t/,$firstLine);
			#print "Score $sp[7]\n";
			#print "6 $sp[6] 7 $sp[7]\n";
			if ($score<=$sp[7]){
				$specialCluster=$cluster;
				}
		}
	$specialCluster=~s/\.input//;
	return $specialCluster;
}
#__________________________________________________________________________
sub cleanFiles{
    #    `rm *.lista`;
        `rm $outname/lista.*`;
        `rm $outname/*.input`;
        if (-e "$outname/*.input2"){`rm $outname/*.input2`;}
        `rm $outname/*.input2`;
        `rm $outname/Core`;
        `rm $outname/PrincipalHits`;
        `rm $outname/PrincipalHits.muscle`;
        `rm $outname/PrincipalHits.muscle-gb`;
        `rm $outname/PrincipalHits.muscle-gb.htm`;
        `rm $outname/*.order`;
        `rm $outname/Core0`;
        `rm -r $outname/OUTSTAR`;
        `rm -r $outname/MINI`;
        `rm -r $outname/*.stockholm`;
        `rm -r $outname/*.faa`;
        `rm -r $outname/*.blast`;
#        `rm -r $outname/*.txt`;
        }
#_____________________________________________________________________________________

sub getDrawInputs{
	my $file=shift;
	my $INPUTS="";
	open (NAMES,$file) or die "Couldnt open $file $!";

	foreach my $line (<NAMES>){
		chomp $line;
		my @spt=split(/_org_|_peg_/,$line);
		$INPUTS.=$spt[2]."_".$spt[1]."\.input,";
		#print "$INPUTS\n";
		}
		my $INPUT=chop($INPUTS);
		#print "!$INPUTS!\n";
	#obtener el numero de organismos
	#pasarselo al script 2.Draw.pl
	return $INPUTS;
	}
#_________________________________________________________________________________

sub printVariables{
	if ($verbose){
		print "set_name $setname\n";
		print "e_core $e_core\n ";
		print "blast $blast\n";
		print "list $lista \n";
		print "number $num \n";
		print "dir $dir \n";
		print  "verbose  $verbose\n";
		print "name folder $name\n";
		#print "bitscore $bitscore\n";
		#print "e_value $e_value\n";
		#print "cluster radio $cluster_radio \n";
		#print "rescale $rescale \n";
		}
	}
######################################################################
sub get_number{
        my $list=shift;
        my $list_all=shift;
        my $rast_ids=shift;
        my $NUM;
        #if ($verbose){print "list #$list# total list: $list_all\n";}

        if ($list eq ""){
                $NUM = `wc -l < $rast_ids`;
                chomp $NUM;
                $NUM=int($NUM);
                if ($verbose ){print "Every genome on data base would explored\n"; }
                }
        else {
                if($list_all =~ /not/){
                        print "Your genome list must be numbers separated by , \n you can select intervals using ':'\n Example: 2,3,4:7,9 means 2,3,4,5,6,7,9";
                        print"$list_all\n";
                        exit;
                        }
                else{
                        if ($verbose){print "You will explore genomes $list_all\n";}
                        my @st=split(",",$list_all);
                        $NUM=scalar @st;
                        }
                }
                print "You will explore $NUM genomes\n";
                return $NUM;
        }

#________________________________________________________________________

sub get_lista{
        my $list=shift;
        my $verbose=shift;
        my $rast_ids=shift;
        my $result;
        my $bool=1;
        my @all;

        if ($list eq ""){
                print "\nAll genomes would be procesed\n";
                if (-e $rast_ids){
                        @all=`cut -f1 $rast_ids`;
                        for my $genome (@all){chomp $genome;}
                        @all = grep { $_ ne "" } @all;
                        $result=join(',',@all);
                        if($verbose) { for my $genome (@all){print  "#$genome#\t";}}
                        }
                else {
                        print "$rast_ids file is needed\n$!";
                        exit;
                        }
                }
        else{
                my @split_list= split(",",$list);
                #if ($verbose){print "your list is $list\n";}
                foreach my $st (@split_list){
                        #if($verbose){print "With elements: #$st#\n";}
                        if($st=~/\:/){  ## If st is a range of numbers
                                #if ($verbose){print "Range $st on list\n";}                            
                                my @range=split(/\:/,$st);
                                my $init=$range[0];
                                my $end=$range[1];
                                my $bool_init=0;
                                my $bool_end=0;

                                if($init=~/^\d+\z$/) {$bool_init=1;}
				if($end=~/^\d+\z$/) {$bool_end=1;}

                                if ($init > $end){
                                        print "ERROR: You selected interval $init:$end, at an interval, you must be sure that initial number is lower than end number\n";
                                        exit;
                                        }

                                if($bool_init==1 and $bool_end==1 ){
                                        for (my $element=$init;$element<=$end;$element++){
                                                push(@all,$element);
        #                                       if ($verbose){print "Adding element $element to list\n";}

                                                }
                                        }
                                else{
                                        $bool=$bool_init*$bool_end;
                                        }
                                }
                        else{
                                ## if st is a number
                                my $bool_st=0;
                                if($st=~/^\d+\z$/) {$bool_st=1;}
                                if ($bool_st==1){
                                        push(@all,$st);
                #                       if ($verbose){print "Adding element $st to list\n";}
                                        }
                                else{$bool=$bool*$bool_st;}
                                }
                        }
                if($bool==1){
                        my @sorted_numbers = sort { $a <=> $b } @all;
                            my @unique = do { my %seen; grep { !$seen{$_}++ } @sorted_numbers };
                        $result=join(',',@unique)}
                        else{$result="This is not a list number acepted format, please use only , and : to separate numbers\n";}

                #if ($verbose) {print "Whole list #$result#\n";}
        }
        return $result;
        }

