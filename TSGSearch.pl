use strict;
my $dir="/home/gsc/houston/upload/tsg";

my($abstract,$aliasfile)=@ARGV;

sub textming_extract{
my $file=shift @_;
my $outfile1=">$file.out1.txt";
my $outfile2=">$file.out2.txt";
my $outfile3=">$file.out3.txt";
chdir $dir or die"can't into $dir $!";
my $n=-1;
my $m=0;
my ($Uper,%TSG,%TSG2,%TSG3,%lastall,%lastall_reverse,$TSGtxt,@list,@list2);
my $species="Human|Drosophila";
my $number="As|II|Some|Several|Many|One|Other|Three|The|This|These|Nine|Two|Whereas|It|Specific|three|Notably|Although|Silencing|Similar|Potential|Putative|Silence|Two|Three";
my $bio="DNA|RNA|CCC";
my $condition="As|B-Cell|Adenomatous|Certain|TSG|Hypermethylation|SCC-RELATED|SILENCING|UV|Lindau|RA-regulated|QPCR";
my $paper="OBJECTIVE|PURPOSE|BACKGROUND|Targeting|Identifying|Imprinted|Tazarotene-induced|Ras-related|XIAP-associated|X-linked|SCC-related|NPC-related|NPC-associated|Multiple|Kruppel-like|Known|Hypermethylation-mediated|Disabled|Different|DNMT1-mediated|Conversely|Candidate|Cancer|CLL-associated|CC-related|C-terminal|N-terminal";
my $cancer="Cancer|LOH|Hippel-Lindau|HCC";
my $cancertype="hepatocarcinoma|liver cancer|hepatocelluar carcinoma|liver tumor";

open IN, $file or die"no such $file or can't open $file";
open OUT1,$outfile1;
while (<IN>){
	if (/$cancertype|./){
    if (/\(([A-Z]([\w-])*\w)\s*\), a tumor suppressor gene/){
      $n+=1;
#      print "$1\n";
      $list[$n]=$1;   
	}	
	elsif(/(\(([A-Z]([\w-])*\w)\))\s*(is|as|,)?\s*(a)?\s*(known|potential|candidate|)?\s*tumor\ssuppressor\sgene/){ #REIC/Dickkopf-3 (Dkk-3), a tumor suppressor gene
       	 $n+=1;
#       print "$2\n";
      $list[$n]=$2;
	}
		elsif(/tumor\ssuppressor\sgene\s*(\w*\s)*(\(([A-Z]([\w-])*\w)\))/){ #REIC/Dickkopf-3 (Dkk-3), a tumor suppressor gene
       	 $n+=1;
#       print "$3\n";
      $list[$n]=$3;
	}
	elsif(/tumor suppressor gene ([A-Z]([\w-])*\w)/){
       	 $n+=1;
#       print "$1\n";
      $list[$n]=$1;
      	}	
	elsif(/tumor suppressor gene\s*\(([A-Z]([\w-])*\w)\)/){
      	 $n+=1;
       print "$1\n";
      $list[$n]=$1;
      
	}	
	elsif(/([A-Z]([\w-])*\w)\s*\((.|\s)*\)\s*(is|as|,)?\s*(a)?\s*(known|potential|candidate|putative)?\s*(\w+)\s*tumor\ssuppressor/){ #PDCD4 (programmed cell death 4), a tumor suppressor gene
       	 $n+=1;
#       print "$1\n";
      $list[$n]=$1;
      
	}	
	elsif(/tumor suppressor ([A-Z]([\w-])*\w)/){
      	 $n+=1;
#       print "$1\n";
      $list[$n]=$1;
      
	}
    elsif(/([A-Z][A-Z]([\w-])*\w).*(functions|acts|,) as a tumor\ssuppressor/){
       	 $n+=1;
#       print "$1\n";
      $list[$n]=$1;
	}	
    elsif(/([A-Z]([\w-])*\w)\s*\(a tumor suppressor gene\)/){
       	 $n+=1;
#       print "$1\n";
      $list[$n]=$1;
      
	}		
	elsif(/([A-Z]([\w-\.])*\w)\s*(is|as|,)?\s*(a)?\s*(known|potential|candidate|putative)?\s*tumor\ssuppressor/){  #REIC/Dickkopf-3 (Dkk-3), a tumor suppressor gene
       	 $n+=1;
#       print "$1\n"; 
      $list[$n]=$1;
	}	
}
	}	
    print OUT1"tumor suppressor gene\t frequency\n";   
    
foreach (@list){       
if (/($number|$bio|$paper|$condition|$cancer|$species)/) {	
}
else {	   	
    print  "$_\t";
    my $Uper=uc($_);
    print OUT1 "$Uper\t";
   	$TSG2{$Uper}++; 
   	print OUT1 "$TSG2{$Uper}\n";
  # print  "$TSG2{$Uper}\n";
    } 	
  }
 close IN;
 close OUT1;
 my @uniquegenesymbol;
 open OUT2,$outfile2;
 foreach my $TSGU(sort keys %TSG2){
 	my $value=$TSG2{$TSGU};
 	#print "$TSGU\t";
 	push @uniquegenesymbol,$TSGU;
  # print "$value\n";
 	print OUT2 "$TSGU\t";
 	print OUT2 "$value\n";
 }
    return(@uniquegenesymbol); 
 }


sub alias_unique{
	my @gene0=@_;	
	my ($geneout,@geneout);
	open F,$aliasfile;
	my (%gene);
	while(<F>){
	my @line=split /\t/;
	foreach my $i (0..$#line){
	$gene{$line[$i]}=$line[0];
	}
	}
	
	foreach my $gene0(@gene0){
	$geneout=$gene{$gene0};
	push @geneout, $geneout;
	}
	return(@geneout);
}
 
my @sortunique_genelist= textming_extract($abstract);
my @aliasunique_genelist= alias_unique(@sortunique_genelist) ;

foreach my $uniquegen(@aliasunique_genelist){
	print "$uniquegen\n" unless ($uniquegen eq "") ;	
} 




