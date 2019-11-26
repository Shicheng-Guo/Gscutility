
# /gpfs/home/guosa/hpc/project/TCGA/mutation

my @file=glob("*.txt");
my $data;
my $symbol;
foreach my $file(@file){
open F, $file;
while(<F>){
next if !/Missense|Nonsense/;
my ($id,$symbol)=split/\s+/;
$data{$id}{$symbol}++;
$symbol{$symbol}=$symbol;
}
}

foreach my $gene (sort keys %symbol){
push @genes,$gene;
}
my $genes=join("\t",@genes);

open OUT,">maf2matrix.tab";
print OUT "\t$genes\n";

foreach my $id(sort keys %data){
    print OUT "$id";
    foreach my $symbol(sort keys %symbol){
    if(defined $data{$id}{$symbol}){
    print OUT "\t$data{$id}{$symbol}";
	}else{
	print OUT "\t0";
	}
	}
	print OUT "\n";
}
