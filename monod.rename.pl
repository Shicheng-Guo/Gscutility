use strict;
my @file=glob("6-*.bam");
my $i;
my $j;
foreach my $file(@file){
        $i++;
        $j++;
        my ($cancer,$type,$id)=split /[.-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file CRC-$type-$id2.bam");
        print "$id\t$id2\n";
}


use strict;
my @file=glob("PC-*.bam");
my $i;
foreach my $file(@file){
        $i++;
        $j++;
        my ($cancer,$type,$id)=split /[.-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file PC-$type-$id2.bam");
        print "$id\t$id2\n";
}


my @file=glob("7-*.hapInfo.txt");
my $i;
foreach my $file(@file){
        $i++;
        $j++;
        my ($cancer,$type,$id)=split /[.-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file LC-$type-$id2.bam");
        print "$id\t$id2\n";
}

my @file=glob("NC-*.bam");
my $i;
foreach my $file(@file){
        $i++;
        $j++;
        my ($cancer,$type,$id)=split /[.-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file NC-$type-$id2.bam");
        print "$id\t$id2\n";
}


my @file=glob("RRBS-6P*.bam");
my $i;
foreach my $file(@file){
        $i++;
        $j++;
        my (undef,$cancer,$type,$id)=split /[.P-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file CRC-P-$id2.bam");
        print "$id\t$id2\n";
}

my @file=glob("RRBS-7P*.bam");
my $i;
foreach my $file(@file){
        $i++;
        $j++;
        my (undef,$cancer,$type,$id)=split /[.P-]/g,$file;
        my $id2 = sprintf("%03d",$id);
        system("cp $file LC-P-$id2.bam");
        print "$id\t$id2\n";
}


