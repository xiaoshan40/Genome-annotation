#!/usr/bin/perl -w
use strict;

my $timer;

open(exonerate,$ARGV[0]) or die $!;
while(<exonerate>){
  chomp();
  if($_=~/(\S+)\t(\S+)\t(gene)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/){
    $timer++;
    print "$1\t$2\tmRNA\t$4\t$5\t"."\."."\t$7\t$8\tID=exonerate$timer;\n";
  }elsif($_=~/(\S+)\t(\S+)\t(exon)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/){
    print "$1\t$2\tCDS\t$4\t$5\t"."\."."\t$7\t$8\tParent=exonerate$timer;\n";
  }else{
    next;
  }
}
close exonerate;
