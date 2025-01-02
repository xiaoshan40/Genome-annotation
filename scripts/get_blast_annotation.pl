#!/usr/bin/perl -w
use strict;

my $last_gene='';
my (@gene_list,@hit_list);
my %link;

open my $BLAST_OUTPUT, "<", $ARGV[0] or die "Cannot find the specified fasta file!";
while (<$BLAST_OUTPUT>) {
###Gqin00002  XP_021184173.1  75.2  214  53  0  37  250  37  250  6.1e-76  293.9
###Gqin00002  XP_026733160.1  75.1  217  54  0  34  250  34  250  2.3e-75  292.0
###Gqin00003  XP_026487295.1  64.4  907  290  12  24  915  1  889  0.0e+00  1144.0
  chomp();
  my @temp = split /\t/, $_;
  if ($temp[0] ne $last_gene){
    push @gene_list,$temp[0];
    push @hit_list,$temp[1];
    $link{$temp[0]}=$temp[1];
    $last_gene=$temp[0];
  }else{next;}
  }
close $BLAST_OUTPUT;

my %count;
my @hit_list_uniq = grep{++$count{$_}<2}@hit_list;
my $out1=join("\n",@hit_list_uniq);
print "$out1\n";

my %Annotation;
open my $NR, "<", $ARGV[1] or die "Cannot find the specified fasta file!";
while (<$NR>) {
  chomp();
  if ($_=~/^>(\S+)\s+(\S.*?\])/ ){
    if (grep /$1/, @hit_list_uniq ){
      $Annotation{$1}=$2;
    }else{next;}
  }else{next;}
}
close $NR;

my $result;
foreach (@gene_list){
  print "$_\t$Annotation{$link{$_}}\n"};
