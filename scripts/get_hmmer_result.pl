#!/usr/bin/perl -w
use strict;

my %Gene_domain;
my @Gene_list;

open my $HMMSCAN, "<", $ARGV[0] or die "Cannot find the specified fasta file!";

while (<$HMMSCAN>) {
  chomp();
  if ($_=~/^#/ ){
    next;
  }else{
    my @temp=split(/\s+/);
    $Gene_domain{$temp[3]}.="$temp[0]\t";
    push @Gene_list,$temp[3];
  }
}
close $HMMSCAN;

my %count;
my @Gene_list_uniq = grep{++$count{$_}<2}@Gene_list;

my $result;
foreach(@Gene_list_uniq){
  my @domain_list = split /\t/, $Gene_domain{$_};
  my %domain_count;
  my @domain_list_uniq = grep{++$domain_count{$_}<2}@domain_list;
  my $out=join("\;",@domain_list_uniq);
  $result.="$_\t$out\n";
  }

open OUT1, ">", "Gqin.OGS.pfam.txt";
print OUT1 $result;
