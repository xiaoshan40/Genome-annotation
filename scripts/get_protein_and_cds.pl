#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my (%Prot_length,%Prot_seq,$PROTEIN_ID);
open my $PROTEIN, "<", $ARGV[0] or die "Cannot find the specified fasta file!";
while (<$PROTEIN>) {
  chomp();
  if ($_=~/>(\S+)/){
    $PROTEIN_ID=$1;
  }else{
    $Prot_length{$PROTEIN_ID}+=length($_);
    $Prot_seq{$PROTEIN_ID}.="$_";
  }
}
close $PROTEIN;

my (%cds_seq,$prot_ID);
open my $CDS, "<", $ARGV[1] or die "Cannot find the specified fasta file!";
while (<$CDS>) {
  chomp();
  if ($_=~/>\S+_cds_(\S+)_\S/){ #>lcl|XM_021335028.1_cds_XP_021190703.1_1 
    $prot_ID=$1;
  }else{
    $cds_seq{$prot_ID}.="$_";
  }
}
close $CDS;

my @list;
open my $GFF, "<", $ARGV[2] or die "Cannot find the specified fasta file!";
while (<$GFF>) {
  chomp();
  if ($_=~/\S+\tCDS\t.*ID\=cds-(\S*?);Parent\=rna-(\S*?);.*gene\=(\S*?);.*product\=(.*?);/){
#...CDS...ID=cds-XP_0215.1;Parent=rna-XM_0210.1;...gene=LOC1103...product=protein 2;...
    push @list,[$3,$1,$Prot_length{$1},$2,$4];
  }
}
close $GFF;

my @sorted_list = sort { $a->[0] cmp $b->[0] || $b->[2] <=> $a->[2] } @list;

my ($protein_result,$cds_result,$annotation_result,%output,$prefix);

GetOptions(
  'prefix=s' => \$prefix,
);

for (0..$#sorted_list){
  if (exists $output{$sorted_list[$_][0]}){
    next;
  }else{
    $output{$sorted_list[$_][0]}=1;
    $protein_result.=">$prefix"."_$sorted_list[$_][1]\n$Prot_seq{$sorted_list[$_][1]}\n";
    $cds_result.=">$prefix"."_$sorted_list[$_][1]\n$cds_seq{$sorted_list[$_][1]}\n";
    $annotation_result.="$prefix"."_$sorted_list[$_][1]\t$sorted_list[$_][4]\n";
  }
}

open OUT1, ">", "$prefix"."\.protein.fasta";
print OUT1 $protein_result;

open OUT2, ">", "$prefix"."\.cds.fasta";
print OUT2 $cds_result;

open OUT3, ">", "$prefix"."\.annotation.txt";
print OUT3 $annotation_result;
