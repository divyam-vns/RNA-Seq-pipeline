#!/usr/bin/perl
use strict;

($#ARGV==0) or die "Usage: $0 annot_gtf < in_results_file > out_results_file\n";

my $annotFile = shift;

my %GeneName; # indexed by gene_id
open(F, "<$annotFile") or die "Fatal Error: Could not open GTF gene annotation file $annotFile.\n";
while (<F>) {
   next if /^#/;

   /gene_name (\S+);/ or die "Fatal error: Incorrect format for GTF line (gene_name). $_";
   my $gene_name = $1;
   if ($gene_name =~/^\"(\S+)\"/) { $gene_name = $1; }

   /gene_id (\S+);/ or die "Fatal error: Incorrect format for GTF line (gene_id). $_";
   my $gene_id = $1;
   if ($gene_id =~/^\"(\S+)\"/) { $gene_id = $1; }

   $GeneName{$gene_id} = $gene_name;
}
close(F);

while (<>) {
  chomp;
  my @elems = split ',', $_;
  my $x = $elems[0]; 
  if ($x=~/^\"(\S+)\"$/) { $x = $1; }
  if (defined($GeneName{$x})) { print "\"",$GeneName{$x}, "\","; }
  else { print "\"\","; }
  print $_, "\n";
}

exit;
