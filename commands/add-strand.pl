#!/usr/bin/perl
use strict;

($#ARGV==0 || $#ARGV==1) or die "Usage: $0 annot_gtf [-r] < in_gtf_file > out_gtf_file\n" .
                   "          -r     remove no-strand genes and transcripts\n";

my $annotFile = shift;
my $rFlag = shift;

my %Strand; # indexed by gene_name
open(F, "<$annotFile") or die "Fatal Error: Could not open GTF gene annotation file $annotFile.\n";
while (<F>) {
   next if /^#/;

   /\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t/ or die "Fatal Error: Incorrect format for GTF line (strand). $_";
   my $strand = $1;
   /gene_name (\S+);/ or die "Fatal error: Incorrect format for GTF line (gene_name). $_";
   my $gene_name = $1;
   if ($gene_name =~/^\"(\S+)\"/) { $gene_name = $1; }
   $Strand{$gene_name} = $strand;
}
close(F);

while (<>) {
   if (m/^#/) {
      print $_; next;
   }

   chomp;
   /(\S+\t\S+\t\S+\t\S+\t\S+\t\S+)\t(\S+)\t(.*)$/ or die "Fatal Error: Incorrect format for input GTF file. $_"; 
   my ($left,$strand,$right) = ($1,$2,$3);
   if (($strand eq "+") || ($strand eq "-")) {
      print $_, "\n";
   } elsif ($strand eq ".") {
      if ($right=~/gene_name (\S+);/) {
	 my $gene_name = $1;
         if ($gene_name=~/^\"(\S+)\"/) { $gene_name = $1; }
         if (defined($Strand{$gene_name})) {
	    print "$left\t", $Strand{$gene_name}, "\t$right\n";
	 } else {
	    defined($rFlag) or print $_, "\n";
	 }
      } else {
         defined($rFlag) or print $_, "\n";
      }
   } else {
      die "Fatal Error: Incorrect line format (unrecognized strand). $_";
   }
}

exit(0);
      
