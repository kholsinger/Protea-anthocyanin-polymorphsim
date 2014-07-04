#!/usr/bin/perl -w

use strict;

use constant FALSE => 0;
use constant TRUE => 1;

open(INFILE, "<results.txt") or
  die "Could not open results.txt: $!";
my $analysis;
my %tabulate;
my $ct = -1;
while (<INFILE>) {
  if (/^power-analysis\s*$/) {
    $analysis = "Categorical";
    $ct += 1;
  } elsif (/^power-analysis-logistic\s*$/) {
    $analysis = "Logistic";
  }
  if (/^\s+beta/) {
    chomp;
    my $coeff = $_;
    $coeff =~ s/^\s*beta\.(.*)\[.+$/$1/;
    $tabulate{$coeff}{$analysis}{$ct} = 1;
  }
}

my %sum;
foreach my $beta (sort keys %tabulate) {
  $sum{$beta}{"Categorical"} = 0;
  $sum{$beta}{"Logistic"} = 0;
  foreach my $analysis (sort keys %{$tabulate{$beta}}) {
    foreach my $ct (sort keys %{$tabulate{$beta}{$analysis}}) {
      $sum{$beta}{$analysis} += $tabulate{$beta}{$analysis}{$ct}
    }
  }
}

$ct += 1;   # since counting started at 0
print "Total number of simulations: $ct\n\n";

my ($coeff, $categorical, $logistic);
foreach $coeff (sort keys %sum) {
  $categorical = $sum{$coeff}{"Categorical"};
  $logistic = $sum{$coeff}{"Logistic"};
  write;
}

format STDOUT_TOP =
                    Analysis
Coefficient  Categorical  Logistic
-----------  -----------  --------
.

format STDOUT =
@>>>>>>>>>>  @##########  @#######
$coeff,      $categorical, $logistic
.
