#!/usr/bin/perl -w

use strict;

use constant FALSE => 0;
use constant TRUE => 1;

my $asymm = shift;

my $inCoefficients = FALSE;
while (<>) {
  next unless m/^-----/ || $inCoefficients;
  last if m/^Posterior/;
  $inCoefficients = TRUE;
  next unless m/(alpha|beta)/;
  my ($coefficient, $lo, $hi) = split;
  $hi =~ s/\*//;
  next unless keep($lo, $hi, $asymm);
  print "$coefficient: ($lo, $hi)";
  if (($lo > 0.0) || ($hi < 0.0)) {
    print "*";
  }
  print "\n";
}

sub keep {
  my $lo = shift;
  my $hi = shift;
  my $asymm = shift;
  my $retval;
  if ($hi <= 0.0) {
    $retval = TRUE;
  } elsif ($lo >= 0.0) {
    $retval = TRUE;
  } elsif ($asymm > 0.0) {
    my $loTest = min(abs($lo), abs($hi));
    my $hiTest = max(abs($lo), abs($hi));
    $retval = ($hiTest/$loTest >= $asymm) ? TRUE : FALSE;
  }
  return $retval;
}

sub min {
  my $x = shift;
  my $y = shift;
  return ($x < $y) ? $x : $y;
}

sub max {
  my $x = shift;
  my $y = shift;
  return ($x > $y) ? $x : $y;
}
