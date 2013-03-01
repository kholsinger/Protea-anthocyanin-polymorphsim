#!/usr/bin/perl -w

use strict;

use constant FALSE => 0;
use constant TRUE => 1;

my $asymm = shift;

my $inSlopes = FALSE;
while (<>) {
  next unless m/^Environmental coefficients/ || $inSlopes;
  $inSlopes = TRUE;
  next unless m/^beta/;
  my ($coefficient, $mean, $int) = split;
  my ($lo, $hi) = split /,/, $int;
  $coefficient =~ s/://;
  $lo =~ s/\(//;
  $hi =~ s/\)\**//;
  next unless keep($lo, $hi, $asymm);
  print "$coefficient: $mean ($lo, $hi)\n";
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
  } else {
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
