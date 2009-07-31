#!/usr/bin/perl
# Used to see how we did.
use strict;
use warnings;

my $reposfile = "../../input/repos.txt";
my $datafile = "../../input/data.txt";
my $langfile = "../../input/lang.txt";

unless ($ARGV[0]) {
    print STDERR "Usage: $0 USER:repo1,repo2,...\n";
    print STDERR "E.g.: $0 1654:654,226,616,58,76,8,84,70,650,81\n";
    exit 1;
}

my ($user, $repos) = split(/:/, $ARGV[0]);

open(my $f, $datafile);
my @x = ();
while (<$f>) {
    push(@x, $1) if (/^$user:(\d+)/);
}
close($f);


my $j = join("|", @x);
my $r = qr/^($j):/o;
my %x = ();
open($f, $langfile);
while (<$f>) {
    chomp;
    if ($_ =~ $r) {
        my ($x, $y) = split(/:/, $_);
        $x{$x} = $y;
    }
}
close($f);

print "USER REPOS\n==========\n";
open($f, $reposfile);
while (<$f>) {
    chomp;
    if ($_ =~ $r) {
        print "${_}::" . ($x{$1} || '')."\n";
    }
}
close($f);

print "\n\nRECOMMENDED REPOS\n===================\n";
@x = ();
foreach (split(/,/, $repos)) {
    my ($y) = split(/;/);
    push(@x, $y);
}
$repos = join("|", @x);
$r = qr/^($repos):/;

%x = ();
open($f, $langfile);
while (<$f>) {
    chomp;
    if ($_ =~ $r) {
        my ($x, $y) = split(/:/, $_);
        next unless ($y);
        $x{$x} = $y;
    }
}
close($f);

open($f, $reposfile);
my $i = 0;
while (<$f>) {
    chomp;
    if ($_ =~ $r) {
        print "${_}::" . ($x{$1} || '')."\n";
        $i++;
        last if ($i >= 10);
    }

}
close($f);
