#!/usr/bin/perl
use strict;
use warnings;

# Take an input file and indescriminately generate
# a results.txt compatible file.
my $num_results = 0;
my $lines = 0;
while (<>) {
    chomp;
    my ($user, $rest) = split(/:/o, $_, 2);
    my @repos = split(/,/o, $rest);
    my @results = ();

    foreach my $repo (@repos) {
        ($repo) = split(/;/o, $repo, 2);
        $num_results ++;
        push(@results, $repo);
        last if (scalar @results >= 10);
    }

    $lines ++;
    print sprintf("%s:%s\n", $user,
                  join(",", @results));
}

print STDERR "Printed $num_results repos in $lines lines.\n";
