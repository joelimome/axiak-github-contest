#!/usr/bin/perl
while (<>) {
    my ($repo) = split(/:/);
    print "$repo:\n";
}
