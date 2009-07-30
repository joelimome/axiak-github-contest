#!/usr/bin/perl
my $num = 0;
while (<>) {
	my @x = split(/,/);
	$num += scalar @x;
}
print "$num\n";
