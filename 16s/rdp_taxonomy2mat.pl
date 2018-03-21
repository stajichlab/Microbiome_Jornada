#!/usr/bin/perl
use strict;
use warnings;

my $header = <>;
my @levels = qw(k p c o f g s);
print join("\t","",qw(Kingdom Phylum Class Order Family Genus Species)),"\n";

while(<> ) {
    s/\"//g;
    chomp;
    my ($id,@taxonomystrings) = split(/\t/,$_);
    my %d;
    for my $w ( @taxonomystrings ) {
	(undef,$w) = split(/;/,$w);
	my @plevels = split(/,/,$w);

	for my $l ( @plevels ) {
	    my ($name,$val) = split(/:/,$l);
	    $d{$name} = $val if ! exists $d{$name};
	}
	last
    }
    my $last = 'Unknown';
    print join("\t", $id, map { my $r = $d{$_} || 'NA';
	       $last = $r;
	       $r } @levels), "\n";
}
