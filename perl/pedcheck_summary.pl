#!perl
use strict;
use warnings;
while(<>){
if(/# Pedigree (\d+) removed\.)
print "$1\n"
}