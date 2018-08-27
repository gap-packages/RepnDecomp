#!/usr/bin/perl
use warnings;
use strict;

my $source_code_file = $ARGV[0];


# read the code blocks into a map with the names
open(my $fh, '<', $source_code_file);

my %code_blocks_map = ();

my $current_code_block = "";
my $current_block_name = "NONE";

while(<$fh>) {
  if (/^#%%(.+)$/) {
    $code_blocks_map{$current_block_name} = $current_code_block;
    $current_code_block = "";
    $current_block_name = $1;
  } else {
    $current_code_block .= $_;
  }
}

$code_blocks_map{$current_block_name} = $current_code_block;

# now read in the input source and wherever we see a tag, put the
# corresponding block instead
while(<STDIN>) {
  if (/^%%#(.+)$/) {
    print $code_blocks_map{$1};
  } else {
    print;
  }
}
