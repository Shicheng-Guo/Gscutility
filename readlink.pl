#!/usr/bin/perl
use Cwd 'abs_path';
if (@ARGV < 1 ) {
  print "Usage: $0 dir\n";
  exit;
}
$realfilepath = abs_path("$ARGV[0]");
print "$realfilepath\n";
