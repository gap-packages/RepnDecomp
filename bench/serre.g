Read("common.g");

BenchForSmallGroups(function(r,i) return BlockDiagonalRepresentation(r); end, "serre.txt");
