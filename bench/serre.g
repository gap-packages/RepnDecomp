Read("common.g");

BenchRepForSmallGroups(function(r,i) return BlockDiagonalRepresentation(r); end, "serre.txt");
