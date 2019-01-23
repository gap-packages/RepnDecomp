Read("common.g");

BenchCanonicalSummandForSmallGroups(DecomposeCanonicalSummandFast@RepnDecomp, "canonical_fast.txt");
BenchCanonicalSummandForSmallGroups(DecomposeCanonicalSummand@RepnDecomp, "canonical_serre.txt");

BenchRepForSmallGroups(BlockDiagonalRepresentationFast, "fast.txt");
BenchRepForSmallGroups(function(r,i) return BlockDiagonalRepresentation(r); end, "serre.txt");

BenchDegreeForSmallGroups(100, BlockDiagonalRepresentationFast, "fast_deg.txt");
BenchDegreeForSmallGroups(100, function(r,i) return BlockDiagonalRepresentation(r); end, "serre_deg.txt");
