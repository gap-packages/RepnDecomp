gap> tester := function(rep)
> local A;
> A := TransposedMat(BlockDiagonalBasisOfRepresentation(rep.rep));
> return TestCentralizerBasis@RepnDecomp(rep, List(CentralizerOfRepresentation(rep.rep), M -> A * M * A^-1));
> end;;
gap> TestMany@RepnDecomp(tester, 5);
true