gap> tester := function(rep)
> local A;
> A := BlockDiagonalBasis(rep.rep);
> A := TransposedMat(A);
> return TestCentralizerBasis@RepnDecomp(rep, List(RepresentationCentralizer(rep.rep), M -> A * M * A^-1));
> end;;
gap> TestMany@RepnDecomp(tester, 5);
true