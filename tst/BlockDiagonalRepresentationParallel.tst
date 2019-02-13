gap> tester := rep -> TestBlockDiagonalRepresentation@RepnDecomp(rep, BlockDiagonalRepresentationParallel(rep.rep, 2));;
gap> TestMany@RepnDecomp(tester, 2);
true