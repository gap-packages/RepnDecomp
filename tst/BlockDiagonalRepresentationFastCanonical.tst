gap> tester := rep -> TestBlockDiagonalRepresentation@RepnDecomp(rep, BlockDiagonalRepresentationFastCanonical(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true