gap> tester := rep -> TestBlockDiagonalRepresentation@RepnDecomp(rep, BlockDiagonalRepresentationFast(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true