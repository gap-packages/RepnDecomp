gap> tester := rep -> TestBlockDiagonalRepresentation@RepnDecomp(rep, REPN_ComputeUsingMyMethodCanonical(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true