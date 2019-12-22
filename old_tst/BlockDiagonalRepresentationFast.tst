gap> tester := rep -> TestBlockDiagonalRepresentation@RepnDecomp(rep, REPN_ComputeUsingMyMethod(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true