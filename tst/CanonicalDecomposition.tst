gap> tester := rep -> TestCanonicalDecomposition@RepnDecomp(rep, CanonicalDecomposition(rep.rep));;
gap> TestMany@RepnDecomp(tester, 5);
true