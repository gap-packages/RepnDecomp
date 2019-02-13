gap> tester := rep -> TestInvariantDecomposition@RepnDecomp(rep, IrreducibleDecomposition(rep.rep));;
gap> TestMany@RepnDecomp(tester, 5);
true