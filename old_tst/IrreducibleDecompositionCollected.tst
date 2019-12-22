gap> tester := rep -> TestIrreducibleDecomposition@RepnDecomp(rep, IrreducibleDecompositionCollected(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true
gap> # and with kronecker tricks
gap> tester := rep -> TestIrreducibleDecomposition@RepnDecomp(rep, IrreducibleDecompositionCollected(rep.rep : use_kronecker), l -> List(l, r -> r.space));;
gap> TestMany@RepnDecomp(tester, 2);
true