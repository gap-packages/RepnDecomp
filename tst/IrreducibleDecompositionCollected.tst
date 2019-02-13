gap> tester := rep -> TestIrreducibleDecomposition@RepnDecomp(rep, List(IrreducibleDecompositionCollected(rep.rep).decomp, l -> List(l, r -> r.space)));;
gap> TestMany@RepnDecomp(tester, 5);
true