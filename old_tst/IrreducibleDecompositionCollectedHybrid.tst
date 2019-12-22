gap> tester := rep -> TestIrreducibleDecomposition@RepnDecomp(rep, List(IrreducibleDecompositionCollectedHybrid@RepnDecomp(rep.rep).decomp, l -> List(l, r -> VectorSpace(Cyclotomics, r.basis))));;
gap> TestMany@RepnDecomp(tester, 5);
true
