gap> tester := rep -> TestInvariantDecomposition@RepnDecomp(rep, IrreducibleDecomposition(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true
gap> # also with otimes tricks that use more memory
gap> tester := rep -> TestInvariantDecomposition@RepnDecomp(rep, IrreducibleDecomposition(rep.rep : use_kronecker));;
gap> TestMany@RepnDecomp(tester, 2);
true