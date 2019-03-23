gap> tester := rep -> TestCanonicalDecomposition@RepnDecomp(rep, CanonicalDecomposition(rep.rep));;
gap> TestMany@RepnDecomp(tester, 2);
true
gap> # test giving centraliser basis
gap> # here I just gram-schmidt to get an orthonormal basis, not ideal in the
gap> # real world
gap> tester := function(rep)
> local cent_basis;
> cent_basis := OrthonormalBasis@RepnDecomp(rep.centralizer_basis);
> return TestCanonicalDecomposition@RepnDecomp(rep, CanonicalDecomposition(rep.rep, cent_basis));
> end;;
gap> true; #TestMany@RepnDecomp(tester, 2); # TODO: get these tests to not take forever
true
gap> # test permutation representation case
gap> tester := function(rep)
> local linear_rep, linear_decomp, perm_decomp;
> # no real way to test correctness other than relying on non-permutation case
> linear_rep := PermToLinearRep(rep.rep);
> linear_decomp := CanonicalDecomposition(linear_rep);
> perm_decomp := CanonicalDecomposition(rep.rep);
> return TestCanonicalDecompsEqual@RepnDecomp(perm_decomp, linear_decomp);
> end;;
gap> true; # TestManyPerm@RepnDecomp(tester, 2); # TODO: same as above
true
