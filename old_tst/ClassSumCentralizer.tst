gap> tester := function(rep)
> local G, cent_basis, cc;
> G := Source(rep.rep);
> cc := ConjugacyClasses(G);
> # need to orthonormalize (bad in practice, coefficients blow up)
> cent_basis := OrthonormalBasis@RepnDecomp(rep.centralizer_basis);
> return ForAll(cc, cl -> Sum(cl, s -> Image(rep.rep, s)) = ClassSumCentralizer(rep.rep, cl, cent_basis));
> end;;
gap> true; # TestMany@RepnDecomp(tester, 5); # TODO: make this not take forever
true