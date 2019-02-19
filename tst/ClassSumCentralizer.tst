gap> tester := function(rep)
> local G, cent_basis, cc, A, diag_rep;
> G := Source(rep.rep);
> cc := ConjugacyClasses(G);
> # need to change to nice basis since we only have an ortho basis for C there
> A := TransposedMat(rep.candidate_nice_basis);
> # need to normalize (already ortho)
> cent_basis := List(rep.centralizer_basis, x -> A * x * A^-1);
> cent_basis := List(cent_basis, m -> 1/Sqrt(Trace(m * ConjugateTranspose@RepnDecomp(m))) * m);
> diag_rep := ComposeHomFunction(rep.rep, x -> A * x * A^-1);
> return ForAll(cc, cl -> Sum(cl, s -> Image(diag_rep, s)) = ClassSumCentralizer(diag_rep, cl, cent_basis));
> end;;
gap> TestMany@RepnDecomp(tester, 5);
true