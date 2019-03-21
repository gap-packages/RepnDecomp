gap> tester := function(rep)
> local conds, decomp;
> conds := [];
> decomp := CanonicalDecomposition(rep.rep);
> # if we restrict to a canonical summand, that restricted rep
> # should only have 1 canonical summand itself
> Add(conds, ForAll(decomp, V -> Length(CanonicalDecomposition(RestrictRep@RepnDecomp(rep.rep, Basis(V)))) = 1));
> # if we restrict to canonical summands then take the direct sum
> # we should get the original representation
> Add(conds, AreRepsIsomorphic(rep.rep, DirectSumOfRepresentations(List(decomp, V -> RestrictRep@RepnDecomp(rep.rep, Basis(V))))));
> return ForAll(conds, x->x);
> end;;
gap> TestMany@RepnDecomp(tester, 3);
true