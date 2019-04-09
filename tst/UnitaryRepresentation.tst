gap> deunitarise := function(rho)
> local A;
> A := RandomInvertibleMat(DegreeOfRepresentation(rho));
> return ComposeHomFunction(rho, m -> A^-1 * m * A);
> end;;
gap> # need big group otherwise it'll be too trivial
gap> G := SmallGroup(RandomGroup@RepnDecomp(rec(lo:=100, hi:=200)));;
gap> irreps := IrreducibleRepresentations(G);;
gap> # this is because the reps probably already are unitary
gap> irreps := List(irreps, deunitarise);;
gap> # try to deunitarise 3 irreps (there might be too many to test them all)
gap> irreps := List([1..3], x -> Random(irreps));;
gap> ForAll(irreps, i -> IsUnitaryRepresentation(UnitaryRepresentation(i).unitary_rep));
true