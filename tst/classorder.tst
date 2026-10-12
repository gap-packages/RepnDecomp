gap> START_TEST("classorder.tst");

# The classes of CharacterTable(G) need not be ordered like
# ConjugacyClasses(G). Here the group lists classes 2 and 3 of the
# table in the other order.
gap> S4tbl := CharacterTable(SymmetricGroup(4));;
gap> S4reps := List(ConjugacyClasses(S4tbl), Representative){[1,3,2,4,5]};;
gap> S4 := SymmetricGroup(4);;
gap> SetConjugacyClasses(S4, List(S4reps, x -> ConjugacyClass(S4, x)));
gap> # CharacterTableWithStoredGroup wants a table that stores no group
> S4tbl := ConvertToCharacterTable(rec(
>     UnderlyingCharacteristic := 0,
>     OrdersClassRepresentatives := OrdersClassRepresentatives(S4tbl),
>     SizesCentralizers := SizesCentralizers(S4tbl),
>     Irr := List(Irr(S4tbl), ValuesOfClassFunction)));;
gap> S4tbl := CharacterTableWithStoredGroup(S4, S4tbl);;
gap> IsIdenticalObj(S4tbl, CharacterTable(S4));
true
gap> IdentificationOfConjugacyClasses(S4tbl);
[ 1, 3, 2, 4, 5 ]

# The natural permutation character is the sum of the trivial character
# and one of degree 3.
gap> S4proj := ProjIrr@RepnDecomp(CohCfgFromPermGroup@RepnDecomp(S4));;
gap> SortedList(TransposedMat([S4proj.dims, S4proj.blocksizes]));
[ [ 1, 1 ], [ 3, 1 ] ]

# Decompositions do not look at the table when the irreps are given.
gap> S4irreps := IrreducibleRepresentations(S4);;
gap> S4rho := DirectSumOfRepresentations([S4irreps[4], S4irreps[5]]);;
gap> S4dec := REPN_ComputeUsingMyMethod(S4rho : irreps := S4irreps);;
gap> Length(S4dec.centralizer_basis);
2
gap> List(S4dec.decomposition, Vs -> List(Vs, Dimension));
[ [ 3 ], [ 3 ] ]
gap> AreRepsIsomorphic(S4rho, S4dec.diagonal_rep : irreps := S4irreps);
true
gap> List(CanonicalDecomposition(S4rho : irreps := S4irreps), Dimension);
[ 3, 3 ]

# Characters of the table can be passed along with the irreps. Results
# are stored in the representation, so this needs a new one.
gap> S4chars := List(S4irreps, r -> First(Irr(S4),
>         chi -> ForAll(S4reps, g -> g^chi = Trace(Image(r, g)))));;
gap> S4rho := DirectSumOfRepresentations([S4irreps[4], S4irreps[5]]);;
gap> S4withchars := REPN_ComputeUsingMyMethod(S4rho : irreps := S4irreps,
>                                                    irr_chars := S4chars);;
gap> List(S4withchars.decomposition, Vs -> List(Vs, Dimension));
[ [ 3 ], [ 3 ] ]
gap> S4serrechars := REPN_ComputeUsingSerre(S4rho : irreps := S4irreps,
>                                                  irr_chars := S4chars);;
gap> Length(S4serrechars.centralizer_basis);
2

# The degree of an irrep is not read off the first class of
# ConjugacyClasses(G): here that list ends with the identity class.
gap> S4reps := List(ConjugacyClasses(SymmetricGroup(4)), Representative);;
gap> S4reps := Concatenation(S4reps{[2..5]}, S4reps{[1]});;
gap> S4 := SymmetricGroup(4);;
gap> SetConjugacyClasses(S4, List(S4reps, x -> ConjugacyClass(S4, x)));
gap> S4irreps := IrreducibleRepresentations(S4);;
gap> S4rho := DirectSumOfRepresentations([S4irreps[4], S4irreps[5]]);;
gap> S4my := REPN_ComputeUsingMyMethod(S4rho : irreps := S4irreps);;
gap> Length(S4my.centralizer_basis);
2
gap> List(S4my.decomposition, Vs -> List(Vs, Dimension));
[ [ 3 ], [ 3 ] ]
gap> S4serre := REPN_ComputeUsingSerre(S4rho : irreps := S4irreps);;
gap> Length(S4serre.centralizer_basis);
2
gap> List(Filtered(S4serre.decomposition, l -> Length(l) > 0),
>         l -> List(l, r -> Length(r.basis)));
[ [ 3 ], [ 3 ] ]

#
gap> STOP_TEST("classorder.tst", 1);
