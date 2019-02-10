# Correctness checks, so we can see if our decompositions, centralizer
# bases, etc are correct

# Generate a random representation of a random group
RandomRepresentation@ := function()
    # smallgrp has groups of order at most 2000 except 1024
    size := Random(Concatenation([1..1023], [1025..2000]));
    id := Random([1..NrSmallGroups(size)]);
    G := SmallGroup(size, id);
    irreps := IrreducibleRepresentations(G);

    # we pick 2 of the same irrep and 1 other for some variety
    irrep1 := Random(irreps);
    irrep2 := Random(irreps);

    rho := DirectSumRepList([irrep1, irrep1, irrep2]);
    n := DegreeOfRepresentation(rho);

    # scramble the basis to avoid it being too easy
    A := RandomInvertibleMat(n);
    rho := ComposeHomFunction(rho, x -> A^-1 * x * A);

    # We know some things about rho without needing to compute them
    # the long way. We return them for use in testing.

    # centralizer basis
    centralizer_basis := [rec(dimension := DegreeOfRepresentation(irrep1),
                              nblocks := 2),
                          rec(dimension := DegreeOfRepresentation(irrep2),
                              nblocks := 1)];
    centralizer_basis := List(SizesToBlocks@(centralizer_basis), BlockDiagonalMatrix);

    # put them in the scrambled basis
    centralizer_basis := List(centralizer_basis, M -> A^-1 * M * A);

    isomorphism_type := [rec(rep := irrep1, m := 2),
                         rec(rep := irrep2, m := 1)];

    return rec(rep := rho,
               isomorphism_type := isomorphism_type,
               centralizer_basis := centralizer_basis,
               G := G);
end;

# checks if the given centralizer basis is correct
TestCentralizerBasis@ := function(rep, cent_basis)
    local correct_basis, C1, C2;

    correct_basis := rep.centralizer_basis;

    C1 := VectorSpace(Cyclotomics, correct_basis);
    C2 := VectorSpace(Cyclotomics, cent_basis);

    return C1 = C2;
end;

# checks if decomp is a decomposition into G-invariant subspaces
TestInvariantDecomposition@ := function(rep, decomp)
    local rho, G;

    rho := rep.rep;
    G := Source(rho);

    return ForAll(decomp, V -> IsGInvariant(G, rho, V));
end;

# checks some necessary conditions for decomp to be the full
# irreducible (collected) decomposition of rho
TestIrreducibleDecomposition@ := function(rep, decomp)
    local rho, G, conds;

    rho := rep.rep;

    conds := [];

    # must be an invariant decomposition
    Add(conds, TestInvariantDecomposition@(rep, Flat(decomp)));

    # we know the degrees and isomorphism type, so we know the
    # dimensions and the lengths of lists we should have

    # same number of irrep types found
    Add(conds, Length(decomp) = Length(rep.isomorphism_type));

    # same number of irreps in total, with same dimension
    Add(conds, SortedList(List(Flat(decomp), Dimension)) =
               SortedList(Flat(List(rep.isomorphism_type,
                                    t -> Replicate@(DegreeOfRepresentation(t.rep),
                                                    t.m)))));


    # we can't really verify correctness fully without doing the
    # calculation ourselves I think

    # TODO: add some more necessary conditions

    return ForAll(conds, x->x);
end;

# same as irreducible decomp, these are only necessary conditions
# for correctness
TestCanonicalDecomposition@ := function(rep, decomp)
    local conds;

    conds := [];

    # needs to be a G-invariant decomp
    Add(conds, TestInvariantDecomposition@(rep, decomp));

    # same number of irrep types found
    Add(conds, Length(decomp) = Length(rep.isomorphism_type));

    # right dimensions
    Add(conds, SortedList(List(decomp, Dimension)) =
               SortedList(List(rep.isomorphism_type,
                               t -> DegreeOfRepresentation(t.rep)*t.m)));

    return ForAll(conds, x->x);
end;
