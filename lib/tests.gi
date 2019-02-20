# Correctness checks, so we can see if our decompositions, centralizer
# bases, etc are correct

# Generate a random representation of a random group
RandomRepresentation@ := function()
    local size, id, G, irreps, irrep1, irrep2, rho, n, A, centralizer_basis, isomorphism_type;

    # smallgrp has groups of order at most 2000 except 1024
    #size := Random(Concatenation([1..1023], [1025..2000]));
    size := Random([1..100]); # don't want tests to take forever
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

    # if the irreps ended up the same, we have to correct the info
    if irrep1 = irrep2 then
        centralizer_basis := [rec(dimension := DegreeOfRepresentation(irrep1),
                                  nblocks := 3)];
        isomorphism_type := [rec(rep := irrep1, m := 3)];
    else
        centralizer_basis := [rec(dimension := DegreeOfRepresentation(irrep1),
                                  nblocks := 2),
                              rec(dimension := DegreeOfRepresentation(irrep2),
                                  nblocks := 1)];

        isomorphism_type := [rec(rep := irrep1, m := 2),
                             rec(rep := irrep2, m := 1)];
    fi;


    centralizer_basis := List(SizesToBlocks@(centralizer_basis), BlockDiagonalMatrix);

    # put them in the scrambled basis
    centralizer_basis := List(centralizer_basis, M -> A^-1 * M * A);

    return rec(rep := rho,
               isomorphism_type := isomorphism_type,
               centralizer_basis := centralizer_basis,
               candidate_nice_basis := TransposedMat(A), # note this is not the unique right answer...
               G := id);
end;

# we return a random (injective!) permutation representation
RandomPermRepresentation@ := function()
    local id, G, degree, H, all_homs, rho;

    id := RandomGroup@();
    G := SmallGroup(id[1], id[2]);

    # This is what we want to do - get all homs then pick out an
    # injective one. But this takes too long.

    # keep degree small otherwise listing homs will take forever
    #degree := 10;
    #H := SymmetricGroup(degree);

    #all_homs := AllHomomorphismClasses(G, H);
    #rho := FuncToHom@(G, g -> ());
    #repeat
    #    rho := Random(all_homs);
    #until IsInjective(rho);

    # So we just get the easiest one GAP uses (might be too easy...)
    rho := IsomorphismPermGroup(G);

    return rec(rep := rho,
               G := id);
end;

# checks if the given centralizer basis is correct
TestCentralizerBasis@ := function(rep, cent_basis, args...)
    local correct_basis, C1, C2;

    correct_basis := rep.centralizer_basis;

    C1 := VectorSpace(Cyclotomics, correct_basis);
    C2 := VectorSpace(Cyclotomics, cent_basis);

    return C1 = C2;
end;

# checks if decomp is a decomposition into G-invariant subspaces
TestInvariantDecomposition@ := function(rep, decomp, args...)
    local rho, G;

    rho := rep.rep;
    G := Source(rho);

    return ForAll(decomp, V -> IsGInvariant@(rho, V));
end;

# checks some necessary conditions for decomp to be the full
# irreducible (collected) decomposition of rho
TestIrreducibleDecomposition@ := function(rep, decomp, args...)
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

    if ForAll(conds, x->x) then
        return true;
    else
        Error("test failed!");
        return false;
    fi;
end;

# same as irreducible decomp, these are only necessary conditions
# for correctness
TestCanonicalDecomposition@ := function(rep, decomp, args...)
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

    if ForAll(conds, x->x) then
        return true;
    else
        Error("test failed!");
        return false;
    fi;
end;

# test full BlockDiagonalRepresentation{Fast, Parallel} info
TestBlockDiagonalRepresentation@ := function(rep, info)
    local conds, A;

    A := TransposedMat(info.basis);

    conds := [];

    Add(conds, TestIrreducibleDecomposition@(rep, info.decomposition));

    Add(conds, TestCentralizerBasis@(rep, List(info.centralizer_basis, blocks -> A * BlockDiagonalMatrix(blocks) * A^-1)));

    if ForAll(conds, x->x) then
        return true;
    else
        Error("test failed!");
        return false;
    fi;
end;

# takes a function f : random representation -> boolean and tests it
# on n random representations
TestMany@ := function(f, n)
    local tested, rep;

    tested := 0;

    repeat
        rep := RandomRepresentation@();

        if not f(rep) then
            Print("FAILED: ", rep, "\n");
            return false;
        fi;

        tested := tested + 1;
    until tested = n;

    return true;
end;

# Same as TestMany, but with permutation representations
TestManyPerm@ := function(f, n)
    local tested, rep;

    tested := 0;

    repeat
        rep := RandomPermRepresentation@();

        if not f(rep) then
            Print("FAILED: ", rep, "\n");
            return false;
        fi;

        tested := tested + 1;
    until tested = n;

    return true;

end;
