# Correctness checks, so we can see if our decompositions, centralizer
# bases, etc are correct

RandomGroup@ := function(args...)
    local size, id;
    if Length(args) >= 1 then
        size := Random([args[1].lo..args[1].hi]);
    else
        size := Random([2..100]); # don't want tests to take forever
    fi;
    id := Random([1..NrSmallGroups(size)]);
    return [size, id];
end;

# Generate a random representation of a random group
RandomRepresentation@ := function(args...)
    local x, id, G, irreps, chosen_irreps, rho, n, A, centralizer_basis, isomorphism_type, diag_rho, opt, degree_so_far, new_irrep;

    if Length(args) >= 1 then
        opt := args[1];
    else
        opt := rec(lo := 2,
                   hi := 50,
                   num_irreps := 2,
                   min_multiplicity := 1,
                   max_multiplicity := 2,
                   max_total_degree := 4,
                   restrict_small_degree := true,
                   small_degree := 4);
    fi;

    id := RandomGroup@(opt);

    G := SmallGroup(id[1], id[2]);

    irreps := IrreducibleRepresentations(G);

    # in some cases we might want to restrict to only small degree
    if opt.restrict_small_degree then
        irreps := Filtered(irreps, irrep -> DegreeOfRepresentation(irrep) <= opt.small_degree);
    fi;

    degree_so_far := 0;
    chosen_irreps := [];

    for x in [1..opt.num_irreps] do
        # if the degree total is too high, bail out!
        new_irrep := Random(irreps);
        Add(chosen_irreps, new_irrep);
        degree_so_far := degree_so_far + DegreeOfRepresentation(new_irrep);

        if degree_so_far >= opt.max_total_degree then
            break;
        fi;
    od;

    # eliminate duplicates to avoid counting issues
    chosen_irreps := Set(chosen_irreps);

    # add some random multiplicity

    chosen_irreps := List(chosen_irreps, irrep -> Replicate@(irrep, Random([opt.min_multiplicity..opt.max_multiplicity])));

    diag_rho := DirectSumOfRepresentations(Flat(chosen_irreps));
    n := DegreeOfRepresentation(diag_rho);

    # scramble the basis to avoid it being too easy
    A := RandomInvertibleMat(n);
    rho := ComposeHomFunction(diag_rho, x -> A^-1 * x * A);

    # We know some things about rho without needing to compute them
    # the long way. We return them for use in testing.

    isomorphism_type := List(chosen_irreps, irrep_list -> rec(rep := irrep_list[1], m := Length(irrep_list)));
    centralizer_basis := List(chosen_irreps, irrep_list -> rec(dimension := DegreeOfRepresentation(irrep_list[1]),
                                                               nblocks := Length(irrep_list)));

    centralizer_basis := List(SizesToBlocks(centralizer_basis), BlockDiagonalMatrix);

    # put them in the scrambled basis
    centralizer_basis := List(centralizer_basis, M -> A^-1 * M * A);

    return rec(rep := rho,
               diag_rep := diag_rho,
               isomorphism_type := isomorphism_type,
               centralizer_basis := centralizer_basis,
               candidate_nice_basis := TransposedMat(A), # note this is not the unique right answer...
               irreps := irreps, # so we can cheat if we want to for e.g. benchmarks
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

# Tests if two canonical decompositions are "the same"
TestCanonicalDecompsEqual@ := function(decomp1, decomp2)
    local conds;

    conds := [];

    # Need to have same number of summands
    Add(conds, Length(decomp1) = Length(decomp2));

    # Decomps should be the same up to reordering and dimension
    Add(conds, SortedList(List(decomp1, Dimension)) = SortedList(List(decomp2, Dimension)));

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

    # want to change back to the old basis to verify centralisers are the same
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
        Reset(GlobalMersenneTwister, NanosecondsSinceEpoch());;
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
        Reset(GlobalMersenneTwister, NanosecondsSinceEpoch());;
        rep := RandomPermRepresentation@();

        if not f(rep) then
            Print("FAILED: ", rep, "\n");
            return false;
        fi;

        tested := tested + 1;
    until tested = n;

    return true;

end;

# benchmarks a function f : random representation -> anything for n
# random representations, printing output in the format:
#
# group size | group id | nr conjugacy classes | degree | time taken
#
# output is printed to the file with name out
BenchMany@ := function(f, out, n)
    local rep, size, id, num_classes, degree, t0, t1, time_taken, tested, opt;

    tested := 0;

    # clear file
    PrintTo(out);

    # this set of options to the random representation generator is
    # designed to avoid trivial cases
    opt := rec(lo := 1500,
               hi := 2000,
               num_irreps := 1,
               min_multiplicity := 1,
               max_multiplicity := 1,
               restrict_small_degree := true,
               small_degree := 10);

    # We deliberately avoid resetting the GlobalMersenneTwister since
    # we want the same reps to come up when you bench. This doesn't
    # really matter, but it's nicer.
    repeat
        rep := RandomRepresentation@(opt);
        size := rep.G[1];
        id := rep.G[2];
        num_classes := Length(ConjugacyClasses(Source(rep.rep)));
        degree := DegreeOfRepresentation(rep.rep);

        # do the bench
        t0 := Runtime();
        f(rep);
        t1 := Runtime();
        time_taken := t1 - t0;

        AppendTo(out, size, " ", id, " ", num_classes, " ", degree, " ", time_taken, "\n");

        tested := tested + 1;
        Print("done ", tested, "\n");
    until tested = n;
end;
