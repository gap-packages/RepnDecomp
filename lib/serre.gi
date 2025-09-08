# These are the implementations of Serre's formulas from his book
# Linear Representations of Finite Groups.

# Divides mat into nxn blocks, returns the matrix corresponding to the
# block in the ath row, bth column
ExtractBlock@ := function(mat, a, b, n)
    local result, row, i, j;
    result := [];
    for i in [1..n] do
        row := [];
        for j in [1..n] do
            Add(row, mat[i+(a-1)*n][j+(b-1)*n]);
        od;
        Add(result, row);
    od;
    return result;
end;

MatrixImage@ := function(p, V)
    local F;

    # F is the base field of V
    F := LeftActingDomain(V);

    # We return the span of the images of the basis under p, which
    # gives p(V)
    return VectorSpace(F,
                       List(Basis(V), v -> p * v),
                       Zero(V));
end;

# Converts rho to a matrix representation if necessary
ConvertRhoIfNeeded@ := function(rho)
    # We want rho to be a homomorphism to a matrix group since this
    # algorithm works on matrices. We convert a permutation group into
    # an isomorphic matrix group so that this is the case. If we don't
    # know how to convert to a matrix group, we just fail.

    if IsFiniteGroupLinearRepresentation(rho) then
        return rho;
    fi;

    if IsFiniteGroupPermutationRepresentation(rho) then
        return PermToLinearRep(rho);
    fi;

    return fail;
end;

# assumes rho is faithful permutation representation, calculates centralizer
RepresentationCentralizerPermRep@ := function(rho)
    local G, H, T, two_orbit_reps;
    G := Source(rho);
    H := Group(List(GeneratorsOfGroup(G), g -> Image(rho, g))); # is perm group
    T := CohCfgFromPermGroup@(H); # computes conjugacy classes and orbitals
    two_orbit_reps := CCTransversal@(T); # list of reps of 2-orbits (pairs)

    # the orbital matrices themselves, this gives a basis for the
    # centraliser
    return List(two_orbit_reps, rep -> OrbitalMatrix@(H, rep));
end;

# calculates the matrix with blocks p_ab
PMatrix@ := function(rho, irrep)
    local G, irrep_dual, tau, p;
    G := Source(rho);
    irrep_dual := FuncToHom@(G, g -> TransposedMat(Image(irrep, g^-1)));
    tau := KroneckerProductOfRepresentations(irrep_dual, rho);
    p := GroupSumBSGS(G, g -> Image(tau, g));
    return p;
end;

# Gives the canonical summand (basis) corresponding to irrep
IrrepCanonicalSummand@ := function(rho, irrep)
    local G, degree, V, cent_basis, opt, p, projection, character,
          cc, summand, canonical_summand;

    G := Source(rho);

    degree := DegreeOfRepresentation(rho);

    # vector space rho(g) acts on
    V := Cyclotomics^degree;

    # if we are given an orthonormal basis for the centralizer of rho,
    # then we can use it to speed up class sum computation
    cent_basis := ValueOption("centralizer_basis");

    # sometimes we may want to disable all optimisations
    opt := ValueOption("no_optimisations") <> true;

    # In Serre's text, irrep is called W_i, the character is chi_i

    # Now we calculate the projection map from V to irrep using
    # Theorem (Serre)

    # first we check if we are allowed to use the kronecker products,
    # this lets us use the fact that p_i = \sum_\alpha
    # p_{\alpha\alpha}, we can calculate these now and pass these to
    # the canonical summand breakdown algorithm later
    if opt and ValueOption("use_kronecker") = true then
        p := PMatrix@(rho, irrep);
        projection := Sum([1..DegreeOfRepresentation(irrep)],
                          alpha -> ExtractBlock@(p, alpha, alpha, degree));
    elif opt and cent_basis <> fail then
        # if we are given a basis for the centralizer, we assume the
        # user knows what they're doing - the rep is unitary and the
        # basis is orthonormal, this is not checked (too expensive).
        character := g -> Trace(Image(irrep, g));
        cc := ConjugacyClasses(G);
        projection := (degree/Order(G)) * Sum(cc,
          cl -> ComplexConjugate(character(Representative(cl))) * ClassSumCentralizerNC(rho, cl, cent_basis));
    else
        # Lastly, given no special info at all we just have to sum over G
        character := g -> Trace(Image(irrep, g));

        # This maps t to the summand from Serre's formula
        summand := t -> ComplexConjugate(character(t)) * Image(rho, t);

        projection := (degree/Order(G)) * Sum(G, summand);
    fi;

    # Calculate V_i, the canonical summand
    canonical_summand := MatrixImage@(projection, V);
    canonical_summand := Basis(canonical_summand);

    # if we calculated using kronecker products, we calculated the
    # p_ab, let's return it so the later stages can use it
    if ValueOption("use_kronecker") = true then
        return rec(space := canonical_summand, p := p);
    else
        return canonical_summand;
    fi;
end;

# filters out the irreps that don't appear in the direct sum
# decomposition of rho
RelevantIrreps@ := function(rho, irreps)
    local G, irr_chars, char_rho, relevant_indices;

    G := Source(rho);
    irr_chars := IrrWithCorrectOrdering@(G : irreps := irreps);
    char_rho := IrrVectorOfRepresentation@(rho, irr_chars);

    # indices of irreps that appear in rho
    relevant_indices := Filtered([1..Length(irreps)],
                                 i -> char_rho[i] > 0);

    return irreps{relevant_indices};
end;

InstallGlobalFunction( CanonicalDecomposition, function(rho)
    local G, irreps;

    # The group we are taking representations of
    G := Source(rho);

    # The list of irreps W_i of G over F appearing in rho
    irreps := ValueOption("irreps");

    if irreps = fail then
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    # We might need to convert here, since this function needs a
    # linear rep
    irreps := RelevantIrreps@(ConvertRhoIfNeeded@(rho),
                              irreps);

    # if there's only 1 irrep, the canonical summand is just the whole
    # space!
    if Length(irreps) = 1 then
        return [Cyclotomics^DegreeOfRepresentation(rho)];
    fi;

    # otherwise do the calculation per irrep

    # if given a basis for centralizer, it lives in the option stack
    # and gets used
    return List(irreps,
                irrep -> VectorSpace(Cyclotomics,
                                     IrrepCanonicalSummand@(rho, irrep),
                                     List([1..DegreeOfRepresentation(rho)], x -> 0)));
end );

# Decomposes the representation V_i into a direct sum of some number
# (maybe zero) of spaces, all isomorphic to W_i. W_i is the space
# corresponding to the irrep : G -> GL(W_i). rho is the "full"
# representation that we're decomposing.
DecomposeCanonicalSummand@ := function(rho, irrep, V_i)
    local projections, p_11, V_i1, basis, n, step_c, G, H, F, V, m, p, tau, irrep_dual, myspace;

    G := Source(irrep);

    # This is a subgroup of Aut(some space), we don't really know or
    # care what the space actually is
    H := Range(irrep);

    # This gives the dimension of the space, degree of irrep
    n := Length(H.1);
    F := Cyclotomics;

    # if the V_i is a record, then this means we have been given the
    # matrix with the blocks p_\alpha\beta
    if IsRecord(V_i) then
        p := V_i.p;
        myspace := V_i.space;
    else
        p := fail;
        myspace := V_i;
    fi;

    # it is possible that V_i (or V_i.space) was given as a list of
    # basis vectors, so we might need to convert it to a space
    if not IsVectorSpace(myspace) then
        myspace := VectorSpace(F, myspace, List([1..n], x -> 0));
    fi;

    # First compute the projections p_ab. We only actually use
    # projections with a=1..n and b=1, so we can just compute
    # those. projections[a] is p_{a1} from Serre. Here we can use a
    # neat trick to sum over a BSGS, but this requires some possibly
    # big matrices that might not fit in memory so we don't do it by
    # default.
    if ValueOption("use_kronecker") = true then
        # p = \sum_g (\rho_i^* \otimes \rho) so p_ab is the matrix we
        # get if we divide that sum into deg rho square blocks and pick the _ab
        # one
        if p = fail then
            p := PMatrix@(rho, irrep);
        fi;
        projections := List([1..n], a -> ExtractBlock@(p, a, 1, DegreeOfRepresentation(rho)));
    else
        projections := List([1..n], function(a)
                               local summand;
                               summand := t -> Image(irrep,t^-1)[1][a]*Image(rho,t);
                               return (n/Order(G)) * Sum(G, summand);
                           end );
    fi;

    p_11 := projections[1];
    V_i1 := MatrixImage@(p_11, myspace);
    basis := Basis(V_i1);

    # Now we define the map taking x to W(x), a subrepresentation of
    # V_i isomorphic to W_i. (This is step (c) of Proposition 8)
    step_c := function(x1)
        # This is the list of basis vectors for W(x1)
        return List([1..n], alpha -> projections[alpha] * x1);
    end;

    # If x1^1 .. x1^m is a basis for V_i1 (this is in the `basis`
    # variable), then V_i decomposes into the direct sum W(x1^1)
    # ... W(x1^m), each isomorphic to W_i.
    #
    # We return a list of lists of (vector space, basis) pairs where
    # the basis (TODO: confirm this?) has the special property
    return List(basis, function(x)
                   local b;
                   b := step_c(x);
                   # we don't need the actual space since we can
                   # recover it from the basis
                   return rec(#space := VectorSpace(F, b, Zero(V)),
                              basis := b);
               end);
end;

# Decompose rho into irreducible representations with the reps that
# are isomorphic collected together. This returns a list of lists of
# vector spaces (L) with each element of L being a list of vector
# spaces arising from the same irreducible.
InstallGlobalFunction( IrreducibleDecompositionCollected, function(rho)
    # user is probably only interested in the nonempty summands, but we
    # store all of them for consistency
    return Filtered(ComputeUsingMethod@(rho).decomposition, l -> Length(l) > 0);
end );

InstallGlobalFunction( IrreducibleDecomposition, function(rho)
    return Flat(IrreducibleDecompositionCollected(rho));
end );

# Uses REPN_ComputeUsingMyMethod to decompose a canonical
# summand. Ideally canonical summands are small compared to the whole
# rep, so could be faster than Serre's formula
DecomposeCanonicalSummandAlternate@ := function(rho, irrep, V_i)
    local basis, restricted_rho, block_diag_info, space_list,
          full_basis, big_space, bas_index, decomp, space,
          current_basis, _;

    if IsRecord(V_i) then
        basis := V_i.space;
    else
        basis := V_i;
    fi;

    # TODO: make this less confusing
    # V_i might be a vector space, we want a basis
    if IsVectorSpace(basis) then
        basis := Basis(basis);
    fi;

    restricted_rho := RestrictRep@(rho, basis);

    # we know in advance that restricted_rho only consists of direct
    # sums of irrep
    block_diag_info := REPN_ComputeUsingMyMethod(restricted_rho : irreps := [irrep]);

    # This is a list of spaces, with vectors written as coefficients
    # of the basis for V_i. We know this list only has 1 element,
    # there is only one irrep.
    space_list := block_diag_info.decomposition[1];

    # this is the nice basis, we can read off the basis for each
    # space_list[i] from this list, in order
    full_basis := block_diag_info.basis;

    # Convert back to big space vectors
    big_space := v -> Sum([1..Length(v)], i -> v[i] * basis[i]);
    full_basis := List(full_basis, big_space);

    bas_index := 1;

    decomp := [];

    for space in space_list do
        current_basis := [];
        for _ in [1..Dimension(space)] do
            Add(current_basis, full_basis[bas_index]);
            bas_index := bas_index + 1;
        od;
        Add(decomp, rec(basis := current_basis));
    od;

    # same format as DecomposeCanonicalSummand
    return decomp;
end;

# Uses Serre's formula to get canonical decomposition, then "fast"
# method to get irreducibles from that. Uses orthonormal basis for
# C_rho if given.
IrreducibleDecompositionCollectedHybrid@ := function(rho)
    local irreps, G, irred_decomp, do_decompose, parallel, ParListByFork;
    G := Source(rho);
    irreps := ValueOption("irreps");
    if irreps = fail then
        irreps := IrreducibleRepresentations(G);
    fi;
    irreps := RelevantIrreps@(rho, irreps);

    do_decompose := irrep -> DecomposeCanonicalSummandAlternate@(rho,
                                                                 irrep,
                                                                 IrrepCanonicalSummand@(rho,
                                                                                        irrep));
    parallel := ValueOption("parallel");
    if IsBoundGlobal("ParListByFork") then
        ParListByFork := ValueGlobal("ParListByFork");
    elif parallel <> fail then
        Error("The GAP package IO must be loaded to use the parallel option!");
    fi;
 
    if IsInt(parallel) then
        irred_decomp := ParListByFork(irreps, do_decompose, rec(NumberJobs := parallel));
    elif parallel <> fail then
        irred_decomp := ParListByFork(irreps, do_decompose, rec(NumberJobs := 4));
    else
        irred_decomp := List(irreps, do_decompose);
    fi;

    return rec(decomp := irred_decomp, used_rho := rho);
end;

InstallMethod( REPN_ComputeUsingSerre, "for linear reps", [ IsGroupHomomorphism ], function(rho)
    local irreps, irr_chars, centralizer_basis, irred_decomp, new_bases, basis, basis_change, diag_rho, char_rho_basis, all_sizes, sizes, centralizer_blocks, G, parallel, do_decompose, ParListByFork;

    G := Source(rho);

    irreps := ValueOption("irreps");
    if irreps = fail then
        irreps := IrreducibleRepresentations(G);
    fi;

    irr_chars := ValueOption("irr_chars");
    if irr_chars = fail then
        irr_chars := IrrWithCorrectOrdering@(G : irreps := irreps);
    fi;


    centralizer_basis := ValueOption("centralizer_basis");

    if centralizer_basis <> fail and not IsOrthonormalSet(centralizer_basis, InnerProduct@) then
        Error("<centralizer_basis> is not orthonormal!");
    fi;

    do_decompose := function(irrep)
        local canonical;
        canonical := IrrepCanonicalSummand@(rho, irrep : centralizer_basis := centralizer_basis);
        return DecomposeCanonicalSummand@(rho,
                                          irrep,
                                          canonical);
    end;

    parallel := ValueOption("parallel");
    if IsBoundGlobal("ParListByFork") then
        ParListByFork := ValueGlobal("ParListByFork");
    elif parallel <> fail then
        Error("The GAP package IO must be loaded to use the parallel option!");
    fi;

    if IsInt(parallel) then
        irred_decomp := ParListByFork(irreps, do_decompose, rec(NumberJobs := parallel));
    elif parallel <> fail then
        # we default the number of jobs to 4 since everyone
        # has 4 threads, at least
        irred_decomp := ParListByFork(irreps, do_decompose, rec(NumberJobs := 4));
    else
        irred_decomp := List(irreps, do_decompose);
    fi;

    new_bases := List(irred_decomp,
                      rec_list -> List(rec_list, r -> r.basis));
    basis := Concatenation(Concatenation(new_bases));
    basis_change := TransposedMat(basis);
    diag_rho := ComposeHomFunction(rho, x -> basis_change^-1 * x * basis_change);

    char_rho_basis := IrrVectorOfRepresentation@(rho, irr_chars);

    # Calculate sizes based on the fact irr_char[1] is the degree
    all_sizes := List([1..Size(irr_chars)],
                      i -> rec(dimension := irr_chars[i][1],
                               nblocks := char_rho_basis[i]));

    # Now we remove all of the ones with nblocks = 0 (doesn't affect
    # end result)
    sizes := Filtered(all_sizes, r -> r.nblocks > 0);

    centralizer_blocks := SizesToBlocks(sizes);

    return rec(basis := basis,
               diagonal_rep := diag_rho,
               decomposition := irred_decomp,
               centralizer_basis := centralizer_blocks);
end );
