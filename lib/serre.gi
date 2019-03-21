LoadPackage("IO");

# These are the implementations of Serre's formulas from his book
# Linear Representations of Finite Groups.

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

# Gives the canonical summand corresponding to irrep
IrrepCanonicalSummand@ := function(rho, irrep)
    local G, V, character, degree, projection, canonical_summand, H, T, cc, serre_class_contribution, two_orbit_reps, orbitals, cent_basis, summand;

    G := Source(rho);
    #ker := KernelOfMultiplicativeGeneralMapping(arg_rho);
    #quotient_hom := NaturalHomomorphismByNormalSubgroupNC(G, ker);
    #rho := GroupHomomorphismByFunction(FactorGroup(G, ker),
    #                                   Range(arg_rho), gclass -> Image(arg_rho, PreImagesRepresentative(quotient_hom, gclass)));

    degree := DegreeOfRepresentation(rho);

    # vector space rho(g) acts on
    V := Cyclotomics^degree;

    # if we are given an orthonormal basis for the centralizer of rho,
    # then we can use it to speed up class sum computation
    cent_basis := ValueOption("centralizer_basis");

    # In Serre's text, irrep is called W_i, the character is chi_i

    # Calculate the projection map from V to irrep using Theorem 8 (Serre)
    if false and cent_basis <> fail then
        # First, if we are given a basis for the centralizer
        character := g -> Trace(Image(irrep, g));
        cc := ConjugacyClasses(G);
        projection := (degree/Order(G)) * Sum(cc,
          cl -> ComplexConjugate(character(Representative(cl))) * ClassSumCentralizer(rho, cl, cent_basis));
    # I want to check if rho is injective - otherwise this doesn't
    # work, but checking that checks each element of G has a different
    # image - slow.
    #elif IsPermGroup(Range(rho)) and IsInjective(rho) then
    elif IsPermGroup(Range(rho)) then
        # Then if we are given a perm group, but not a basis for the
        # centralizer, can calculate a basis for C from from scratch
        character := h -> Trace(Image(irrep, PreImagesRepresentative(rho, h)));
        H := Range(rho); # is perm group
        T := CohCfgFromPermGroup@(H); # computes conjugacy classes and orbitals
        cc := ConjugacyClasses(H);

        # the orbital matrices give a basis for the centralizer
        orbitals := RepresentationCentralizerPermRep@(rho);

        serre_class_contribution := function(class)
            local rep, v, A;
            rep := Representative(class);

            # coeffs of orbitals making up the class sum
            v := ClassSum@(T, class);

            # the class sum itself
            A := Sum([1..Length(v)], i -> v[i] * orbitals[i]);
            return ComplexConjugate(character(rep)) * A;
        end;

        projection := (degree/Order(G)) * Sum(cc, serre_class_contribution);
    else
        # Lastly, given no special info at all we just have to sum over G

        # Doing this with the StabChain method might actually be quite
        # fast, maybe faster than the other methods, need to investigate

        character := g -> Trace(Image(irrep, g));

        # This maps t to the summand from Serre's formula
        summand := t -> ComplexConjugate(character(t)) * Image(rho, t);

        projection := (degree/Order(G)) * Sum(G, summand);
    fi;

    # Calculate V_i, the canonical summand
    canonical_summand := MatrixImage@(projection, V);
    return canonical_summand;
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
    local G, F, n, V, irreps, chars, char_to_proj, canonical_projections, canonical_summands;

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
    return List(irreps, irrep -> IrrepCanonicalSummand@(rho, irrep));
end );

# Decomposes the representation V_i into a direct sum of some number
# (maybe zero) of spaces, all isomorphic to W_i. W_i is the space
# corresponding to the irrep : G -> GL(W_i). rho is the "full"
# representation that we're decomposing.
DecomposeCanonicalSummand@ := function(rho, irrep, V_i)
    local projections, p_11, V_i1, basis, n, step_c, G, H, F, V, m;

    G := Source(irrep);

    # This is the general linear group of some space, we don't really
    # know or care what the space actually is
    H := Range(irrep);

    # This gives the dimension of the space of which W is the general
    # linear group (the size of the matrices representing the maps),
    # is also the degree.
    n := Length(H.1);

    F := Cyclotomics;
    V := F^n;

    # First compute the projections p_ab. We only actually use projections with
    # a=1..n and b=1, so we can just compute those. projections[a] is p_{a1}
    # from Serre.
    projections := List([1..n], function(a)
                           local summand;
                           summand := t -> Image(irrep,t^-1)[1][a]*Image(rho,t);
                           return (n/Order(G)) * Sum(G, summand);
                       end );

    p_11 := projections[1];
    V_i1 := MatrixImage@(p_11, V_i);
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
                   return rec(space := VectorSpace(F, b, Zero(V)), basis := b);
               end);
end;

# Decompose rho into irreducible representations with the reps that
# are isomorphic collected together. This returns a list of lists of
# vector spaces (L) with each element of L being a list of vector
# spaces arising from the same irreducible.
InstallGlobalFunction( IrreducibleDecompositionCollected, function(rho)
    return ComputeUsingMethod@(rho).decomposition;
end );

# Uses REPN_ComputeUsingMyMethod to decompose a canonical
# summand. Ideally canonical summands are small compared to the whole
# rep, so could be faster than Serre's formula
DecomposeCanonicalSummandAlternate@ := function(rho, irrep, V_i)
    local basis, restricted_rho, block_diag_info, space_list, big_space;

    basis := Basis(V_i);
    restricted_rho := RestrictRep@(rho, basis);

    # we know in advance that restricted_rho only consists of direct
    # sums of irrep
    block_diag_info := REPN_ComputeUsingMyMethod(restricted_rho : irreps := [irrep]);

    # This is a list of spaces, with vectors written as coefficients
    # of the basis for V_i. We know this list only has 1 element,
    # there is only one irrep.
    space_list := block_diag_info.decomposition[1];

    # Convert back to big space vectors
    big_space := function(space)
        local orig, new;
        orig := Basis(space);
        new := List(orig, v -> Sum([1..Length(v)], i -> v[i] * basis[i]));
        return VectorSpace(Cyclotomics, new, Zero(V_i));
    end;

    space_list := List(space_list, big_space);

    # same format as DecomposeCanonicalSummand
    return List(space_list,
                space -> rec(basis := List(Basis(space))));
end;

# Uses Serre's formula to get canonical decomposition, then "fast"
# method to get irreducibles from that. Uses orthonormal basis for
# C_rho if given.
IrreducibleDecompositionCollectedHybrid@ := function(rho)
    local irreps, G, irred_decomp, do_decompose, parallel;
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
    if parallel = true then
        irred_decomp := ParListByFork(irreps, do_decompose);
    else
        irred_decomp := List(irreps, do_decompose);
    fi;

    return rec(decomp := irred_decomp, used_rho := rho);
end;

# TODO: move all calculations here and make the function access this record
InstallMethod( REPN_ComputeUsingSerre, "for linear reps", [ IsGroupHomomorphism ], function(rho)
    local irreps, irr_chars, centralizer_basis, irred_decomp, new_bases, basis, basis_change, diag_rho, char_rho_basis, all_sizes, sizes, centralizer_blocks, G, parallel, do_decompose;

    G := Source(rho);

    irreps := ValueOption("irreps");
    if irreps = fail then
        irreps := IrreducibleRepresentations(G);
    fi;

    irr_chars := ValueOption("irr_chars");
    if irr_chars = fail then
        irr_chars := IrrWithCorrectOrdering@(G);
    fi;


    centralizer_basis := ValueOption("centralizer_basis");

    do_decompose := function(irrep)
        local canonical;
        canonical := IrrepCanonicalSummand@(rho, irrep : centralizer_basis := centralizer_basis);
        return DecomposeCanonicalSummand@(rho,
                                          irrep,
                                          canonical);
    end;

    parallel := ValueOption("parallel");
    if parallel = true then
        irred_decomp := ParListByFork(irreps, do_decompose);
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
