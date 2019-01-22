# These are the implementations of Serre's formulas from his book
# Linear Representations of Finite Groups.

# This is used for speed in some special cases
LoadPackage("cohcfg");

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
    local gens, ims, high, new_ims, new_range, new_rho, G;

    # We want rho to be a homomorphism to a matrix group since this
    # algorithm works on matrices. We convert a permutation group into
    # an isomorphic matrix group so that this is the case. If we don't
    # know how to convert to a matrix group, we just fail.

    if IsFiniteGroupLinearRepresentation(rho) then
        return rho;
    fi;

    if IsFiniteGroupPermutationRepresentation(rho) then
        G := Source(rho);
        gens := GeneratorsOfGroup(G);
        ims := List(gens, g -> Image(rho, g));
        high := LargestMovedPoint(ims);
        new_ims := List(ims, i -> PermutationMat(i, high));
        new_range := Group(new_ims);
        return GroupHomomorphismByImagesNC(G, new_range, gens, new_ims);
    fi;

    return fail;
end;

# Gives the canonical summand corresponding to irrep
IrrepCanonicalSummand@ := function(rho, irrep)
    local G, V, character, degree, projection, canonical_summand, H, T, cc, serre_class_contribution, two_orbit_reps, orbitals;

    G := Source(rho);

    degree := DegreeOfRepresentation(rho);

    # vector space rho(g) acts on
    V := Cyclotomics^degree;

    # In Serre's text, irrep is called W_i, this character is chi_i
    character := g -> Trace(Image(irrep, g));

    # Calculate the projection map from V to irrep using Theorem 8 (Serre)
    if not IsPermGroup(Range(rho)) or not IsInjective(rho) then
        # Given as a matrix, using Serre's formula directly, p_i is:
        projection := (degree/Order(G)) * Sum(G, t -> ComplexConjugate(character(t)) * Image(rho, t));
    else
        # TODO: check this trick works...
        H := Range(rho); # is perm group
        T := CohCfgFromPermGroup(H); # computes conjugacy classes and orbitals
        cc := ConjugacyClasses(H);
        two_orbit_reps := CCTransversal(T); # list of reps of 2-orbits (pairs)

        # the orbital matrices themselves
        orbitals := List(two_orbit_reps, OrbitalMatrix@);

        serre_class_contribution := function(class)
            local rep, v, A;
            rep := Representative(class);

            # this uses the fact the each orbital can be seen as the
            # adjacency matrix of a graph with appropriate isomorphism
            # group to sum the conjugacy class

            # coeffs of orbitals making up the class sum
            v := ClassSum(H, class);
            A := Sum([1..Length(v)], i -> v[i] * orbitals[i]);
            return ComplexConjugate(character(rep)) * A;
        end;

        projection := (degree/Order(G)) * Sum(cc, serre_class_contribution);
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
    irr_chars := IrrWithCorrectOrdering@(G, irreps);
    char_rho := DecomposeCharacter@(rho, irr_chars);

    # indices of irreps that appear in rho
    relevant_indices := Filtered([1..Length(irreps)],
                                 i -> char_rho[i] > 0);

    return irreps{relevant_indices};
end;

InstallMethod( CanonicalDecomposition, [ IsGroupHomomorphism ], function(rho)
    local G, F, n, V, irreps, chars, char_to_proj, canonical_projections, canonical_summands;

    # The group we are taking representations of
    G := Source(rho);

    # The list of irreps W_i of G over F appearing in rho
    # We need to convert here, since this function needs a linear rep
    irreps := RelevantIrreps@(ConvertRhoIfNeeded@(rho),
                              IrreducibleRepresentations(G));

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
    projections := List([1..n], a -> (n/Order(G)) * Sum(G, t -> Image(irrep,t^-1)[1][a]*Image(rho,t)));

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
InstallMethod( IrreducibleDecompositionCollected, "for linear representations", [ IsGroupHomomorphism ], function(arg_rho)
    local irreps, N, canonical_summands, full_decomposition, G, F, n, V, gens, ims, high, new_ims, new_range, rho;

    rho := ConvertRhoIfNeeded@(arg_rho);

    G := Source(rho);

    F := Cyclotomics;
    n := Length(Range(rho).1);
    V := F^n;

    irreps := RelevantIrreps@(rho, IrreducibleRepresentations(G));

    N := Size(irreps);

    # This gives a list of lists of vector spaces, each a
    # decomposition of a canonical summand into irreducibles.
    full_decomposition := List(irreps,
                               irrep -> DecomposeCanonicalSummand@(rho,
                                                                   irrep,
                                                                   IrrepCanonicalSummand@(rho,
                                                                                          irrep)));

    # Here we return the rho we actually used i.e. after we convert to
    # an isomorphic rep that goes to a matrix group (not a permutation
    # group)
    return rec(decomp := full_decomposition, used_rho := rho);
end );

# Gives the list of vector spaces in the direct sum decomposition of
# rho : G -> GL(V) into irreducibles.
InstallMethod( IrreducibleDecomposition, "for linear representations", [ IsGroupHomomorphism ], function(arg_rho)
    local rho;
    rho := ConvertRhoIfNeeded@(arg_rho);
    # We only want to return the vector spaces here
    return Flat(List(IrreducibleDecompositionCollected(rho).decomp,
                     rec_list -> List(rec_list, r -> r.space)));
end );
