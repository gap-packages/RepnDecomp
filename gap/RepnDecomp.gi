#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# Implementations
#

# Calculates p(V) for p a linear map (given as a matrix in the
# standard basis) and a vector space V
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

InstallGlobalFunction( DecomposeRepresentationCanonical, function(rho)
    local G, F, n, V, irreps, chars, char_to_proj, canonical_projections, canonical_summands;

    # The group we are taking representations of
    G := Source(rho);

    # The field we are working over: it's always the Cyclotomics
    F := Cyclotomics;

    # The dimension of the V in rho : G -> GL(V). Since we have the
    # images of rho as matrices, this is just the width or height of
    # any image of any generator of G.
    n := Length(Range(rho).1);

    # The vector space that the linear maps act on
    V := F^n;

    # The full list of irreps W_i of G over F
    irreps := IrreducibleRepresentations(G, F);

    # The characters chi_i of each irrep W_i
    chars := List(irreps,
                  irrep -> GroupHomomorphismByFunction(G, F,
                                                       g -> Trace(Image(irrep, g))));

    # Given a character chi_i, calculate the projection onto V_i using Theorem 8
    # This is given as a matrix
    char_to_proj := function(char)
        local degree;

        # The degree n_i of char
        degree := Image(char, One(G));
        return (degree/Order(G)) * Sum(G,
                                       t -> ComplexConjugate(Image(char, t)) * Image(rho, t));
    end;

    # The list of the p_i in matrix form
    canonical_projections := List(chars, char_to_proj);

    # The list of the V_i
    canonical_summands := List(canonical_projections, p -> MatrixImage@(p, V));

    return canonical_summands;
end );

# Decomposes the representation V_i into a direct sum of some number
# (maybe zero) of spaces, all isomorphic to W_i. W_i is the space
# corresponding to the irrep : G -> GL(W_i). rho is the "full"
# representation that we're decomposing.
DecomposeCanonicalSummand@ := function(rho, irrep, V_i)
    local projection, p_11, V_i1, basis, n, step_c, G, H, F, V, m;

    G := Source(irrep);

    # This is the general linear group of some space, we don't really
    # know or care what the space actually is
    H := Range(irrep);

    # This gives the dimension of the space of which W is the general
    # linear group (the size of the matrices representing the maps)
    n := Length(H.1);

    m := Length(Range(rho).1);
    F := Cyclotomics;
    V := F^m;

    # First compute the projections p_ab
    projection := function(a, b)
        return (n/Order(G))*Sum(Elements(G),
                                t -> Image(irrep,t^-1)[b][a]*Image(rho,t));
    end;

    p_11 := projection(1, 1);
    V_i1 := MatrixImage@(p_11, V_i);
    basis := Basis(V_i1);

    # Now we define the map taking x to W(x), a subrepresentation of
    # V_i isomorphic to W_i. (This is step (c) of Proposition 8)
    step_c := function(x1)
        # This is the list of basis vectors for W(x1)
        return List([1..n],
                    alpha -> projection(alpha, 1) * x1);
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
DecomposeIsomorphicCollected@ := function(orig_rho)
    local irreps, N, canonical_summands, full_decomposition, G, F, n, V, gens, ims, high, new_ims, new_range, rho;

    rho := orig_rho;
    G := Source(rho);

    # We want rho to be a homomorphism to a matrix group since this
    # algorithm works on matrices. We convert a permutation group into
    # an isomorphic matrix group so that this is the case. If we don't
    # know how to convert to a matrix group, we just fail.
    if not IsMatrixGroup(Range(rho)) then
        if IsPermGroup(Range(rho)) then
            gens := GeneratorsOfGroup(G);
            ims := List(gens, g -> Image(rho, g));
            high := LargestMovedPoint(ims);
            new_ims := List(ims, i -> PermutationMat(i, high));
            new_range := Group(new_ims);
            rho := GroupHomomorphismByImages(G, new_range, gens, new_ims);
        else
            Error("rho is not a matrix or permutation group!");
        fi;
    fi;

    F := Cyclotomics;
    n := Length(Range(rho).1);
    V := F^n;

    irreps := IrreducibleRepresentations(G, F);

    N := Size(irreps);

    # This gives a list of vector spaces, each a canonical summand
    canonical_summands := DecomposeRepresentationCanonical(rho);

    # This gives a list of lists of vector spaces, each a
    # decomposition of a canonical summand into irreducibles.
    full_decomposition := List([1..N],
                               i -> DecomposeCanonicalSummand@(rho, irreps[i], canonical_summands[i]));

    return full_decomposition;
end;

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
#
# TODO: Make it so that blocks from isomorphic representations are the
# same. Or maybe they already are, further testing required.
BlockDiagonalizeRepresentation@ := function(rho)
    local decomp, new_bases, new_basis;

    # First decompose rho, keeping the information about which irreps
    # are isomorphic and which is the special basis to use
    decomp := DecomposeIsomorphicCollected@(rho);

    # Extract the basis vectors, this is now a list of lists of bases
    # (each basis is a list of vectors)
    new_bases := List(decomp,
                      rec_list -> List(rec_list, r -> r.basis));

    # List of new basis row vectors
    new_basis := Concatenation(Concatenation(new_bases));

    # Base change matrix from new basis to standard basis
    return TransposedMat(new_basis);
end;

# This gives the list of vector spaces in the direct sum
# decomposition of rho : G -> GL(V) into irreducibles.
InstallGlobalFunction( DecomposeRepresentationIrreducible, function(orig_rho)
    # We only want to return the vector spaces here
    return Flat(List(DecomposeIsomorphicCollected@(orig_rho),
                     rec_list -> List(rec_list, r -> r.space)));
end );
