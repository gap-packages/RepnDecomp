# Block diagonalizing representations

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
BaseChangeMatrix@ := function(rho, decomp)
    local new_bases, new_basis;

    # Extract the basis vectors, this is now a list of lists of bases
    # (each basis is a list of vectors)
    new_bases := List(decomp,
                      rec_list -> List(rec_list, r -> r.basis));

    # List of new basis row vectors
    new_basis := Concatenation(Concatenation(new_bases));

    # Base change matrix from new basis to standard basis
    return TransposedMat(new_basis);
end;

# Takes a representation going to a matrix group and gives you an
# isomorphic representation where the images are block-diagonal with
# each block corresponding to an irreducible representation
InstallGlobalFunction( BlockDiagonalizeRepresentation, function(rho, arg...)
    local decomp, A, G, gens, imgs, range;

    # Use precomputed decomposition, if available
    if Size(arg) > 0 then
        decomp := arg[1];
    else
        decomp := DecomposeIsomorphicCollected@(rho);
    fi;

    A := BaseChangeMatrix@(rho, decomp.decomp);
    G := Source(rho);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, g -> A^-1 * Image(decomp.used_rho, g) * A);

    range := Group(imgs);

    return GroupHomomorphismByImages(G, range, gens, imgs);
end );

# Calculates a matrix A such that X = A Y A^-1
# TODO: Make this work
BasisChangeMatrixSimilar@ := function(X, Y)
    local new_eigenspaces, old_eigenspaces;
    # TODO: Do the Jordan blocks come in the same order?
    new_eigenspaces := GeneralizedEigenspaces(Cyclotomics, X);
    old_eigenspaces := GeneralizedEigenspaces(Cyclotomics, Y);
    return 1;
end;

# Calculate a basis change matrix that diagonalizes rho (without using
# Serre's formulas)
BasisChangeMatrixAlternate@ := function(rho, args...)
    local G, char_rho_basis, irreps, isomorphic_collected, summands, new_rho_f, new_img, g;

    G := Source(rho);

    # Write the character of rho in the basis of irreducible characters
    char_rho_basis := DecomposeCharacter@(rho);

    if Size(args) > 0 then
        irreps := args[1];
    else
        irreps := IrreducibleRepresentations(G, Cyclotomics);
    fi;

    # TODO: Check if the ordering is safe to rely on!
    # Relying on the ordering of the basis, make a list of irreps in
    # the decomposition of rho.
    # The list of summands with isomorphic summands collected: we just
    # repeat irrep[i] the number of times given by the coefficient of
    # its character in char_rho
    isomorphic_collected := List([1..Size(char_rho_basis)],
                                 i -> Replicate(irreps[i], char_rho_basis[i]));

    summands := Flat(isomorphic_collected);

    new_rho_f := function(g)
        # Take the image with each direct summand and just glue them together
        return BlockDiagonalMatrix(List(summands, irrep -> Image(irrep, g)));
    end;

    # We don't know the basis that the new_rho(g) are written in, but
    # since the representations are isomorphic, there is a basis
    # change matrix A such that new_rho(g) = A * rho(g) * A^-1

    # To calculate A, we use a random element of G
    g := G.1;

    # TODO: Make this work!
    return BasisChangeMatrixSimilar@(new_rho_f(g), Image(rho, g));
end;
