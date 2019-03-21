# Block diagonalizing representations

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
BasisChangeMatrix@ := function(rho)
    # Base change matrix from new basis to standard basis, acts on the left
    return TransposedMat(BlockDiagonalBasisOfRepresentation(rho));
end;

# Gives the nice basis for rho
InstallGlobalFunction( BlockDiagonalBasisOfRepresentation, function(arg_rho)
    local new_bases, new_basis, rho;

    rho := ConvertRhoIfNeeded@(arg_rho);

    # Extract the basis vectors, this is now a list of lists of bases
    # (each basis is a list of vectors)
    new_bases := List(IrreducibleDecompositionCollected(rho).decomp,
                      rec_list -> List(rec_list, r -> r.basis));

    # List of new basis as row vectors
    return Concatenation(Concatenation(new_bases));
end );

# Takes a representation going to a matrix group and gives you an
# isomorphic representation where the images are block-diagonal with
# each block corresponding to an irreducible representation
InstallGlobalFunction( BlockDiagonalRepresentation, function(arg_rho)
    local decomp, A, G, gens, imgs, range, rho;

    rho := ConvertRhoIfNeeded@(arg_rho);
    A := BasisChangeMatrix@(rho);

    return ComposeHomFunction(rho, x -> A^-1 * x * A);
end );

# Calculates a matrix P such that X = P^-1 Y P
BasisChangeMatrixSimilar@ := function(X, Y)
    local A, B;

    # We find the rational canonical form conjugation matrices
    A := RationalCanonicalFormTransform(X);
    B := RationalCanonicalFormTransform(Y);

    # Now A^-1 X A = B^-1 Y B, so P = BA^-1
    return B * A^-1;
end;
