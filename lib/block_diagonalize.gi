# Block diagonalizing representations

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
BasisChangeMatrix@ := function(rho)
    # Base change matrix from new basis to standard basis, acts on the left
    return TransposedMat(BlockDiagonalBasisOfRepresentation(rho));
end;

# Gives the nice basis for rho
InstallGlobalFunction( BlockDiagonalBasisOfRepresentation, function(rho)
    return ComputeUsingMethod@(rho).basis;
end );

# Takes a representation going to a matrix group and gives you an
# isomorphic representation where the images are block-diagonal with
# each block corresponding to an irreducible representation
InstallGlobalFunction( BlockDiagonalRepresentation, function(rho)
    return ComputeUsingMethod@(rho).diagonal_rep;
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

InstallGlobalFunction( BasisChangeMatrixSimilar, BasisChangeMatrixSimilar@ );
