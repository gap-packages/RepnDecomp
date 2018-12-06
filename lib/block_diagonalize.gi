# Block diagonalizing representations

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
BasisChangeMatrix@ := function(rho)
    # Base change matrix from new basis to standard basis, acts on the left
    return TransposedMat(BlockDiagonalBasis(rho));
end;

# Gives the nice basis for rho
InstallMethod( BlockDiagonalBasis, "for linear reps", [ IsGroupHomomorphism ], function(arg_rho)
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
InstallMethod( BlockDiagonalRepresentation, "for linear reps", [ IsGroupHomomorphism ], function(arg_rho)
    local decomp, A, G, gens, imgs, range, rho;

    rho := ConvertRhoIfNeeded@(arg_rho);

    A := BasisChangeMatrix@(rho);
    G := Source(rho);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, g -> A^-1 * Image(rho, g) * A);

    range := Group(imgs);

    return GroupHomomorphismByImages(G, range, gens, imgs);
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

# Calculate a basis change matrix that diagonalizes rho (without using
# Serre's formulas)
BasisChangeMatrixAlternate@ := function(rho, args...)
    local G, char_rho_basis, irreps, isomorphic_collected, summands, new_rho_f, new_img, g, basis_change, basis, full_space_list, current_space_list, chars, new_rho, irrep_list, r;

    G := Source(rho);

    # Write the character of rho in the basis of irreducible characters
    char_rho_basis := DecomposeCharacter@(rho);

    irreps := IrreducibleRepresentationsDixon(G);

    chars := Irr(G);

    # TODO: Check if the ordering is safe to rely on!
    # Relying on the ordering of the basis, make a list of irreps in
    # the decomposition of rho.
    # The list of summands with isomorphic summands collected: we just
    # repeat irrep[i] the number of times given by the coefficient of
    # its character in char_rho
    isomorphic_collected := List([1..Size(char_rho_basis)],
                                 i -> Replicate@(rec(rep := irreps[i],
                                                     dim := chars[i][1]),
                                                 char_rho_basis[i]));

    summands := List(Flat(isomorphic_collected), r -> r.rep);

    new_rho_f := function(g)
        # Take the image with each direct summand and just glue them together
        return BlockDiagonalMatrix(List(summands, irrep -> Image(irrep, g)));
    end;

    new_rho := GroupHomomorphismByImages(G, Group(List(GeneratorsOfGroup(G), new_rho_f)));

    # We don't know the basis that the new_rho(g) are written in, but
    # since the representations are isomorphic, there is a basis
    # change matrix A such that new_rho(g) = A^-1 * rho(g) * A

    # To calculate A, we use a random element of G
    g := G.1;

    basis_change := BasisChangeMatrixSimilar@(new_rho_f(g), Image(rho, g));

    basis := TransposedMat(basis_change);

    # The basis is in the right order, it just needs to be collected
    # into bases for the irrep spaces
    full_space_list := [];
    for irrep_list in isomorphic_collected do
        current_space_list := [];
        for r in irrep_list do
            Add(current_space_list, VectorSpace(Cyclotomics, Take@(basis, r.dim)));
            basis := Drop@(basis, r.dim);
        od;
        Add(full_space_list, current_space_list);
    od;

    return rec(basis := basis,
               diagonal_rep := new_rho,
               decomposition := full_space_list);
end;
