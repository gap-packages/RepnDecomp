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

InstallGlobalFunction( BlockDiagonalRepresentationFast, function(rho, args...)
    local G, char_rho_basis, irreps, isomorphic_collected, summands, new_rho_f, new_img, g, basis_change, basis, full_space_list, current_space_list, chars, new_rho, irrep_list, r, F, ret_basis, all_sizes, centralizer_blocks;

    G := Source(rho);

    # If we are given a list of irreps, we use them and assume that it
    # is complete
    irreps := [];
    if Length(args) > 0 then
        irreps := args[1];
    else
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    F := Cyclotomics;
    # We can also be given a field to work over
    if Length(args) > 1 then
        F := args[2];
    fi;

    # We could just use Irr(G) here, but the ordering of Irr and
    # irreps don't necessarily match up. Also it may be that we want
    # to exclude some irreps for some reason (e.g. avoiding
    # cyclotomics).
    chars := IrrWithCorrectOrdering@(G, irreps);

    # Write the character of rho in the basis of irreducible characters
    char_rho_basis := DecomposeCharacter@(rho, chars);

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

    # We make a copy since we're going to return this one.
    ret_basis := ShallowCopy(basis);

    # The basis is in the right order, it just needs to be collected
    # into bases for the irrep spaces
    full_space_list := [];
    for irrep_list in isomorphic_collected do
        current_space_list := [];
        for r in irrep_list do
            Add(current_space_list, VectorSpace(F, Take@(basis, r.dim)));
            basis := Drop@(basis, r.dim);
        od;
        Add(full_space_list, current_space_list);
    od;

    # We can also compute the basis for the centralizer ring, since we
    # know the block sizes and relevant dimensions
    all_sizes := List([1..Size(chars)], i -> rec(dimension := chars[i][1],
                                                 nblocks := char_rho_basis[i]));

    # Don't use the blocks that don't appear
    centralizer_blocks := SizesToBlocks@(Filtered(all_sizes, r -> r.nblocks > 0));

    return rec(basis := ret_basis,
               diagonal_rep := new_rho,
               decomposition := full_space_list,
               centralizer_basis := centralizer_blocks);
end );
