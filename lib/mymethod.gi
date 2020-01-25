# Functions for calling the methods implemented my original algorithm.

InstallMethod( REPN_ComputeUsingMyMethod, "for linear reps", [ IsGroupHomomorphism ], function(rho)
    local G, char_rho_basis, irreps, isomorphic_collected, summands, new_img, g, basis_change, basis, full_space_list, current_space_list, chars, new_rho, irrep_list, r, ret_basis, all_sizes, centralizer_blocks, rho_cent_basis, new_rho_cent_basis;

    G := Source(rho);

    # If we are given a list of irreps, we use them and assume that it
    # is complete (or, as complete as we need)
    irreps := ValueOption("irreps");
    if irreps = fail then
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    # we might also be given a basis for the centraliser of rho, which
    # we can use to speed up some calculations
    rho_cent_basis := ValueOption("centralizer_basis");

    # We could just use Irr(G) here, but the ordering of Irr and
    # irreps don't necessarily match up. Also it may be that we want
    # to exclude some irreps for some reason (e.g. avoiding
    # cyclotomics).
    chars := ValueOption("irr_chars");
    if chars = fail then
        chars := IrrWithCorrectOrdering@(G : irreps := irreps);
    fi;

    # Write the character of rho in the basis of irreducible characters
    char_rho_basis := IrrVectorOfRepresentation@(rho, chars);

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

    # We can also compute the basis for the centralizer ring, since we
    # know the block sizes and relevant dimensions
    all_sizes := List([1..Size(chars)], i -> rec(dimension := chars[i][1],
                                                 nblocks := char_rho_basis[i]));

    # Don't use the blocks that don't appear
    centralizer_blocks := SizesToBlocks(Filtered(all_sizes, r -> r.nblocks > 0));

    new_rho_cent_basis := List(centralizer_blocks, BlockDiagonalMatrix@);

    # This is the block diagonal rep with minimal size blocks
    new_rho := DirectSumOfRepresentations(summands);

    # We don't know the basis that the new_rho(g) are written in, but
    # since the representations are isomorphic, there is a basis
    # change matrix A such that new_rho(g) = A^-1 * rho(g) * A for all g.
    #
    # This is an intertwining operator for rho and new_rho, or
    # representation isomorphism.

    # Note: this is where the heavy lifting of the function is

    basis_change := LinearRepresentationIsomorphism(new_rho, rho,
                                                    new_rho_cent_basis,
                                                    rho_cent_basis);

    basis := TransposedMat(basis_change);

    # We make a copy since we're going to return this one.
    ret_basis := ShallowCopy(basis);

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


    return rec(basis := ret_basis,
               diagonal_rep := new_rho,
               decomposition := Filtered(full_space_list,
                                         l -> Size(l) > 0), # don't use blocks that don't appear
               centralizer_basis := centralizer_blocks);
end );

# The same as REPN_ComputeUsingMyMethod but we first split into
# canonical summands which could be faster (as always, it's only
# actually faster if you have a basis for the centraliser)
InstallMethod( REPN_ComputeUsingMyMethodCanonical, "for linear reps", [ IsGroupHomomorphism ], function(rho)
    local G, cent_basis, decomp, spaces_collected, block_sizes, new_basis, block_diag_basis, base_change_mat, diag_rho, irreps, new_bases;
    G := Source(rho);

    irreps := ValueOption("irreps");
    if irreps = fail then
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    cent_basis := ValueOption("centralizer_basis");

    decomp := IrreducibleDecompositionCollectedHybrid@(rho : irreps := irreps).decomp;
    spaces_collected := List(decomp, rec_list -> List(rec_list, r -> VectorSpace(Cyclotomics, r.basis)));
    block_sizes := List(spaces_collected, l -> rec(dimension := Dimension(l[1]),
                                                   nblocks := Length(l)));

    new_bases := List(decomp, rec_list -> List(rec_list, r -> r.basis));

    # List of new basis as row vectors
    block_diag_basis := Concatenation(Concatenation(new_bases));

    base_change_mat := TransposedMat(block_diag_basis);

    diag_rho := ComposeHomFunction(rho, A -> base_change_mat^-1 * A * base_change_mat);

    return rec(basis := block_diag_basis,
               diagonal_rep := diag_rho,
               decomposition := spaces_collected,
               centralizer_basis := SizesToBlocks(block_sizes));
end );
