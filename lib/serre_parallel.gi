LoadPackage("IO");

# Same as IrreducibleDecompositionCollectedHybrid, but does the
# calculation for each irrep in parallel.
InstallGlobalFunction( IrreducibleDecompositionCollectedParallel, function(rho, num_jobs, args...)
    local irreps, G, full_decomposition;
    G := Source(rho);

    if Length(args) > 0 then
        irreps := RelevantIrreps@(rho, args[1]);
    else
        irreps := RelevantIrreps@(rho, IrreducibleRepresentations(G));
    fi;

    # This is the full info each process needs to calculate a
    # canonical summand and decompose it for an irrep
    full_decomposition := ParListByFork(irreps,
                                        irrep -> DecomposeCanonicalSummandFast@(rho,
                                                                                irrep,
                                                                                IrrepCanonicalSummand@(rho,
                                                                                                       irrep)),
                                        rec(NumberJobs := num_jobs));

    return rec(decomp := full_decomposition, used_rho := rho);

end );

InstallGlobalFunction( BlockDiagonalRepresentationParallel, function(rho, num_jobs, args...)
    local irreps, G, decomp, spaces_collected, block_sizes, new_bases, block_diag_basis, base_change_mat, diag_rho;
    G := Source(rho);
    if Length(args) > 0 then
        irreps := args[1];
    else
        irreps := IrreducibleRepresentations(G);
    fi;

    decomp := IrreducibleDecompositionCollectedParallel(rho, num_jobs, irreps).decomp;
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
               centralizer_basis := SizesToBlocks@(block_sizes));

end );
