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
