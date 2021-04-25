InstallGlobalFunction( ConjugateTransformationNorm,
function(rho)
    local G, rho_conjugate, X, mu_I, mu;
    G := Source(rho);
    rho_conjugate := ComposeHomFunction(rho, M -> ComplexConjugate(M));

    if not AreRepsIsomorphic(rho, rho_conjugate) then
        Error("<rho> is not isomorphic to <rho> conjugate!");
    fi;

    X := LinearRepresentationIsomorphism(rho, rho_conjugate);

    # According to a paper or something, this will be \mu I
    mu_I := X * ComplexConjugate(X);

    mu := mu_I[1][1];

    # Now we just need to find theta such that N(theta) = mu.
    return mu;
end );
