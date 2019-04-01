InstallGlobalFunction( IsUnitaryRepresentation, function(rho)
    local G;
    # note we only need to check for generators since a product of
    # unitary matrices is unitary
    G := Source(rho);
    return ForAll(GeneratorsOfGroup(G), g -> Image(rho, g^-1) = ConjugateTranspose@(Image(rho, g)));
end );

InstallGlobalFunction( LDLDecomposition, function(A)
    local n, L, D, i, j;

    n := Length(A);

    L := IdentityMat(n);

    # This stores the diagonal of D, not the full matrix
    D := List([1..n], x -> 0);

    for i in [1..n] do
        D[i] := A[i][i] - Sum([1..(i-1)], j -> L[i][j]*ComplexConjugate(L[i][j])*D[j]);
        for j in [i+1..n] do
            L[j][i] := (A[j][i] - Sum([1..(i-1)], k -> L[j][k]*ComplexConjugate(L[i][k])*D[k]))/D[i];
        od;
    od;

    return rec(L := L, D := D);
end );

InstallGlobalFunction( UnitaryRepresentation, function(rho)
    local G, S, decomp, L, D;

    G := Source(rho);

    # If rho is already unitary, then this will be cI for some
    # c>0. Otherwise, since we know rho can be made unitary, we know S
    # can be diagonalised and that transformation also unitarises rho.

    # TODO: do this sum quickly
    S := Sum(G, g -> Image(rho, g) * ConjugateTranspose@(Image(rho, g)));

    decomp := LDLDecomposition(S);
    L := decomp.L;
    D := DiagonalMat(decomp.D);

end );
