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
    local G, S, n, prod, T, S2, i, j, decomp, L, D, s, unitary_rep;

    G := Source(rho);

    # If rho is already unitary, then this will be cI for some
    # c>0. Otherwise, since we know rho can be made unitary, we know S
    # can be diagonalised and that transformation also unitarises rho.

    # This is the slow version of the sum, we don't use it
    # S := Sum(G, g -> Image(rho, g) * ConjugateTranspose@(Image(rho, g)));

    n := DegreeOfRepresentation(rho);

    # fast version
    prod := g -> KroneckerProduct(Image(rho, g), ComplexConjugate(Image(rho, g)));
    T := GroupSumBSGS(G, prod);
    S := NullMat(n,n);

    for i in [1..n] do
        for j in [1..n] do
            S[i][j] := Sum([1..n], k -> ExtractBlock@(T, i, k, n)[j][k]);
        od;
    od;

    decomp := LDLDecomposition(S);

    L := decomp.L;
    D := decomp.D;

    s := List(D, x -> 1);

    # We want to scale D so it's |G|I
    #s := List(D, c -> c/Size(G));

    # then sqrt so that (Ls)(|G|I)(Ls)^*
    #s := List(s, Sqrt);

    return rec(L := L * DiagonalMat(s),
               D := DiagonalMat(s)^2 * decomp.D,
               unitary_rep := ComposeHomFunction(rho, m -> decomp.L^-1 * m * decomp.L));
end );
