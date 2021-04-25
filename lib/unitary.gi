InstallGlobalFunction( IsUnitaryRepresentation, function(rho)
    local G;
    # note we only need to check for generators since a product of
    # unitary matrices is unitary
    G := Source(rho);
    return ForAll(GeneratorsOfGroup(G), g -> Image(rho, g^-1) = ConjugateTranspose@(Image(rho, g)));
end );

InstallGlobalFunction( LDLDecomposition, function(A)
    local n, L, D, i, j;

    if A <> ConjugateTranspose@(A) then
        Error("<A> is not conjugate symmetric!");
        return fail;
    fi;

    n := Length(A);

    L := IdentityMat(n);

    # This stores the diagonal of D, not the full matrix
    D := List([1..n], x -> 1);

    for i in [1..n] do
        D[i] := A[i][i] - Sum([1..(i-1)], j -> L[i][j]*ComplexConjugate(L[i][j])*D[j]);
        for j in [i+1..n] do
            L[j][i] := (A[j][i] - Sum([1..(i-1)], k -> L[j][k]*ComplexConjugate(L[i][k])*D[k]))/D[i];
        od;
    od;

    return rec(L := L, D := D);
end );

InstallGlobalFunction( UnitaryRepresentation, function(rho)
    local G, n, prod, T, S, i, j, decomp, L, D, C, unitary_rep;


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
    C := L * DiagonalMat(List(D, Sqrt));

    return rec(basis_change := C,
               unitary_rep := ComposeHomFunction(rho, m -> C^-1 * m * C));
end );

# finds a nonscalar matrix in the centraliser of rho. This method is
# due to Dixon and involves summing over G, so may not be the
# fastest. Also it only works for unitary matrices.
FindCentraliser@ := function(rho)
    local G, n, A, H, r, s;
    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    A := function(r, s)
        local B;
        B := NullMat(n, n);
        B[r][s] := 1;
        return B;
    end;

    H := function(r, s)
            if r = s then
                return A(r,r);
            elif r > s then
                return A(r,s) + A(s,r);
            else
                return E(4) * (A(r,s) - A(s,r));
            fi;
    end;

    # find a non-scalar H in the centraliser
    for r in [1..n] do
        for s in [1..n] do
            H := 1/Order(G) * Sum(G, g -> ConjugateTranspose@(Image(rho, g)) * H(r,s) * Image(rho, g));

            # if it's nonscalar, return it!
            if H[1][1] * IdentityMat(n) <> H then
                return H;
            fi;
        od;
    od;

    return fail;
end;

InstallGlobalFunction( IrreducibleDecompositionDixon, function(rho)
    local G, n, H;

    G := Source(rho);
    n := DegreeOfRepresentation(rho);

    if not IsUnitaryRepresentation(rho) then
        Error("<rho> is not unitary!");
    fi;

    H := FindCentraliser@(rho);

    if H = fail then
        # this means rho is irreducible, only scalar matrices commute
        # with it, so there's only one component to the decomposition
        return [Cyclotomics^n];
    fi;

    # H is hermitian, so diagonalisable. we can find the jordan
    # decomp. of H then orthonormalise it

    # TODO: implement this maybe
    Error("not implemented yet");

end );
