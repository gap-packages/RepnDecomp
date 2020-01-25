# Returns list of first n elems from list
Take@ := function(list, n)
    return List([1..n], i -> list[i]);
end;

# Returns list of all but first n elems from list
Drop@ := function(list, n)
    return List([n+1..Length(list)], i -> list[i]);
end;

# Returns a list consisting of n copies of elem
Replicate@ := function(elem, n)
    if IsList(elem) then
        return List([1..n], i -> ShallowCopy(elem));
    else
        return List([1..n], i -> elem);
    fi;
end;

# Takes a list of blocks (possibly different sizes) and constructs a
# block diagonal matrix with those blocks.
BlockDiagonalMatrix@ := function(blocks)
    local combine_blocks, result, block;

    # Combines two blocks into a block diagonal matrix
    combine_blocks := function(b1, b2)
        local len1, len2, new_b1, new_b2;
        len1 := Length(b1);
        len2 := Length(b2);

        # Add len2 zeroes to the end of each row in b1
        new_b1 := List(b1, row -> Concatenation(row, Replicate@(0, len2)));

        # Add len1 zeroes to the start of each row in b2
        new_b2 := List(b2, row -> Concatenation(Replicate@(0, len1), row));

        return Concatenation(new_b1, new_b2);
    end;

    result := [];

    for block in blocks do
        result := combine_blocks(result, block);
    od;

    return result;
end;

# Takes 2 lists of reps of G and H and gives a full list of reps of
# GxH you can get from tensor products of them
InstallGlobalFunction( TensorProductRepLists, function(reps1, reps2)
    local G, H, GxH, gens, imgs, pi1, pi2, result, rep1, rep2, new_rep;
    G := Source(reps1[1]);
    H := Source(reps2[1]);
    GxH := DirectProduct(G, H);
    gens := GeneratorsOfGroup(GxH);
    pi1 := Projection(GxH, 1);
    pi2 := Projection(GxH, 2);

    result := [];
    for rep1 in reps1 do
        for rep2 in reps2 do
            imgs := List(gens, gxh -> KroneckerProduct(Image(rep1, Image(pi1, gxh)),
                                                       Image(rep2, Image(pi2, gxh))));
            new_rep := GroupHomomorphismByImagesNC(GxH, Group(imgs), gens, imgs);
            Add(result, new_rep);
        od;
    od;
    return result;
end );

# Returns a homomorphism g s.t. g(x) = func(hom(x))
InstallGlobalFunction( ComposeHomFunction, function(hom, func)
    local G, gens, imgs, H;
    G := Source(hom);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, gen -> func(Image(hom, gen)));
    H := Group(imgs);
    return GroupHomomorphismByImagesNC(G, H, gens, imgs);
end );

# Returns homomorphism from G to f(G)
FuncToHom@ := function(G, f)
    return ComposeHomFunction(IdentityMapping(G), f);
end;

# Returns a block diagonal representation isomorphic to the direct sum
# of the list of reps
InstallGlobalFunction( DirectSumOfRepresentations, function(reps)
    local G, gens, imgs, H;
    G := Source(reps[1]);
    return FuncToHom@(G, g -> BlockDiagonalMatrix@(List(reps, rep -> Image(rep, g))));
end );

# Takes the inner product of two characters, given as rows of the
# character table
InnerProductOfCharacters@ := function(chi1, chi2, G)
    local classes;

    # We avoid summing over the whole group
    classes := ConjugacyClasses(G);

    return (1/Size(G)) * Sum(List([1..Size(classes)],
                                  i -> Size(classes[i]) * chi1[i] * ComplexConjugate(chi2[i])));
end;

# Gives the row of the character table corresponding to irrep
CharacterOfRepresentation@ := function(irrep)
    local G;
    G := Source(irrep);
    return List(ConjugacyClasses(G),
                class -> Trace(Image(irrep, Representative(class))));
end;

# Irr(G), but guaranteed to be ordered the same as
# IrreducibleRepresentationsDixon (or the list of irreps given)
IrrWithCorrectOrdering@ := function(G)
    local irreps;
    irreps := ValueOption("irreps");

    if irreps = fail then
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    return List(irreps, irrep -> CharacterOfRepresentation@(irrep));
end;

# Writes the character of rho as a vector in the basis given by the
# irreducible characters (if chars are not given, use Dixon's method)
IrrVectorOfRepresentation@ := function(rho, irr_chars)
    local G, classes, char_rho, char_rho_basis;

    G := Source(rho);

    # Otherwise, we just compute using characters
    classes := ConjugacyClasses(G);

    char_rho := List(classes, class -> Trace(Image(rho, Representative(class))));

    # Write char_rho in the irr_chars basis for class functions
    char_rho_basis := List(irr_chars,
                           irr_char -> InnerProductOfCharacters@(char_rho, irr_char, G));

    return char_rho_basis;
end;

InstallGlobalFunction( DegreeOfRepresentation, function(rep)
    if IsPermGroup(Range(rep)) then
        return LargestMovedPoint(Range(rep));
    else
        return Trace(Image(rep, One(Source(rep))));
    fi;
end );

# Restricts a matrix to a given space. Resulting rep is given in the
# given basis. This only makes sense if V = span(basis) is
# G-invariant, otherwise results will be nonsense.
RestrictRep@ := function(rho, basis)
    local G, restricted_rho;

    G := Source(rho);

    restricted_rho := function(g)
        local imgs;

        # where g sends the basis vectors (row vects)
        imgs := List(basis, v -> Coefficients(basis, Image(rho, g) * v));

        # transpose to get column vectors
        return TransposedMat(imgs);
    end;

    return FuncToHom@(G, restricted_rho);
end;

# Calculates orbital matrix for G given representative
OrbitalMatrix@ := function(G, representative)
    local n, orbit, orbital, pair;

    # representative is assumed to be a pair [i, j] where i, j <= n

    if not IsPermGroup(G) then
        Error("G is not a permutation group!");
    fi;

    n := LargestMovedPoint(G);

    # We apply all elements of G to [i,j] and mark all reached
    # points with a 1 in the n x n zero matrix.

    # The action of G on pairs is just g(i, j) = (gi, gj)
    # i.e. OnTuples
    orbit := Orbit(G, representative, OnTuples);

    orbital := NullMat(n, n);

    for pair in orbit do
        orbital[pair[1]][pair[2]] := 1;
    od;

    return orbital;
end;

InstallGlobalFunction( PermToLinearRep, function(rho)
    local n;

    if not IsPermGroup(Range(rho)) then
        return fail;
    fi;

    n := DegreeOfRepresentation(rho);

    return ComposeHomFunction(rho, perm -> PermutationMat(perm, n));
end );

# Check a subspace is really G-invariant (action given by rho)
IsGInvariant@ := function(rho, in_space)
    local G, v, g, space;

    space := in_space;

    G := Source(rho);

    if not IsVectorSpace(space) then
        space := VectorSpace(Cyclotomics, space);
    fi;

    for v in Basis(space) do
        for g in GeneratorsOfGroup(G) do
            if not Image(rho, g) * v in space then
                return false;
            fi;
        od;
    od;
    return true;
end;

ConjugateTranspose@ := function(mat)
    return TransposedMat(List(mat, row -> List(row, ComplexConjugate)));
end;

# This is the matrix inner product used throughout the package
InnerProduct@ := function(A, B)
    return Trace(A*ConjugateTranspose@(B));
end;
# v is a list of matrices that are the basis of a space
OrthonormalBasis@ := function(v)
    local prod, proj, N, u, e, k;

    prod := InnerProduct@;

    proj := function(u, v)
        return (prod(u, v) / prod(u, u)) * u;
    end;

    N := Length(v);

    # these lists will be replaced so e[k] is the kth orthonormal basis vector
    u := [1..N];
    e := [1..N];

    for k in [1..N] do
        u[k] := v[k] - Sum([1..k-1], j -> proj(u[j], v[k]));
        e[k] := (1/Sqrt(prod(u[k], u[k]))) * u[k];
    od;

    return e;
end;

InstallGlobalFunction( IsOrthonormalSet, function(S, prod)
    return ForAll(S, v1 -> ForAll(S, function(v2)
                                     if v1 = v2 then
                                         return prod(v1, v2) = 1;
                                     else
                                         return prod(v1, v2) = 0;
                                     fi;
                                 end ));
end );

KroneckerList@ := function(reps)
    local prod, G;
    prod := function(mats)
        local current, mat;
        current := [[1]];
        for mat in mats do
            current := KroneckerProduct(current, mat);
        od;
        return current;
    end;
    G := Source(reps[1]);
    return FuncToHom@(G, g -> prod(List(reps, rep -> Image(rep, g))));
end;
