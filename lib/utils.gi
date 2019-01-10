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
InstallGlobalFunction( BlockDiagonalMatrix, function(blocks)
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
end );

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
            new_rep := GroupHomomorphismByImages(GxH, Group(imgs), gens, imgs);
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
    return GroupHomomorphismByImages(G, H, gens, imgs);
end );

# Returns a block diagonal representation isomorphic to the direct sum
# of the list of reps
InstallGlobalFunction( DirectSumRepList, function(reps)
    local G, gens, imgs, H;
    G := Source(reps[1]);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, g -> BlockDiagonalMatrix(List(reps, rep -> Image(rep, g))));
    H := Group(imgs);
    return GroupHomomorphismByImages(G, H, gens, imgs);
end );

# Takes the inner product of two characters, given as rows of the
# character table
CharacterInnerProduct@ := function(chi1, chi2, G)
    local classes;

    # We avoid summing over the whole group
    classes := ConjugacyClasses(G);

    return (1/Size(G)) * Sum(List([1..Size(classes)],
                                  i -> Size(classes[i]) * chi1[i] * ComplexConjugate(chi2[i])));
end;

# Gives the row of the character table corresponding to irrep
IrrepToCharacter@ := function(irrep)
    local G;
    G := Source(irrep);
    return List(ConjugacyClasses(G),
                class -> Trace(Image(irrep, Representative(class))));
end;

# Irr(G), but guaranteed to be ordered the same as
# IrreducibleRepresentationsDixon (or the list of irreps given)
IrrWithCorrectOrdering@ := function(G, args...)
    local irreps;
    irreps := [];

    if Length(args) > 0 then
        irreps := args[1];
    else
        irreps := IrreducibleRepresentationsDixon(G);
    fi;

    return List(irreps, irrep -> IrrepToCharacter@(irrep));
end;

# Writes the character of rho as a vector in the basis given by the
# irreducible characters (if chars are not given, use Dixon's method)
DecomposeCharacter@ := function(rho, args...)
    local G, classes, irr_chars, char_rho, char_rho_basis;

    G := Source(rho);

    # Otherwise, we just compute using characters
    classes := ConjugacyClasses(G);

    # If we are given chars, just use those
    irr_chars := [];
    if Length(args) > 0 then
        irr_chars := args[1];
    else
        # We could use Irr(G) here, but we want to keep all ordering
        # consistent with IrreducibleRepresentations
        irr_chars := IrrWithCorrectOrdering@(G);
    fi;
    char_rho := List(classes, class -> Trace(Image(rho, Representative(class))));

    # Write char_rho in the irr_chars basis for class functions
    char_rho_basis := List(irr_chars,
                           irr_char -> CharacterInnerProduct@(char_rho, irr_char, G));

    return char_rho_basis;
end;

# Tells you if two representations of the same group are isomorphic by
# examining characters
InstallGlobalFunction( AreRepsIsomorphic, function(rep1, rep2)
    local G, irr_chars;

    if Source(rep1) <> Source(rep2) then
        return false;
    fi;

    G := Source(rep1);
    irr_chars := IrrWithCorrectOrdering@(G);

    # Writes the characters in the irr_chars basis, they are the same
    # iff they are isomorphic
    return DecomposeCharacter@(rep1, irr_chars) = DecomposeCharacter@(rep2, irr_chars);
end );
