#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# Implementations
#

# Calculates p(V) for p a linear map (given as a matrix in the
# standard basis) and a vector space V
MatrixImage@ := function(p, V)
    local F;

    # F is the base field of V
    F := LeftActingDomain(V);

    # We return the span of the images of the basis under p, which
    # gives p(V)
    return VectorSpace(F,
                       List(Basis(V), v -> p * v),
                       Zero(V));
end;

InstallGlobalFunction( DecomposeRepresentationCanonical, function(rho, arg...)
    local G, F, n, V, irreps, chars, char_to_proj, canonical_projections, canonical_summands;

    # The group we are taking representations of
    G := Source(rho);

    # The field we are working over: it's always the Cyclotomics
    F := Cyclotomics;

    # The dimension of the V in rho : G -> GL(V). Since we have the
    # images of rho as matrices, this is just the width or height of
    # any image of any generator of G.
    n := Length(Range(rho).1);

    # The vector space that the linear maps act on
    V := F^n;

    # The full list of irreps W_i of G over F
    #
    # If we are given a list of irreps, we use that instead of
    # calculating it
    if Size(arg) > 0 then
        irreps := arg[1];
    else
        irreps := IrreducibleRepresentations(G, F);
    fi;

    return List(irreps, function (irrep)
                   local character, degree, projection, canonical_summand;

                   # In Serre's text, irrep is called W_i, this character is chi_i
                   character := GroupHomomorphismByFunction(G, F, g -> Trace(Image(irrep, g)));
                   degree := Image(character, One(G));

                   # Calculate the projection map from V to irrep using Theorem 8 (Serre)
                   # Given as a matrix, p_i
                   projection := (degree/Order(G)) * Sum(G, t -> ComplexConjugate(Image(character, t)) * Image(rho, t));

                   # Calculate V_i, the canonical summand
                   canonical_summand := MatrixImage@(projection, V);
                   return canonical_summand;
               end );
end );

# Decomposes the representation V_i into a direct sum of some number
# (maybe zero) of spaces, all isomorphic to W_i. W_i is the space
# corresponding to the irrep : G -> GL(W_i). rho is the "full"
# representation that we're decomposing.
DecomposeCanonicalSummand@ := function(rho, irrep, V_i)
    local projections, p_11, V_i1, basis, n, step_c, G, H, F, V, m;

    G := Source(irrep);

    # This is the general linear group of some space, we don't really
    # know or care what the space actually is
    H := Range(irrep);

    # This gives the dimension of the space of which W is the general
    # linear group (the size of the matrices representing the maps)
    n := Length(H.1);

    m := Length(Range(rho).1);
    F := Cyclotomics;
    V := F^m;

    # First compute the projections p_ab. We only actually use projections with
    # a=1..n and b=1, so we can just compute those. projections[a] is p_{a1}
    # from Serre.
    projections := List([1..n], a -> (n/Order(G)) * Sum(G, t -> Image(irrep,t^-1)[1][a]*Image(rho,t)));

    p_11 := projections[1];
    V_i1 := MatrixImage@(p_11, V_i);
    basis := Basis(V_i1);

    # Now we define the map taking x to W(x), a subrepresentation of
    # V_i isomorphic to W_i. (This is step (c) of Proposition 8)
    step_c := function(x1)
        # This is the list of basis vectors for W(x1)
        return List([1..n], alpha -> projections[alpha] * x1);
    end;

    # If x1^1 .. x1^m is a basis for V_i1 (this is in the `basis`
    # variable), then V_i decomposes into the direct sum W(x1^1)
    # ... W(x1^m), each isomorphic to W_i.
    #
    # We return a list of lists of (vector space, basis) pairs where
    # the basis (TODO: confirm this?) has the special property
    return List(basis, function(x)
                   local b;
                   b := step_c(x);
                   return rec(space := VectorSpace(F, b, Zero(V)), basis := b);
               end);
end;

# Converts rho to a matrix representation if necessary
ConvertRhoIfNeeded@ := function(rho)
    local gens, ims, high, new_ims, new_range, new_rho, G;
    G := Source(rho);
    # We want rho to be a homomorphism to a matrix group since this
    # algorithm works on matrices. We convert a permutation group into
    # an isomorphic matrix group so that this is the case. If we don't
    # know how to convert to a matrix group, we just fail.
    new_rho := rho;
    if not IsMatrixGroup(Range(rho)) then
        if IsPermGroup(Range(rho)) then
            gens := GeneratorsOfGroup(G);
            ims := List(gens, g -> Image(rho, g));
            high := LargestMovedPoint(ims);
            new_ims := List(ims, i -> PermutationMat(i, high));
            new_range := Group(new_ims);
            new_rho := GroupHomomorphismByImages(G, new_range, gens, new_ims);
        else
            Error("rho is not a matrix or permutation group!");
        fi;
    fi;
    return new_rho;
end;

# Decompose rho into irreducible representations with the reps that
# are isomorphic collected together. This returns a list of lists of
# vector spaces (L) with each element of L being a list of vector
# spaces arising from the same irreducible.
DecomposeIsomorphicCollected@ := function(orig_rho, arg...)
    local irreps, N, canonical_summands, full_decomposition, G, F, n, V, gens, ims, high, new_ims, new_range, rho;

    rho := ConvertRhoIfNeeded@(orig_rho);
    G := Source(rho);

    F := Cyclotomics;
    n := Length(Range(rho).1);
    V := F^n;

    # Use the list of irreps if given
    if Size(arg) > 0 then
        irreps := arg[1];
    else
        irreps := IrreducibleRepresentations(G, F);
    fi;

    N := Size(irreps);

    # This gives a list of vector spaces, each a canonical summand
    # Try to use a precomputed decomposition, if given
    if Size(arg) > 1 then
        canonical_summands := arg[2];
    else
        canonical_summands := DecomposeRepresentationCanonical(rho, irreps);
    fi;

    # This gives a list of lists of vector spaces, each a
    # decomposition of a canonical summand into irreducibles.
    full_decomposition := List([1..N],
                               i -> DecomposeCanonicalSummand@(rho, irreps[i], canonical_summands[i]));

    # Here we return the rho we actually used i.e. after we convert to
    # an isomorphic rep that goes to a matrix group (not a permutation
    # group)
    return rec(decomp := full_decomposition, used_rho := rho);
end;

# Takes a rho that goes to a matrix group only. Returns a basis change
# matrix which, when used on a given rho(g) (matrix), block
# diagonalises rho(g) such that each block corresponds to an irrep.
BaseChangeMatrix@ := function(rho, decomp)
    local new_bases, new_basis;

    # Extract the basis vectors, this is now a list of lists of bases
    # (each basis is a list of vectors)
    new_bases := List(decomp,
                      rec_list -> List(rec_list, r -> r.basis));

    # List of new basis row vectors
    new_basis := Concatenation(Concatenation(new_bases));

    # Base change matrix from new basis to standard basis
    return TransposedMat(new_basis);
end;

# Takes a representation going to a matrix group and gives you an
# isomorphic representation where the images are block-diagonal with
# each block corresponding to an irreducible representation
InstallGlobalFunction( BlockDiagonalizeRepresentation, function(rho, arg...)
    local decomp, A, G, gens, imgs, range;

    # Use precomputed decomposition, if available
    if Size(arg) > 0 then
        decomp := arg[1];
    else
        decomp := DecomposeIsomorphicCollected@(rho);
    fi;

    A := BaseChangeMatrix@(rho, decomp.decomp);
    G := Source(rho);
    gens := GeneratorsOfGroup(G);
    imgs := List(gens, g -> A^-1 * Image(decomp.used_rho, g) * A);

    range := Group(imgs);

    return GroupHomomorphismByImages(G, range, gens, imgs);
end );

# Gives the list of vector spaces in the direct sum
# decomposition of rho : G -> GL(V) into irreducibles.
InstallGlobalFunction( DecomposeRepresentationIrreducible, function(orig_rho, arg...)
    local irreps, canonical_summands, rho;
    rho := ConvertRhoIfNeeded@(orig_rho);

    # Use the list of irreps if given
    if Size(arg) > 0 then
        irreps := arg[1];
    else
        irreps := IrreducibleRepresentations(Source(rho), Cyclotomics);
    fi;

    # Try to use a precomputed decomposition, if given
    if Size(arg) > 1 then
        canonical_summands := arg[2];
    else
        canonical_summands := DecomposeRepresentationCanonical(rho, irreps);
    fi;

    # We only want to return the vector spaces here
    return Flat(List(DecomposeIsomorphicCollected@(rho, irreps, canonical_summands).decomp,
                     rec_list -> List(rec_list, r -> r.space)));
end );

# Returns list of first n elems from list
Take@ := function(list, n)
    local result, count;
    result := [];
    count := 0;
    while count < n do
        Add(result, list[count+1]);
        count := count + 1;
    od;
    return result;
end;

# Returns list of all but first n elems from list
Drop@ := function(list, n)
    local result, count, elem;
    result := [];
    count := 0;
    for elem in list do
        if count >= n then
            Add(result, elem);
        fi;
        count := count + 1;
    od;
    return result;
end;

# Returns a list consisting of n copies of elem
Replicate@ := function(elem, n)
    local result, i;
    result := [];
    for i in [1..n] do
        Add(result, elem);
    od;
    return result;
end;

# Decomposes a block-diagonal matrix into a list of blocks given a
# list of block sizes
#
# WARNING: It's not checked whether the bits thrown away were actually
# all zero, it is just assumed. Make sure the block sizes are correct.
DecomposeMatrixIntoBlocks@ := function(matrix, block_sizes)
    local blocks, my_matrix, block, new_matrix, block_size;

    if Length(matrix) <> Sum(block_sizes) then
        Error("block sizes don't match matrix size");
    fi;

    blocks := [];
    my_matrix := matrix;

    for block_size in block_sizes do
        # First cut the matrix in two, with the desired block above
        # and the rest below
        block := Take@(my_matrix, block_size);
        new_matrix := Drop@(my_matrix, block_size);

        # Then cut off the zero parts
        block := List(block, row -> Take@(row, block_size));
        new_matrix := List(new_matrix, row -> Drop@(row, block_size));

        Add(blocks, block);
        my_matrix := new_matrix;
    od;

    return blocks;
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

# Takes a size of small block (dimension) and size of grid of small
# blocks (nblocks).
#
# Returns a list of big matrices (square of size dimension*nblocks)
# each with exactly 1 small block nonzero, equal to the identity, in
# all possible ways. Result will be nblocks^2 long.
GenerateAllBlocks@ := function(dimension, nblocks)
    local result, coords, coord, i, j, gen;
    result := [];

    # Possible locations of the I block
    coords := Cartesian([1..nblocks], [1..nblocks]);

    for coord in coords do
        i := coord[1];
        j := coord[2];

        # a single block at position (i,j)
        gen := BlockMatrix([[i, j, IdentityMat(dimension)]], nblocks, nblocks);

        Add(result, MatrixByBlockMatrix(gen));
    od;

    return result;
end;

# Returns a list of block matrices where zero_blocks[i]
# has been replaced with each block in blocks
ReplaceBlocks@ := function(i, blocks, zero_blocks)
    local result, block, full_matrix_blocks;
    result := [];
    for block in blocks do
        full_matrix_blocks := ShallowCopy(zero_blocks);
        full_matrix_blocks[i] := block;
        Add(result, full_matrix_blocks);
    od;
    return result;
end;

# Computes the centralizer C of rho, returning generators of C as
# lists of blocks
InstallGlobalFunction( RepresentationCentralizerBlocks, function(rho, arg...)
    local decomp, irrep_lists, used_rho, sizes, possible_blocks, zero_blocks, make_full_matrices, std_gens;

    if Size(arg) > 0 then
        decomp := arg[1];
    else
        decomp := DecomposeIsomorphicCollected@(rho);
    fi;

    irrep_lists := decomp.decomp;
    used_rho := decomp.used_rho;

    # There are two "levels" of blocks. First, the blocks
    # corresponding to each irreducible individually. Second, the
    # blocks that are the isomorphic blocks all grouped together.
    #
    # The centralizer only preserves the second type of block,
    # elements of C are block diagonal only according to the second
    # (larger) blocks.
    #
    # To work out the standard generators, we only need to know the
    # block sizes and collect together the isomorphic blocks.

    sizes := List(irrep_lists,
                  irrep_list -> rec(dimension := Dimension(irrep_list[1].space),
                                    nblocks := Length(irrep_list)));

    # If a list of isomorphic blocks is n long, it gives n^2 standard
    # generators, each with exactly 1 block, in the (i,j) position
    # with X_i isomorphic to X_j (the irreps they correspond to) and
    # the block equal to I_{dim X_i} (for all possible i and j).
    #
    # For each collection of isomorphic blocks, we want all possible
    # nonzero big blocks, a list of lists of blocks
    possible_blocks := List(sizes, size -> GenerateAllBlocks@(size.dimension, size.nblocks));

    # A list of correctly sized zero blocks. Big blocks, not
    # individual small blocks
    zero_blocks := List(sizes, size -> NullMat(size.dimension * size.nblocks,
                                               size.dimension * size.nblocks));

    # All standard generators
    std_gens := Concatenation(List([1..Length(possible_blocks)],
                                   i -> ReplaceBlocks@(i, possible_blocks[i], zero_blocks)));

    return std_gens;
end );

# This does the same as the previous, but uses 1x1 identity blocks always
InstallGlobalFunction( RepresentationCentralizerDecomposed, function(rho, arg...)
    local decomp, irrep_lists, used_rho, sizes, cmp_possible_blocks, cmp_zero_blocks, cmp_std_gens;

    if Size(arg) > 0 then
        decomp := arg[1];
    else
        decomp := DecomposeIsomorphicCollected@(rho);
    fi;

    irrep_lists := decomp.decomp;
    used_rho := decomp.used_rho;
    sizes := List(irrep_lists,
                  irrep_list -> rec(dimension := Dimension(irrep_list[1].space),
                                    nblocks := Length(irrep_list)));
    cmp_possible_blocks := List(sizes, size -> GenerateAllBlocks@(1, size.nblocks));
    cmp_zero_blocks := List(sizes, size -> NullMat(size.nblocks,
                                                   size.nblocks));
    cmp_std_gens := Concatenation(List([1..Length(cmp_possible_blocks)],
                                       i -> ReplaceBlocks@(i, cmp_possible_blocks[i], cmp_zero_blocks)));
    return cmp_std_gens;
end );

# Same as DecomposeCentralizerBlocks but converts to full matrices
InstallGlobalFunction( RepresentationCentralizer, function(rho, arg...)
    local decomp;
    if Size(arg) > 0 then
        decomp := arg[1];
    else
        decomp := DecomposeIsomorphicCollected@(rho);
    fi;
    return List(RepresentationCentralizerBlocks(rho, decomp), BlockDiagonalMatrix@);
end );
