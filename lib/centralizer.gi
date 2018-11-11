# Functions to calculate the centralizer of a representation and its
# various decompositions.

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

# Takes the inner product of two characters, given as rows of the
# character table
CharacterInnerProduct@ := function(chi1, chi2, G)
    local classes;

    # We avoid summing over the whole group
    classes := ConjugacyClasses(G);

    return (1/Size(G)) * Sum(List([1..Size(classes)],
                                  i -> Size(classes[i]) * chi1[i] * ComplexConjugate(chi2[i])));
end;

# Writes the character of rho as a vector in the basis given by the
# irreducible characters
# TODO: Check if we can rely on the ordering
DecomposeCharacter@ := function(rho)
    local G, classes, irr_chars, char_rho, char_rho_basis;

    G := Source(rho);

    # Otherwise, we just compute using characters
    classes := ConjugacyClasses(G);
    irr_chars := Irr(G);
    char_rho := List(classes, class -> Trace(Image(rho, Representative(class))));

    # Write char_rho in the irr_chars basis for class functions
    char_rho_basis := List(irr_chars,
                           irr_char -> CharacterInnerProduct@(char_rho, irr_char, G));

    return char_rho_basis;
end;

# Computes the centralizer C of rho, returning generators of C as
# lists of blocks
InstallGlobalFunction( RepresentationCentralizerBlocks, function(orig_rho)
    local decomp, irrep_lists, rho, possible_blocks, zero_blocks, make_full_matrices, std_gens, classes, irr_chars, char_rho, char_rho_basis, all_sizes, sizes, G;

    rho := ConvertRhoIfNeeded@(orig_rho);
    G := Source(rho);

    irr_chars := Irr(G);
    char_rho_basis := DecomposeCharacter@(rho);

    # Calculate sizes based on the fact irr_char[1] is the degree
    all_sizes := List([1..Size(irr_chars)],
                      i -> rec(dimension := irr_chars[i][1],
                               nblocks := char_rho_basis[i]));

    # Now we remove all of the ones with nblocks = 0
    sizes := Filtered(all_sizes, r -> r.nblocks > 0);

    # We sort so we get the same ordering as Serre's formula's decomposition
    # TODO: Make sure this is actually true!!
    SortBy(sizes, r -> r.dimension);

    # There are two "levels" of blocks. First, the blocks
    # corresponding to each irreducible individually. Second, the
    # blocks that are the isomorphic blocks all grouped together.
    #
    # The centralizer only preserves the second type of block,
    # elements of C are block diagonal only according to the second
    # (larger) blocks.
    #
    # To work out the standard generators, we only need to know the
    # block sizes and collect together the isomorphic blocks, this is
    # what is calculated in "sizes".
    #
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

# This does the same as the previous, but uses 1x1 identity blocks
# always. This is the minimal dimension matrix presentation of the
# centraliser of rho.
InstallGlobalFunction( RepresentationCentralizerDecomposed, function(orig_rho)
    local decomp, irrep_lists, rho, possible_blocks, zero_blocks, make_full_matrices, std_gens, classes, irr_chars, char_rho, char_rho_basis, all_sizes, sizes, G;

    rho := ConvertRhoIfNeeded@(orig_rho);
    G := Source(rho);
    irr_chars := Irr(G);
    char_rho_basis := DecomposeCharacter@(rho);
    all_sizes := List([1..Size(irr_chars)],
                      i -> rec(dimension := irr_chars[i][1],
                               nblocks := char_rho_basis[i]));
    sizes := Filtered(all_sizes, r -> r.nblocks > 0);
    SortBy(sizes, r -> r.dimension);
    possible_blocks := List(sizes, size -> GenerateAllBlocks@(1, size.nblocks));
    zero_blocks := List(sizes, size -> NullMat(size.nblocks,
                                               size.nblocks));
    std_gens := Concatenation(List([1..Length(possible_blocks)],
                                   i -> ReplaceBlocks@(i, possible_blocks[i], zero_blocks)));
    return std_gens;
end );

# Same as DecomposeCentralizerBlocks but converts to full matrices
InstallGlobalFunction( RepresentationCentralizer, function(rho)
    return List(RepresentationCentralizerBlocks(rho), BlockDiagonalMatrix);
end );
