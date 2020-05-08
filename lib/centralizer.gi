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

    if dimension = 0 or nblocks = 0 then
        return [];
    fi;

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

# Converts a list of sizes of blocks in a decomposition to a basis for
# the centralizer ring of the representation
InstallGlobalFunction( SizesToBlocks, function(sizes)
    local possible_blocks, zero_blocks, std_gens;

    # Note: we don't sort by the dimension of the blocks or anything
    # since we want the block order to match with the order of
    # irr_chars

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

# Computes the centralizer C of rho, returning generators of C as
# lists of blocks.
# NOTE: This is written in the nice basis given by BlockDiagonalBasisOfRepresentation
InstallGlobalFunction( CentralizerBlocksOfRepresentation, function(rho)
    return ComputeUsingMethod@(rho).centralizer_basis;
end );

# Same as DecomposeCentralizerBlocks but converts to full matrices
InstallGlobalFunction( CentralizerOfRepresentation, function(rho)
    return List(CentralizerBlocksOfRepresentation(rho), BlockDiagonalMatrix@);
end );

# Given an orthonormal basis for the centralizer (w.r.t the product
# tr(AB^*)), a representative, and a conjugacy class, calculates the
# sum of the conjugacy class matrices. Result is given as a list of
# coefficients in the basis.
ClassSumCentralizerCoeffs@ := function(rho, class, cent_basis)
    local coeff, conj;

    coeff := B -> Size(class) * InnerProduct@(Image(rho, Representative(class)), B);

    return List(cent_basis, coeff);
end;

# Returns the actual class sum, i.e. after summing the cent_basis with
# the coefficients. Requires rho to be unitary and cent_basis to be
# orthonormal!
InstallGlobalFunction( ClassSumCentralizerNC, function(rho, class, cent_basis)
    local coeffs;

    if cent_basis <> fail then
        coeffs := ClassSumCentralizerCoeffs@(rho, class, cent_basis);
        return Sum([1..Length(coeffs)], i -> coeffs[i] * cent_basis[i]);
    else
        return Sum(class, g -> Image(rho, g));
    fi;
end );

InstallGlobalFunction( ClassSumCentralizer, function(rho, class, cent_basis)
    if not IsUnitaryRepresentation(rho) then
        Error("<rho> is not unitary!");
    fi;

    if cent_basis <> fail and not IsOrthonormalSet(cent_basis, InnerProduct@) then
        Error("<cent_basis> is not an orthonormal set!");
    fi;

    return ClassSumCentralizerNC(rho, class, cent_basis);
end );

# Does the group sum (just for convenience, nothing clever)
GroupSumWithCentralizer@ := function(rho, cent_basis)
    local G, classes;
    G := Source(rho);
    classes := ConjugacyClasses(G);
    if cent_basis <> fail then
        return Sum(classes, class -> ClassSumCentralizerNC(rho, class, cent_basis));
    else
        return Sum(G, g -> Image(rho, g));
    fi;
end;
