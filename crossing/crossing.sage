from sage.combinat.symmetric_group_representations import *
import sys, time

# Prints a sage matrix as a list of lists
def to_list_list(mat):
    return [list(row) for row in list(mat)]

# Need to define function to convert from Sage representation to GAP
# representation, i.e. from a SpechtRepresentation to a GAP
# homomorphism. This gives a string you can pass to libgap()
def to_gap_homo(m, irrep):
    G = SymmetricGroup(m)
    gens = G.gens()
    imgs = [irrep.representation_matrix(Permutation(g)) for g in gens]
    imgs = [to_list_list(mat) for mat in imgs]
    H = "Group({})".format(imgs)
    return "GroupHomomorphismByImagesNC(SymmetricGroup({}), {}, {}, {})".format(
        m,
        H,
        str(gens),
        str(imgs))

# Converts the list of irreps of S_m into a GAP string
def gap_all_irreps(m):
    irreps = SymmetricGroupRepresentations(m)
    gapped = [to_gap_homo(m, irrep) for irrep in irreps]
    return "[{}]".format(", ".join(gapped))

# Converts a matrix to one over RR
def real_mat(mat):
    return mat.apply_map(lambda x: x.n().real(), RR)

# This computes alpha_m using the full matrices with NO
# reductions. The dimensions are huge, but there is as little of my
# code involved as possible, so the answer is more likely to be
# correct
def compute_alpha_slow(m):
    Q = matrix(libgap.eval("Qmatrix({}).matrix".format(m)).sage())
    n = Q.ncols()
    J = matrix(n, n, lambda i,j: 1)

    # need to install this from GitHub https://github.com/mghasemi/pycsdp
    from SDP import SemidefiniteProgram

    # I have cvxopt here for now, need to try CSDP
    prog = SemidefiniteProgram(primal=False, solver="cvxopt")

    # X is already symmetric. We want all X_ij >= 0, so for each
    # choice of two indices (i, j) (there are n choose 2 choices), we
    # impose X_ij + X_ji >= 0 by adding a 1x1 block to X and imposing
    # X_kk = X_ij + X_ji where k is the index of the new block.
    #
    # We do this by imposing tr(AX) = 0 where 1 = A_ij = A_ji = -A_kk
    # and all other entries zero.
    pairs = Combinations(range(n), 2).list()

    # the new dimension of the matrices
    new_n = n + len(pairs)

    for k in range(len(pairs)):
        # X_kk is where the new block will be
        [i, j] = pairs[k]
        A = matrix(new_n)
        A[i, j] = 1
        A[j, i] = 1
        A[k, k] = -1
        prog.add_constraint(A, 0)

    # we pad Q and J with n choose 2 1x1 zero blocks
    pad_Q = matrix(new_n, new_n, lambda i, j: Q[i][j] if i<n and j<n else 0)
    pad_J = matrix(new_n, new_n, lambda i, j: J[i][j] if i<n and j<n else 0)

    # the SDP solver always assumes we want to maximise tr(QX) but we
    # really want to minimise it, so we negate
    prog.set_objective(-1*pad_Q)

    # tr(JX) = 1
    prog.add_constraint(pad_J, 1)


    soln = prog.solve()

    print("alpha_{} = {}".format(m, -1*soln["DObj"]))


# Computes the nice basis to use (that block diagonalises the
# centralizer ring), using my GAP package
def compute_alpha(m, print_irreps=False, status=True):
    if status:
        sys.stdout.write("Calculating irreps of S_{} x S_2 in Sage: ".format(m))

    t0 = time.time()
    irreps_s_m = libgap.eval("irreps_s_m := {};".format(gap_all_irreps(m)))
    t1 = time.time()

    if status:
        print(t1-t0)

    # the irreps of s_2
    irreps_s_2 = libgap.eval("irreps_s_2 := [ GroupHomomorphismByImagesNC( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ 1 ] ] ]), [ (1,2) ], [ [ [ 1 ] ] ] ), GroupHomomorphismByImagesNC( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ -1 ] ] ]), [ (1,2) ], [ [ [ -1 ] ] ] ) ];")

    # list of irreps of the whole group
    irreps_G = libgap.eval("irreps_G := TensorProductRepLists(irreps_s_m, irreps_s_2);")

    if print_irreps:
        print(libgap.Display(irreps_G))

    # Calculates the matrices needed for the (reduced version of the)
    # semidefinite program
    if status:
        sys.stdout.write("Block diagonalizing SDP in GAP: ")

    t0 = time.time()
    libgap.eval("sdp := CalculateSDP({}, irreps_G);".format(m))
    t1 = time.time()

    if status:
        print(t1-t0)

    if status:
        sys.stdout.write("Formulating SDP in Sage\n")

    # E_i are the orbital matrices, these are orthogonal
    E = [matrix(mat) for mat in libgap.eval("List(sdp.centralizer_basis);").sage()]

    # B_i are the block diagonalised, normalised E_i
    B = [matrix(mat) for mat in libgap.eval("List(sdp.nice_cent_basis);").sage()]

    # (L_k)_ij = lambda^i_{kj}
    L = [matrix(mat) for mat in libgap.eval("sdp.param_matrices;").sage()]

    # pairs[i] = i^* in GAP, need to adjust indices for Sage (GAP is
    # 1-based)
    pair_map = [y-1 for y in libgap.eval("sdp.pairs").sage()]

    # basis the B are written in, the action_hom is written in the
    # standard basis
    basis_change = matrix(libgap.eval("sdp.nice_change;").sage())
    basis_change_inv = matrix(libgap.eval("sdp.nice_change_inv;").sage())

    # TODO: work out where we need to change basis to make this fast
    # the cost matrix Q_xy is the number of adjacent transpositions needed to get from x to y^-1
    if status:
        sys.stdout.write("Calculating Q: ")

    t0 = time.time()
    Q = matrix(libgap.eval("Qmatrix({}).matrix".format(m)).sage())
    t1 = time.time()

    if status:
        print(t1-t0)

    # from this point, we don't need GAP

    if status:
        sys.stdout.write("Calculating block diagonal matrices: ")

    t0 = time.time()
    Q = basis_change_inv * Q * basis_change
    J = matrix(basis_change.nrows(), basis_change.nrows(), lambda i,j: 1)
    J = basis_change_inv * J * basis_change
    t1 = time.time()

    if status:
        print(t1-t0)

    # Our matrices are all real, but largely irrational. This is ok
    # except the irrationals are represented as sums of
    # cyclotomics. Sage doesn't know they are real. We have to
    # truncate the imaginary parts.

    #if status:
    #    sys.stdout.write("Truncating imaginary parts: ")

    #t0 = time.time()
    #B = [real_mat(mat) for mat in B]
    #Q = real_mat(Q)
    #J = real_mat(J)
    #L = [real_mat(mat) for mat in L]
    #t1 = time.time()

    #if status:
    #    print(t1-t0)

    # See the paper https://homepages.cwi.nl/~lex/files/symm.pdf for
    # some explanation of this program
    if status:
        print("Setting up SDP in Sage")
    prog = SemidefiniteProgram(maximization=False)
    x = prog.new_variable()

    # The number of variables we need is equal to the number of {i,
    # i^*} pairs.
    pairs = []
    seen = set()
    for i in range(len(pair_map)):
        if i not in seen:
            seen.add(i)
            seen.add(pair_map[i])
            pairs.append((i, pair_map[i]))

    # verify the pairs
    if status:
        print("All indices have been paired: {}".format(all(x in seen for x in range(len(B)))))
        print("Pairs are correct: {}".format(all(E[x] == E[y].transpose() for (x, y) in pairs)))
        print("d_centralizer = {}".format(len(B)))
        print("d_reduced = {}".format(len(pairs)))

    # The dimension of the subspace of the centralizer consisting of
    # symmetric matrices
    d = len(pairs)

    if status:
        sys.stdout.write("Adding constraints: ")

    t0 = time.time()

    # the x[j] variable we use for a pair (i, i*) is the minimum since
    # the pairs cover [1..d'] where d' is the full dimension of the
    # centralizer
    prog.set_objective(sum((Q * (B[i] + B[pair_map[i]])).trace()*x[i] for i in range(d)))

    # sum_i x_i L_i >= 0. This is true iff sum_i x_i B_i >= 0 due to a
    # *-isomorphism B_i -> L_i
    constraint0 = sum(x[i] * (L[i] + L[pair_map[i]]) for i in range(d)) >= 0

    one = matrix([[1]])

    constraint1 = sum((J * (B[i] + B[pair_map[i]])).trace()*x[i] for i in range(d))*one == one

    prog.add_constraint(constraint0)
    prog.add_constraint(constraint1)

    # need all x_i >= 0, this specifies X >= 0 where X_ii = x_i,
    # otherwise 0, which is the same thing
    C = [matrix() for i in range(d)]
    for i in range(d):
        diag = [0]*d
        diag[i] = 1
        C[i] = diagonal_matrix(diag)

    prog.add_constraint(sum(C[i]*x[i] for i in range(d)) >= 0)

    t1 = time.time()

    if status:
        print(t1-t0)

    if status:
        sys.stdout.write("Solving SDP: ")

    t0 = time.time()
    alpha_m = prog.solve()
    t1 = time.time()

    print(t1-t0)
    print("alpha_{} = {}".format(m, alpha_m))

# Make sure your gap_cmd etc are set up properly so that my package is
# available
libgap.eval('LoadPackage("RepnDecomp");')

# read the file that computes action_hom, the regular
# representation of the group action of G on the m cycles
libgap.eval('Read("crossing.g");')
