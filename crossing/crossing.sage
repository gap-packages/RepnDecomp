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
    return matrix(RR, mat.nrows(), mat.ncols(), lambda i, j: mat[i][j].n().real())

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

    print("alpha_{} = {}".format(m, soln["DObj"]))


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
        sys.stdout.write("Formulating SDP in Sage: ")

    t0 = time.time()
    B = [matrix(mat) for mat in libgap.eval("sdp.centralizer_basis;").sage()]
    basis = libgap.eval("sdp.nice_basis;").sage()
    Q = matrix(libgap.eval("Qmatrix({}).matrix".format(m)).sage())
    basis_change = matrix(libgap.eval("TransposedMat(sdp.nice_basis)").sage())
    Q_block_diag = basis_change^-1 * Q * basis_change
    J = matrix(len(basis), len(basis), lambda i,j: 1)
    J_block_diag = basis_change^-1 * J * basis_change
    L = [matrix(mat) for mat in libgap.eval("sdp.param_matrices;").sage()]

    # Our matrices are all real, but largely irrational. This is ok
    # except the irrationals are represented as sums of
    # cyclotomics. Sage doesn't know they are real. We have to
    # truncate the imaginary parts.
    B = [real_mat(mat) for mat in B]
    Q_block_diag = real_mat(Q_block_diag)
    J_block_diag = real_mat(J_block_diag)
    L = [real_mat(mat) for mat in L]
    t1 = time.time()

    if status:
        print(t1-t0)

    #import pdb; pdb.set_trace()

    # See the paper https://homepages.cwi.nl/~lex/files/symm.pdf for
    # some explanation of this program
    if status:
        print("Setting up SDP in Sage")
    prog = SemidefiniteProgram(maximization=False)
    x = prog.new_variable()

    d = len(L)

    prog.set_objective(sum((Q_block_diag * B[i]).trace()*x[i] for i in range(d)))

    constraint0 = sum(x[i] * L[i] for i in range(d)) >= 0

    # This is what this means, but I need the matrices to avoid type errors
    # constraint1 = sum((J_block_diag * B[i]).trace()*x[i] for i in range(d)) == 1

    one = matrix([[1]])
    constraint1 = sum((J_block_diag * B[i]).trace()*x[i]*one for i in range(d)) == one

    if status:
        print("Adding constraints")
    prog.add_constraint(constraint0)
    prog.add_constraint(constraint1)

    # need all x_i >= 0
    for i in range(d):
        prog.add_constraint(x[i]*one >= 0)

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
