from sage.combinat.symmetric_group_representations import *

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
    return "GroupHomomorphismByImages(SymmetricGroup({}), {}, {}, {})".format(
        m,
        H,
        str(gens),
        str(imgs))

# Converts the list of irreps of S_m into a GAP string
def gap_all_irreps(m):
    irreps = SymmetricGroupRepresentations(m)
    gapped = [to_gap_homo(m, irrep) for irrep in irreps]
    return "[{}]".format(", ".join(gapped))

# Computes the nice basis to use (that block diagonalises the
# centralizer ring), using my GAP package
def compute_alpha(m):
    #import pdb; pdb.set_trace()
    irreps_s_m = libgap.eval("irreps_s_m := {};".format(gap_all_irreps(m)))

    # the irreps of s_2
    irreps_s_2 = libgap.eval("irreps_s_2 := [ GroupHomomorphismByImages( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ 1 ] ] ]), [ (1,2) ], [ [ [ 1 ] ] ] ), GroupHomomorphismByImages( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ -1 ] ] ]), [ (1,2) ], [ [ [ -1 ] ] ] ) ];")

    # list of irreps of the whole group
    irreps_G = libgap.eval("irreps_G := TensorProductRepLists(irreps_s_m, irreps_s_2);")

    # read the file that computes action_hom, the regular
    # representation of the group action of G on the m cycles
    libgap.eval('Read("perms.g");')

    # Calculates the matrices needed for the (reduced version of the)
    # semidefinite program
    libgap.eval("record := CalculateSDP(irreps_G);")

    B = libgap.eval("record.centralizer_basis;").sage()
    basis = libgap.eval("record.nice_basis;").sage()
    Q = matrix(libgap.eval("Qmatrix({}).matrix".format(m)).sage())
    basis_change = matrix(libgap.eval("record.nice_basis").sage()).transpose()
    Q_block_diag = basis_change^-1 * Q * basis_change
    J = matrix(QQ, len(basis), len(basis), lambda i,j: 1)
    J_block_diag = basis_change^-1 * J * basis_change
    L = libgap.eval("record.param_matrices;").sage()

    # See the paper https://homepages.cwi.nl/~lex/files/symm.pdf for
    # some explanation of this program
    prog = SemidefiniteProgram()
    x = p.new_variable()

    d = len(L)

    prog.set_objective(sum((Q_block_diag * B[i]).trace()*x[i] for i in range(d)))
    prog.add_constraint(sum(x[i] * L[i] for i in range(d)) >= 0)
    prog.add_constraint(sum((J_block_diag * B[i]).trace()*x[i]) == 1)

    return prog.solve()

# Loading this file will compute alpha_5. Make sure your gap_cmd is
# set up properly so the RepnDecomp package is available
libgap.eval('LoadPackage("RepnDecomp");')
compute_alpha(3)
