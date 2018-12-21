from sage.combinat.symmetric_group_representations import *

# Prints a sage matrix as a list of lists
def to_list_list(mat):
    return [list(row) for row in list(mat)]

# Need to define function to convert from Sage representation to GAP
# representation, i.e. from a SpechtRepresentation to a GAP
# homomorphism. This gives a string you can pass to gap()
def to_gap_homo(m, irrep):
    G = SymmetricGroup(m)
    gens = G.gens()
    imgs = [irrep.representation_matrix(Permutation(g)) for g in gens]
    imgs = [to_list_list(m) for m in imgs]
    H = "Group({})".format(imgs)
    return "GroupHomomorphismByImages({}, {}, {}, {})".format(
        str(gap(G)),
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
def compute_nice_basis(m):
    irreps_s_m = gap(gap_all_irreps(m))

    # the irreps of s_2
    irreps_s_2 = gap("[ GroupHomomorphismByImages( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ 1 ] ] ]), Pcgs([ (1,2) ]), [ [ [ 1 ] ] ] ), GroupHomomorphismByImages( SymmetricGroup( [ 1 .. 2 ] ), Group([ [ [ -1 ] ] ]), Pcgs([ (1,2) ]), [ [ [ -1 ] ] ] ) ]")

    # list of irreps of the whole group
    irreps_G = gap("TensorProductRepLists({}, {})".format(irreps_s_m.name(), irreps_s_2.name()))

    # The nice basis that block diagonalizes the centralizer with minimal sized blocks
    nice_basis = gap("BlockDiagonalRepresentationFast(action_hom, {}).basis".format(irreps_G.name()))

    return nice_basis.sage()

def compute_alpha(m):
    basis = compute_nice_basis(m)
    Q = gap("Qmatrix({})").format(m)
    basis_change = matrix(basis).transpose()
    Q_block_diag = basis_change^-1 * Q * basis_change
    J = matrix(QQ, len(basis), len(basis), lambda i,j: 1)

    prog = SemidefiniteProgram()
    X = p.new_variable()
    prog.set_objective((Q_block_diag * X).trace())
    prog.add_constraint(X >= 0)
    prog.add_constraint((J * X).trace() == 1)

    optX = prog.solve()
    return (Q_block_diag * optX).trace()
