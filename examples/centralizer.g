#! @BeginChunk Example_CentralizerBlocksOfRepresentation
#! @BeginExample
G := DihedralGroup(8);;
irreps := IrreducibleRepresentations(G);;
# rho is the sum of two isomorphic degree 1 irreps, and a degree
# 2 irrep.
rho := DirectSumOfRepresentations([irreps[4], irreps[4], irreps[5]]);;
# Compute a basis for the centralizer (in blocks)
cent_basis_blocks := CentralizerBlocksOfRepresentation(rho);;
# Verify that the dimension is the sum of the multiplicities squared,
# in this case 2^2 + 1 = 5.
Length(cent_basis_blocks) = 5;
#! true
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_CentralizerOfRepresentation
#! @BeginExample
# This is the actual basis for the centralizer.
cent_basis := CentralizerOfRepresentation(rho);;
# All matrices in the span should commute with the image of rho.
ForAll(G, g -> ForAll(cent_basis, M -> Image(rho, g)*M = M*Image(rho,g)));
#! true
#! @EndExample
#! @EndChunk

#! @BeginChunk Example_ClassSumCentralizer
#! @BeginExample
# Now we have a basis for the centralizer, we can sum a conjugacy class
# of G.
class := List(ConjugacyClasses(G)[3]);;
# We can do the computation naively, with no centralizer basis given:
sum1 := ClassSumCentralizer(rho, class, fail);;
# Before summing with th centralizer basis given, we need to
# orthonormalize it. It's already orthogonal, but not normal:
orth_basis := OrthonormalBasis@RepnDecomp(cent_basis);;
IsOrthonormalSet(orth_basis, InnerProduct@RepnDecomp);
#! true
# And with the centralizer given, should be more efficient in certain
# cases (small degree, low multiplicities, but very large group)
sum2 := ClassSumCentralizer(rho, class, orth_basis);;
# Should be the same:
sum1 = sum2;
#! true
#! @EndExample
#! @EndChunk


#! @BeginChunk Example_ClassSumCentralizerNC
#! @BeginExample
# The very same as the above, but with no checks on orthonormality.
sum3 := ClassSumCentralizerNC(rho, class, orth_basis);;
sum1 = sum3;
#! true
#! @EndExample
#! @EndChunk
