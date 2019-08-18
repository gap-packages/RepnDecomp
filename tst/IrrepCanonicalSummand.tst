gap> # there are 3 branches in this function, need to test all
gap> # first when we give an orthonormal basis for the centralizer
gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[2], irreps[2], irreps[3]]);;
gap> diag_info := REPN_ComputeUsingMyMethod(rho);;
gap> rho := diag_info.diagonal_rep;;
gap> cc := ConjugacyClasses(G);;
gap> cent_basis := List(diag_info.centralizer_basis, BlockDiagonalMatrix);; # given in block form, need to convert
gap> cent_basis := List(cent_basis, m -> 1/Sqrt(Trace(m * TransposedMat(m))) * m);; # need to normalize (already ortho)
gap> # coeffs[i] = number of times irreps[i] appears in the decomposition
gap> coeffs := IrrVectorOfRepresentation@RepnDecomp(rho, IrrWithCorrectOrdering@RepnDecomp(G : irreps := irreps));;
gap> # check the blocks are the right size when we give the cent_basis
gap> ForAll([1..Length(irreps)], i -> Length(IrrepCanonicalSummand@RepnDecomp(rho, irreps[i] : centralizer_basis:= cent_basis)) = coeffs[i]*DegreeOfRepresentation(irreps[i]));
true
gap> # also when we use the no cent_basis method, still linear
gap> ForAll([1..Length(irreps)], i -> Length(IrrepCanonicalSummand@RepnDecomp(rho, irreps[i])) = coeffs[i]*DegreeOfRepresentation(irreps[i]));
true
gap> # now a permutation representation
gap> h := RegularActionHomomorphism(G);;
gap> rho := GroupHomomorphismByImages(G, Image(h, G));;
gap> linear_rho := ConvertRhoIfNeeded@RepnDecomp(rho);;
gap> # check block sizes are right and spaces are really G-invariant
gap> coeffs := IrrVectorOfRepresentation@RepnDecomp(linear_rho, IrrWithCorrectOrdering@RepnDecomp(G : irreps := irreps));;
gap> ForAll([1..Length(irreps)], function(i) local basis; basis := IrrepCanonicalSummand@RepnDecomp(linear_rho, irreps[i]); return Length(basis) = coeffs[i]*DegreeOfRepresentation(irreps[i]) and IsGInvariant@RepnDecomp(linear_rho, basis); end );
true