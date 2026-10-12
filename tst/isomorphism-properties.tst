# Deterministic inputs for all explicit-intertwiner implementations.

gap> START_TEST("isomorphism-properties.tst");

gap> G := SymmetricGroup(4);; irreps := IrreducibleRepresentations(G);;
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[3], irreps[3]]);;
gap> B := REPN_TEST_UpperUnitriangularMat(5);;
gap> tau := REPN_TEST_ConjugateRepresentation(rho, B);;
gap> AreRepsIsomorphic(rho, tau);
true
gap> IsLinearRepresentationIsomorphism(B^-1, rho, tau);
true
gap> IsLinearRepresentationIsomorphism(NullMat(5,5), rho, tau);
false

gap> iso := LinearRepresentationIsomorphism(rho, tau);;
gap> IsLinearRepresentationIsomorphism(iso, rho, tau);
true
gap> isoK := LinearRepresentationIsomorphism(rho, tau : use_kronecker);;
gap> IsLinearRepresentationIsomorphism(isoK, rho, tau);
true
gap> isoOrbit := LinearRepresentationIsomorphism(rho, tau : use_orbit_sum);;
gap> IsLinearRepresentationIsomorphism(isoOrbit, rho, tau);
true
gap> isoSlow := LinearRepresentationIsomorphismSlow(rho, tau);;
gap> IsLinearRepresentationIsomorphism(isoSlow, rho, tau);
true

gap> AreRepsIsomorphic(irreps[1], irreps[2]);
false
gap> LinearRepresentationIsomorphism(irreps[1], irreps[2]);
fail
gap> AreRepsIsomorphic(rho, irreps[1]);
false

gap> H := CyclicGroup(2);; hrep := IrreducibleRepresentations(H)[1];;
gap> AreRepsIsomorphic(irreps[1], hrep);
false

gap> STOP_TEST("isomorphism-properties.tst", 1);
