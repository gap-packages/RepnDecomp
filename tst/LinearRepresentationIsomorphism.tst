gap> G := SymmetricGroup(4);;
gap> irreps := IrreducibleRepresentations(G);;
gap> # rho, tau, tau2 are isomorphic. We check if we can find the
gap> # isomorphisms between them
gap> tau := DirectSumOfRepresentations([irreps[3], irreps[1], irreps[3]]);;
gap> rho := DirectSumOfRepresentations([irreps[1], irreps[3], irreps[3]]);;
gap> B := RandomInvertibleMat(5);;
gap> tau2 := ComposeHomFunction(tau, x -> B^-1 * x * B);;
gap> M := LinearRepresentationIsomorphism(rho, tau);;
gap> IsRepresentationIsomorphism@RepnDecomp(rho, tau, M);
true
gap> M := LinearRepresentationIsomorphism(rho, tau2);;
gap> IsRepresentationIsomorphism@RepnDecomp(rho, tau2, M);
true
gap> M := LinearRepresentationIsomorphism(tau, tau2);;
gap> IsRepresentationIsomorphism@RepnDecomp(tau, tau2, M);
true