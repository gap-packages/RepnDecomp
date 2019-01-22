gap> # some random group
gap> G := DirectProduct(SymmetricGroup(3), SymmetricGroup(3));;
gap> irreps := IrreducibleRepresentations(G);;
gap> # rho is irreps2 oplus irreps2 oplus irreps2 (rho has degree 3)
gap> rho := DirectSumRepList([irreps[2], irreps[2], irreps[2]]);;
gap> # so the canonical summand is just the whole space
gap> summand := Cyclotomics^3;;
gap> # and each axis is an irreducible G-invariant space
gap> DecomposeCanonicalSummand@RepnDecomp(rho, irreps[2], summand); # should get 3 bits, 1 for each irrep
[ rec( basis := [ [ 1, 0, 0 ] ], space := <vector space over Cyclotomics, with
        1 generators> ),
  rec( basis := [ [ 0, 1, 0 ] ], space := <vector space over Cyclotomics, with
        1 generators> ),
  rec( basis := [ [ 0, 0, 1 ] ], space := <vector space over Cyclotomics, with
        1 generators> ) ]