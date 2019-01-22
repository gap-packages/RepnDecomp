gap> G := SymmetricGroup(4);;
gap> T := CohCfgFromPermGroup(G);;
gap> reps := CCTransversal(T);;
gap> n := LargestMovedPoint(G);;
gap> all_one := Replicate@RepnDecomp(Replicate@RepnDecomp(1, n), n);;
gap> # make sure orbitals include all elements
gap> Sum(reps, rep -> OrbitalMatrix@RepnDecomp(G, rep)) = all_one;
true