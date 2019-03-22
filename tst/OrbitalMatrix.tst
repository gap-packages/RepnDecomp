gap> G := SymmetricGroup(4);;
gap> T := CohCfgFromPermGroup@RepnDecomp(G);;
gap> reps := CCTransversal@RepnDecomp(T);;
gap> n := LargestMovedPoint(G);;
gap> all_one := Replicate@RepnDecomp(Replicate@RepnDecomp(1, n), n);;
gap> orbitals := List(reps, rep -> OrbitalMatrix@RepnDecomp(G, rep));;
gap> # make sure orbitals include all elements
gap> Sum(orbitals) = all_one;
true
gap> # check they are all zero-one matrices
gap> ForAll(orbitals, orbital -> ForAll(orbital, row -> ForAll(row, entry -> entry = 0 or entry = 1)));
true