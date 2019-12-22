gap> G := SmallGroup(24, 3);;
gap> chars := Irr(G);;
gap> irreps := IrreducibleRepresentationsDixon(G);;
gap> # check our conversion is correct relying on the fact that chars and irreps are ordered the same
gap> ForAll([1..Length(irreps)], i -> CharacterOfRepresentation@RepnDecomp(irreps[i]) = List(chars[i]));
true