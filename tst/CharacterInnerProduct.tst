gap> # some random group with lots of characters
gap> G := SmallGroup(56, 5);;
gap> # check all characters have norm 1
gap> ForAll(Irr(G), char -> InnerProductOfCharacters@RepnDecomp(char, char, G) = 1);
true
gap> # check all distinct characters are orthogonal
gap> ForAll(Irr(G), chi1 -> ForAll(Irr(G), chi2 -> chi1 = chi2 or InnerProductOfCharacters@RepnDecomp(chi1, chi2, G) = 0));
true
