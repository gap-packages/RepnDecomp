# Deterministic tests for small matrix, character, and representation helpers.

gap> START_TEST("core-utils.tst");

gap> Take@RepnDecomp([ 1, 2, 3 ], 0) = [];
true
gap> Take@RepnDecomp([ 1, 2, 3 ], 2) = [ 1, 2 ];
true
gap> Drop@RepnDecomp([ 1, 2, 3 ], 0) = [ 1, 2, 3 ];
true
gap> Drop@RepnDecomp([ 1, 2, 3 ], 3) = [];
true

gap> original := [ 1, 2, 3 ];; copies := Replicate@RepnDecomp(original, 2);;
gap> copies[1][1] := 9;; copies[2] = original and original = [ 1, 2, 3 ];
true
gap> Replicate@RepnDecomp(7, 3) = [ 7, 7, 7 ];
true

gap> BlockDiagonalMatrix@RepnDecomp([]) = [];
true
gap> BlockDiagonalMatrix@RepnDecomp([ IdentityMat(2), IdentityMat(3) ]) = IdentityMat(5);
true
gap> BlockDiagonalMatrix@RepnDecomp([ [[2]], [], [[3,0],[0,4]] ]) = DiagonalMat([2,3,4]);
true

gap> M := DiagonalMat([ 1, 2, 3, 4 ]);;
gap> ExtractBlock@RepnDecomp(M, 1, 1, 2) = DiagonalMat([1,2]);
true
gap> ExtractBlock@RepnDecomp(M, 1, 2, 2) = NullMat(2,2);
true
gap> ExtractBlock@RepnDecomp(M, 2, 2, 2) = DiagonalMat([3,4]);
true

gap> V := VectorSpace(Rationals, [[1,0],[0,1]], [0,0]);;
gap> W := MatrixImage@RepnDecomp([[1,0],[0,0]], V);;
gap> Dimension(W) = 1 and [1,0] in W and not [0,1] in W;
true

gap> IsOrthonormalSet([[1,0],[0,1]], function(x,y) return x*y; end);
true
gap> IsOrthonormalSet([[1,0],[1,0]], function(x,y) return x*y; end);
false

gap> G := SymmetricGroup(3);; irreps := IrreducibleRepresentations(G);;
gap> ForAll(irreps, IsFiniteGroupLinearRepresentation);
true
gap> ForAll(irreps, rho -> DegreeOfRepresentation(rho) = Length(Image(rho, One(G))));
true
gap> chars := List(irreps, CharacterOfRepresentation@RepnDecomp);;
gap> ForAll([1..Length(chars)], i -> ForAll([1..Length(chars)], function(j) if i=j then return InnerProductOfCharacters@RepnDecomp(chars[i],chars[j],G)=1; else return InnerProductOfCharacters@RepnDecomp(chars[i],chars[j],G)=0; fi; end));
true

gap> direct := DirectSumOfRepresentations([irreps[1], irreps[3]]);;
gap> DegreeOfRepresentation(direct) = 3;
true
gap> IrrVectorOfRepresentation@RepnDecomp(direct, chars) = [1,0,1];
true

gap> perm := ActionHomomorphism(G, [1..3]);;
gap> IsFiniteGroupPermutationRepresentation(perm);
true
gap> linear := PermToLinearRep(perm);;
gap> IsFiniteGroupLinearRepresentation(linear) and DegreeOfRepresentation(linear) = 3;
true
gap> ForAll(GeneratorsOfGroup(G), g -> Image(linear,g) = PermutationMat(Image(perm,g),3));
true

gap> STOP_TEST("core-utils.tst", 1);
