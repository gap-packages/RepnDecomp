#! @Chapter Miscellaneous useful functions

#! @Section Predicates for representations

#! @Arguments rho

#! @Returns true or false

#! @Description Tells you if <A>rho</A> is a linear representation of
#! a finite group. The algorithms implemented in this package work on
#! these homomorphisms only.
DeclareAttribute( "IsFiniteGroupLinearRepresentation", IsGroupHomomorphism );


#! @Arguments rho

#! @Returns true or false

#! @Description Tells you if <A>rho</A> is a homomorphism from a
#! finite group to a permutation group.
DeclareAttribute( "IsFiniteGroupPermutationRepresentation", IsGroupHomomorphism );
