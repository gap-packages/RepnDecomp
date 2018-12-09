#! @Chapter Useful predicates

#! @Section Types of group representations

#! @Arguments rho

#! @Returns true or false

#! @Description Tells you if <A>rho</A> is a linear representation of
#! a finite group. This is important since Serre's algorithms only
#! work on these.
DeclareAttribute( "IsFiniteGroupLinearRepresentation", IsGroupHomomorphism );


#! @Arguments rho

#! @Returns true or false

#! @Description Tells you if <A>rho</A> is a homomorphism from finite
#! group to a permutation group. Such homomorphisms occur often in
#! applications.
DeclareAttribute( "IsFiniteGroupPermutationRepresentation", IsGroupHomomorphism );
