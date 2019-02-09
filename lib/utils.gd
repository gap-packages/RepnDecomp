#! @Chapter Utility functions

#! @Section Miscellaneous

#! @Arguments blocks

#! @Returns Matrix given by putting the given matrix <A>blocks</A> on
#! the diagonal
DeclareGlobalFunction( "BlockDiagonalMatrix" );

#! @Arguments hom, func

#! @Returns Homomorphism g given by g(x) = func(hom(x)).

#! @Description This is mainly for convenience, since it handles all
#! GAP accounting issues regarding the range, ByImages vs ByFunction,
#! etc.
DeclareGlobalFunction( "ComposeHomFunction" );

#! @Arguments elem, n

#! @Returns List of <A>n</A> copies of <A>elem</A>
DeclareGlobalFunction( "Replicate" );

#! @Section Representation theoretic functions

#! @Arguments list1, list2

#! @Returns All possible tensor products given by $\rho \otimes \tau$
#! where $\rho$ is taken from <A>list1</A> and $\tau$ is taken from
#! <A>list2</A>.
DeclareGlobalFunction( "TensorProductRepLists" );

#! @Arguments list

#! @Returns Direct sum of the list of representations <A>list</A>
DeclareGlobalFunction( "DirectSumRepList" );

#! @Arguments rho

#! @Returns Degree of the representation <A>rho</A>. That is,
#! $\mbox{Tr}(\rho(e_G))$, where $e_G$ is the identity of the group
#! $G$ that <A>rho</A> has as domain.
DeclareGlobalFunction( "DegreeOfRepresentation" );

#! @Arguments rho

#! @Returns Linear representation $\rho$ isomorphic to the permutation
#! representation <A>rho</A>.
DeclareGlobalFunction( "PermToLinearRep" );
