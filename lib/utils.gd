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
#! where $\rho : G \to \mbox{GL}(V)$ is taken from <A>list1</A> and
#! $\tau : H \to \mbox{GL}(W)$ is taken from <A>list2</A>. The result
#! is then a list of representations of $G \times H$.
DeclareGlobalFunction( "TensorProductRepLists" );

#! @Arguments list

#! @Returns Direct sum of the list of representations <A>list</A>
DeclareGlobalFunction( "DirectSumOfRepresentations" );

#! @Arguments rho

#! @Returns Degree of the representation <A>rho</A>. That is,
#! $\mbox{Tr}(\rho(e_G))$, where $e_G$ is the identity of the group
#! $G$ that <A>rho</A> has as domain.
DeclareGlobalFunction( "DegreeOfRepresentation" );

#! @Arguments rho

#! @Returns Linear representation $\rho$ isomorphic to the permutation
#! representation <A>rho</A>.
DeclareGlobalFunction( "PermToLinearRep" );

#! @Arguments S, prod

#! @Returns Whether <A>S</A> is an orthonormal set with respect to the
#! inner product <A>prod</A>.
DeclareGlobalFunction( "IsOrthonormalSet" );

#! @Arguments rho

#! @Returns Whether <A>rho</A> is unitary, i.e. for all $g \in G$,
#! $\rho(g^{-1}) = \rho(g)^*$ (where $^*$ denotes the conjugate
#! transpose).
DeclareGlobalFunction( "IsUnitaryRepresentation" );

#! @Arguments rho

#! @Returns A unitary representation isomorphic to <A>rho</A>
DeclareGlobalFunction( "UnitarizationOfRepresentation" );
