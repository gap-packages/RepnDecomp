#! @Arguments rho

#! @Returns a unitary representation isomorphic to <A>rho</A>
DeclareGlobalFunction( "UnitaryRepresentation" );

#! @Arguments rho

#! @Returns Whether <A>rho</A> is unitary, i.e. for all $g \in G$,
#! $\rho(g^{-1}) = \rho(g)^*$ (where $^*$ denotes the conjugate
#! transpose).
DeclareGlobalFunction( "IsUnitaryRepresentation" );

#! @Arguments A

#! @Returns a pair $[L, D]$ such that $A = L\mbox{diag}(D)L^*$. $D$ is
#! the $1 \times n$ vector which gives the diagonal matrix
#! $\mbox{diag}(D)$ (where <A>A</A> is an $n \times n$ matrix).
DeclareGlobalFunction( "LDLDecomposition" );
