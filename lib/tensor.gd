#! @Chapter Representations of tensor products

#! @Section Space-efficient representation of tensors of matrices

#! Suppose we have representations of $G$, $\rho$ and $\tau$, with
#! degree $n$ and $m$. If we would like to construct the tensor
#! product representation of $G$, $\rho \otimes \tau$, the usual way
#! to do it would be to take the Kronecker product of the
#! matrices. This means we now have to store very large $nm \times nm$
#! matrices for each generator of $G$.

#! This can be avoided by storing the tensor of matrices as pairs,
#! essentially storing $A \otimes B$ as a pair $(A,B)$ and
#! implementing group operations on top of these, along with some
#! representation-theoretic functions.

#! Even though matrices do form a ring, it is only possible to
#! guarantee an economical representation for pure tensors,
#! i.e. matrices of the form $A \otimes B$. These are closed under
#! group operations, so it is natural to define a group structure.
DeclareCategory( "IsTensorProductOfMatricesObj", IsMultiplicativeElementWithInverse );

#! Position $i$ in this representation stores the matrix $A_i$ in the
#! tensor product $A_1 \otimes A_2$.
DeclareRepresentation( "IsTensorProductPairRep", IsPositionalObjectRep, [1, 2] );

DeclareGlobalFunction( "TensorProductOfMatricesObj" );

DeclareGlobalFunction( "TensorProductOfMatricesFamily" );
