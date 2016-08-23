Read("projs.g");
s4_0:=SymmetricGroup(4);
o:=
     [ [ [ 1 ], [ 2 ], [ 3 ] ], [ [ 1 ], [ 2 ], [ 4 ] ], [ [ 1 ], [ 3 ], [ 4 ] ],
     [ [ 2 ], [ 3 ], [ 4 ] ], [ [ 1 ], [ 2, 2 ] ],
      [ [ 1 ], [ 3, 3 ] ], [ [ 1 ], [ 4, 4 ] ], [ [ 1, 1 ], [ 2 ] ], [ [ 1, 1 ], [ 3 ] ],
      [ [ 1, 1 ], [ 4 ] ], [ [ 2 ], [ 3, 3 ] ], [ [ 2 ], [ 4, 4 ] ],
      [ [ 2, 2 ], [ 3 ] ], [ [ 2, 2 ], [ 4 ] ], [ [ 3 ], [ 4, 4 ] ], [ [ 3, 3 ], [ 4 ] ] ];
s4:=Action(s4_0,o,OnSetsTuples);
s4p:=Group(PermutationMat(s4.1,16),PermutationMat(s4.2,16)); # the permutation module we decompose
irrs:=IrreducibleRepresentationsDixon(s4);

# for non-permutation groups, we need such a conversion (not needed in this example)
# irrs:=List(irrs,rr->GroupHomomorphismByImages(s4,Group(Image(rr,s4.1),Image(rr,s4.2)),
#        [s4.1,s4.2],[Image(rr,s4.1),Image(rr,s4.2)]));
rep:=GroupHomomorphismByImages(s4,s4p,[s4.1,s4.2],[s4p.1,s4p.2]);

M:=conjmat(irrs,rep); # M[2] holds the sparsity pattern

cmats:=List([s4p.1,s4p.2],x->M[1]*x*M[1]^-1); # block-decomposed matrices

# they are better used packed, as follows ( so you don't need to call conjmat directly)
cmats_packed:=conjrep(irrs,rep);
# the latter returns list of generators in block-diagonal form; GAP does not allow
# generic sparse matrices, so we just keep each generator as a list of diagonal blocks
