# call agens(n), e.g. for n=7 do
# gap> blah:=agens(7,"fileout");;
# to get \overline{L}_k, etc.
# namely, blah.gens will be these matrices,
# blah.v will be their row sums,
# blah.prime will be indices of the images of generators
# under *.
# n is the 1st argument; 2nd argument is the filename to write the data out 
Read("classmatr.g");

cvec:=function(n,l)
local v,i,x;
v:=List([1..n],x->0);
for i in l do
v[i]:=1;
od;
return v;
end;

p2p:=function(t)
local p,n,z,k;
z:=Concatenation([1],t);
n:=Length(z);
p:=List([1..n]);
for k in [1..n-1] do
  p[z[k]]:=z[k+1];
od;
p[z[n]]:=1;
return PermList(p);
end;

# returns the group 2xS_n in its action on the n-cycles by
# conjugation, and the corresponding n-cycles. 
# The first two generators generate the point
# stabilizer in this action, and the 4th (redundant) is the pairing
# g<->g^-1.
#
grp:=function(n)
local g,cn,x,inv,inv1,gg;
g:=SymmetricGroup(n);
if Size(Group(g.1))<>n or g.2<>(1,2) then 
        Error("g is generated in an unexpected way\n");
return;
fi;
inv:=PermList(Concatenation([1],List([1..n-1],x->n+1-x)));
cn:=List(PermutationsList([2..n]),p2p);
inv1:=PermList(List(cn,x->Position(cn,x^-1)));
gg:=List([inv,g.1,g.2],x->Permutation(x,cn));
return rec(g:=GroupWithGenerators([gg[1]*inv1,gg[2],gg[3],inv1]), c:=cn);
end;

# algebra generators (non-normalized)
# n is the 1st argument; 
# then also adjacency matrices for each graph are returned
#
agens:=function(n,fname)
local printarr, k, givemats,adists,c,graphs,y,x,d,L,g,h,a,ogcm,G,z,
	T,diam,p,v,oreps,gc,hdists,horbs;
gc:=grp(n);
g:=gc.g;
c:=gc.c;
h:=Subgroup(g,[g.1,g.2]);
T:=TwoOrbitNumbers(g,[h]);
StructureConstantsTwoOrbitsAlgebra(T);
G:=NullGraph(g);
AddEdgeOrbit(G,[1,1^g.3],h);
horbs:=Orbits(h,[1..G.order]);
adists:=List(horbs,x->Distance(G,1,x[1]^g.4));

v:=List([1..Length(horbs)],x->StructureConstantTwoOrbitsAlgebra(T,1,x,x^T.pairing));

printarr:=function(a)
local k, s;
s:=false;
for k in a do 
  if s then AppendTo(fname, " "); fi;
  AppendTo(fname,  k); 
  s:=true;
od;
AppendTo(fname,"\n");
return 1;
end;

PrintTo(fname); # clean the output file
printarr(adists);
printarr(Concatenation(ListPerm(T.pairing),
  	[LargestMovedPoint(T.pairing)+1..TwoOrbitAlgebraDimension(T)]));
printarr(v);
for k in [1..T.dimension] do
 for z in [1..T.dimension] do
   for x in T.P[k][z] do
	AppendTo(fname, z, " ", k, " ", x[1], " ", x[2], "\n");
   od;
 od;
od;

return 1;
end;
