LoadPackage("grape");

TwoOrbitNumbers:=function(arg)
# G=arg[1] - the group, 
# H=arg[2] (optional) list of point stabilisers, one for each orbit
local os, x, ctr, it, deg,
   gpairs, ttorbs, topair,
   len, onum, gens, w, onums, n, G, H, O, O2, d, nG, I, i, J, j, j0, sv;

G:=arg[1];

if not IsPermGroup(G) then 
   Error("usage: TwoOrbitNumbers( <PermGroup> [, stabilisers] )");
fi;
n:= LargestMovedPoint(G);
O:=GRAPE_OrbitNumbers(G,n);
# orbs:=List(O.representatives,w->Orbit(G,w));

d:=Length(O.representatives); # number of orbits of G
if IsBound(arg[2]) then
  if IsList(arg[2]) then H:=arg[2];
  else H:=[arg[2]];
  fi;
  if Length(H)<>d then 
     Error("wrong number of point stabilisers supplied");
  fi;
  for i in [1..d] do
  if ForAny(GeneratorsOfGroup(H[i]),x->O.representatives[i]^x<>O.representatives[i]) then
      Error("H[",i,"] does not fix the point ",O.representatives[i]);
  fi;
  od;
else 
  H:=List(O.representatives, i->Stabilizer(G,i));
fi;

onums:=[];
ttorbs:=0;
for I in [1..d] do 
# Add(onums, List([1..d], rec(orbnums:=GRAPE_OrbitNumbers(H[I],n), paired:=[])));
 w:=GRAPE_OrbitNumbers(H[I],n);
 # need to localize orbit representatives per G-orbit.
 ctr:=List([1..d],i->0);
 len:=Length(w.representatives);
 onum:=List([1..len],i->0);
 deg:=List([1..len],i->0);
# computing subdegrees
 for i in w.orbitNumbers do
   deg[i]:=deg[i]+1;
 od;
 os:=List([1..d],i->[]);
 for i in [1..len] do
    x:=O.orbitNumbers[w.representatives[i]]; # index of the corr. G-orbit
    Add(os[x],w.representatives[i]);         # index of the corr. G-orbit
    ctr[x]:=ctr[x]+1; # the per-block number of this H[I] - orbit
    onum[i]:=ctr[x];
 od;
 # internally, 2-orbits are indexed by pairs (G-orbit, H[I]-orbit)
 Add(onums, rec(globalTwoOrbitNumbers:=List(w.orbitNumbers, i->i+ttorbs),
 		stabiliser:=w,
                representatives:=os,
		perBlockOrbitNumber:=onum,
		numberTwoOrbits:=Length(Flat(os)),
		paired:=List([1..d],x->[]),
		subdegree:=deg));
 ttorbs:=ttorbs+len;
od;

nG:=NullGraph(G);
sv:=nG.schreierVector;
gens:=GeneratorsOfGroup(G);

# establishing pairings
gpairs:=List([1..ttorbs],i->0);

for I in [1..d] do
# deal with the I-th block row, and the global row O.representatives[I]
 i := O.representatives[I]; # global row number (fixed by H[I])
 for J in [1..d] do 
# deal with the J-th block column, and the global columns 
#    onums[I].representatives[J]
# onums[I].paired[J]:=...; #must be paired to something in onums[J]
   for j0 in [1..Length(onums[I].representatives[J])] do
     j:=onums[I].representatives[J][j0]; # global column number
     topair:=onums[I].globalTwoOrbitNumbers[j];
# need to canonise (j,i), i.e. find g in G that maps j to O.representatives[J].
#  (using the Schreier vector sv in the G's nullgraph for this)
# then g(i) should be located in the orbits of the stabilizer of 
# O.representatives[J], i.e. in onums[J][I].stabiliser.orbitNumbers.
# So onums[I].paired[J][j0] must get the corresponding entry assigned.
     w:=sv[j];
     it:=i;
     while w>0 do
        j:=j/gens[w];
	it:=it/gens[w];
	w:=sv[j];
     od;

     gpairs[topair]:=onums[J].globalTwoOrbitNumbers[it];

     onums[I].paired[J][j0]:=
     	onums[J].perBlockOrbitNumber[onums[J].stabiliser.orbitNumbers[it]];
#Error("2");
     if J>I then # ??? or onums[J].stabiliser.orbitNumbers[i]<>j0 then 
       onums[J].paired[I][
         onums[J].perBlockOrbitNumber[onums[J].stabiliser.orbitNumbers[it]]]:=
	     onums[J].stabiliser.orbitNumbers[onums[J].representatives[I][j0]];
     fi;
   od;
  od;
od;
return rec(isTwoOrbitNumbers:=true, group:=G,
   n:=n, orbitNumbers:=O.orbitNumbers, P:=[],
   orbitNumberRepresentatives:=O.representatives, 
   orbitLengths:=List(O.representatives, w->OrbitLength(G,w)),
   twoOrbitNumbers:=onums, schreierVector:=sv,
   dimension:=ttorbs,
   pairing:=PermList(gpairs));
end;

###########################################################################
# standard function to test whether this is the consistent record
###########################################################################
IsTwoOrbitNumbers := function(T)
  return IsRecord(T) and IsBound(T.isTwoOrbitNumbers) and T.isTwoOrbitNumbers;
end;

###########################################################################
# underlying group 
###########################################################################
TwoOrbitsAlgebraGroup := function(T)
  if not IsTwoOrbitNumbers(T) then 
    Error("TwoOrbitsAlgebraGroup: this is not a coherent configuration! "); 
  fi;
  return T.group;
end;

###########################################################################
# size of the domain operated on by the underlying group 
###########################################################################
TwoOrbitsAlgebraDegree := function(T)
  return T.n;
end;

###########################################################################
# orbit lengths of the underlying group 
###########################################################################
TwoOrbitsAlgebraOrbitLengths := function(T)
  return T.orbitLengths;
end;

NumberTwoOrbitsInBlock := function(TN,I)
return TN.twoOrbitNumbers[I].numberTwoOrbits;
end;

StabiliserOrbitNumberInBlock := function(TN, onum, blk)
local b;
for b in [1..blk-1] do
  onum:=onum-Length(TN.twoOrbitNumbers[b].stabiliser.representatives);
od;
return onum;
end;

# convert the rec(orbitNumber,twoOrbitNumber) into the global number
GlobalTwoOrbitNumber := function(TN,lnum)
local j, s, I;
I:=lnum.orbitNumber;
s:=lnum.twoOrbitNumber;
for j in [1..I-1] do
 s:=s+NumberTwoOrbitsInBlock(TN,j);
od;
return s;
end;

# and the reverse conversion
LocalTwoOrbitNumber := function(TN, i)
local I, lnum;
I:=1;
while i>NumberTwoOrbitsInBlock(TN,I) do 
 i:=i-NumberTwoOrbitsInBlock(TN,I);
 I:=I+1;
od;
lnum:=rec(orbitNumber:=I, twoOrbitNumber:=i);
return lnum;
end;

########################################################################
# for an arc (a,b), find the (global) number of the two-orbit it lies in
########################################################################
TwoOrbitNumber := function(TN, a, b)
local sv, w, gens, I, lnum, globalNumber;
if not IsTwoOrbitNumbers(TN) or  a>TwoOrbitsAlgebraDegree(TN) or b>TwoOrbitsAlgebraDegree(TN) or 
       not IsPosInt(a) or not IsPosInt(b) then
  Error("usage: TwoOrbitNumber(TwoOrbitNumbers struct, int, int)");
fi;
sv:=TN.schreierVector;
gens:=GeneratorsOfGroup(TwoOrbitsAlgebraGroup(TN));
w:=sv[a];
while w>0 do
 a:=a/gens[w];
 b:=b/gens[w];
 w:=sv[a];
od;
I:=TN.orbitNumbers[a];
return TN.twoOrbitNumbers[I].globalTwoOrbitNumbers[b];
end;

NumberOfBlocks := function(TN)
if not IsTwoOrbitNumbers(TN) then
  Error("usage: NumberOfBlocks(TwoOrbitNumbers struct)");
fi;
return Length(TN.orbitNumberRepresentatives);
end;

########################################################################
# dimension of the algebra
########################################################################
TwoOrbitAlgebraDimension:= function(TN)
local x;
if not IsTwoOrbitNumbers(TN) then
  Error("usage: TwoOrbitAlgebraDimension(TwoOrbitNumbers struct)");
fi;
return TN.dimension;
end;

# compute the G-0rbit number
GlobalOrbitNumberOfTwoOrbit:=function(TN,I,t)
local j;
j:=TN.twoOrbitNumbers[I].stabiliser.representatives[t];
return TN.orbitNumbers[j];
end;


#########################################################################
# list a pair for each 2-orbit
#########################################################################
TwoOrbitRepresentatives := function(T)
local i,x;
return Concatenation(List([1..Length(T.orbitNumberRepresentatives)], i->
	List(Flat(T.twoOrbitNumbers[i].representatives),x->
		[T.orbitNumberRepresentatives[i],x])));
end;

#########################################################################
# a naive implementation ignoring block structures
# 6.12.08 - added sparsity handling
#########################################################################
MakeStructureConstantsData:=function(T, ci, cj, I, J, k)
local t,p,i, M;
M:=List([1..T.dimension],i->Set([]));
for i in [1..TwoOrbitsAlgebraDegree(T)] do 
  p:=PositionSet(List(M[ci[i]],t->t[1]),cj[i]);
  if p=fail then 
    AddSet(M[ci[i]],[cj[i],1]);
  else 
    M[ci[i]][p][2] := M[ci[i]][p][2]+1;
  fi;
od;
#M:=NullMat(T.dimension, T.dimension);
#for i in [1..TwoOrbitsAlgebraDegree(T)] do 
#  M[ci[i]][cj[i]] := M[ci[i]][cj[i]]+1;
#od;
T.P[k]:=M;
return 1;
end;


#########################################################################
# for a given k, compute p^k_{ij} for all i, j
#########################################################################
StructureConstantsTwoOrbitAlgebra := function(TN, k)
local word, w, sv, gens, x, iii, ci, I, twoonum, Ir,
gamma, cj, J, M, i,j; 

Ir:=LocalTwoOrbitNumber(TN, k);
I:=Ir.orbitNumber;
twoonum:=Ir.twoOrbitNumber;

# locate the arc (i,j) from the (I,J)-block
j:=TN.twoOrbitNumbers[I].stabiliser.representatives[twoonum];

# the k-th 2-orbit is located in (I,J)-block
J:=GlobalOrbitNumberOfTwoOrbit(TN,I,twoonum);
i:=TN.orbitNumberRepresentatives[I];
# now we have to compute the "scalar product" of i-th row and j-th column
# of the "big matrix"
# and store it in M. More precisely, M[a][b] will get the number of points
# p such that the arc (i,p) has colour a and the arc (p,j) - colour b.
#
# the i-th row is already there, it is 
# ci=TN.twoOrbitNumbers[I].globalTwoOrbitNumbers;
#
# the j-th column need to be computed from 
# TN.twoOrbitNumbers[J].globalTwoOrbitNumbers by applying the suitable
# permutation and then the pairing.

ci:=TN.twoOrbitNumbers[I].globalTwoOrbitNumbers;

sv:=TN.schreierVector;
gens:=GeneratorsOfGroup(TwoOrbitsAlgebraGroup(TN));
w:=sv[j];

# need to compute the word in generators,

word:=[];
while w>0 do
 j:=j/gens[w];
 Add(word,w);
 w:=sv[j];
od;

# and apply it to cj
cj:=List(TN.twoOrbitNumbers[J].globalTwoOrbitNumbers, w->
	w^TN.pairing);
for w in Reversed(word) do
 cj:=Permuted(cj,gens[w]);
od;
MakeStructureConstantsData(TN, ci, cj, I, J, k);
return 1;
end;

StructureConstantsTwoOrbitsAlgebra := function(TN)
local k;
if TN.P=[] then # compute them
	List([1..TwoOrbitAlgebraDimension(TN)],k-> 
		StructureConstantsTwoOrbitAlgebra(TN,k));
	return 1;
fi;
return 2; # they were already computed
end;

##########################################################################
# computing p^a_{bc}
##########################################################################
StructureConstantTwoOrbitsAlgebra := function(TN,a,b,c)
local p,x;
#return TN.P[a][b][c]; # naive, non-sparse implementation
p:=PositionSet(List(TN.P[a][b],x->x[1]),c);
if p=fail then return 0;
else return TN.P[a][b][p][2];
fi;
end;

###########################################################################
# computing the isomorphism to the regular representation, i.e.
# A_b -> (L_b)_{ac}:=p^a_{bc}, for a,b,c=1,...,algebra dimension
###########################################################################
RegularRepresentationMatrix := function(TN,b)
local L, d, a, c;
d := TwoOrbitAlgebraDimension(TN);
L := NullMat(d,d);
for a in [1..d] do
  for c in [1..d] do
     L[a][c]:=StructureConstantTwoOrbitsAlgebra(TN,a,b,c);
  od;
od;
return L;
end;

SparseRegularRepresentationMatrix := function(TN,b)
local v,a,c,d, L;
d := TwoOrbitAlgebraDimension(TN);
L:=List([1..d],i->Set([]));
for a in [1..d] do
  for c in [1..d] do
     v:=StructureConstantTwoOrbitsAlgebra(TN,a,b,c);
     if v<>0 then
       AddSet(L[a],[c,v]);
     fi;
  od;
od;
return L;
end;


SparseMatrixTranspose := function(A)
local p, d, i,j,L;
d:=Length(A);
L:=List([1..d],i->Set([]));
for i in [1..d] do
  for p in A[i] do
     AddSet(L[p[1]],[i,p[2]]);
  od;
od;
return L;
end;


