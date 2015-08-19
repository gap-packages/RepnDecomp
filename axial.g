#
# Expanded version of the package for computation 
# with axial algebras
# 19/08/2015
#

LoadPackage("grape");

AllTauMaps:=function(G)
 local k,group,i,map,j,m,good,N,o,p,H,graph,M;
 if not IsPermGroup(G) then
  return fail;
 fi;
 k:=LargestMovedPoint(G);
 group:=rec(G:=G);
 group.orbs:=Orbits(G,[1..k]);
 group.reps:=List(group.orbs,o->o[1]);
 group.r:=Length(group.orbs);
 group.maps:=Cartesian(
   List(group.reps,r->
     Filtered(Center(Stabilizer(G,r)),z->z^2=One(G))));
 group.pairs:=Orbits(G,Combinations([1..k],2),OnSets); 
 group.comps:=[];
 group.shs:=[];
 for i in [1..Length(group.maps)] do
  group.comps[i]:=[];
  group.shs[i]:=[];
  map:=[];
  for j in [1..group.r] do
   for m in group.orbs[j] do
    map[m]:=
      group.maps[i][j]^RepresentativeAction(G,group.reps[j],m);
   od;
  od;
  good:=true;
  N:=[];
  o:=[];
  for j in [1..Length(group.pairs)] do
   p:=group.pairs[j][1];
   H:=Subgroup(G,map{p});
   o[j]:=[Set(Orbit(H,p[1])),Set(Orbit(H,p[2]))];
   if o[j][1]=o[j][2] then
    if (Length(o[j][1]) in [3,5]) 
         and Size(H)=2*Length(o[j][1]) then
     N[j]:=Length(o[j][1]);
    else
     good:=false;
     break;
    fi;
   else
    if Length(o[j][1])=Length(o[j][2]) 
         and (Length(o[j][2]) in [1,2,3]) 
         and (Size(H) in 
           [1,2*Length(o[j][1]),4*Length(o[j][1])]) then
     N[j]:=2*Length(o[j][1]);
    else
     good:=false;
     break;
    fi;
   fi;
  od;
  if good then
   graph:=NullGraph(Group(()),Length(group.pairs));
   for j in [1..Length(group.pairs)] do
    if N[j]=4 then
     for m in [1..Length(group.pairs)] do
      if N[m]=2 then
       if o[j][1] in group.pairs[m] then
        break;
       fi;
      fi;
     od;
     if not(m in Adjacency(graph,j)) then
      AddEdgeOrbit(graph,[j,m]);
      AddEdgeOrbit(graph,[m,j]);
     fi;
     for m in [1..Length(group.pairs)] do
      if N[m]=2 then
       if o[j][2] in group.pairs[m] then
        break;
       fi;
      fi;
     od;
     if not(m in Adjacency(graph,j)) then
      AddEdgeOrbit(graph,[j,m]);
      AddEdgeOrbit(graph,[m,j]);
     fi;
    fi;
    if N[j]=6 then
     for m in [1..Length(group.pairs)] do
      if N[m]=3 then
       if o[j][1]{[1,2]} in group.pairs[m] then
        break;
       fi;
      fi;
     od;
     if not(m in Adjacency(graph,j)) then
      AddEdgeOrbit(graph,[j,m]);
      AddEdgeOrbit(graph,[m,j]);
     fi;
     for m in [1..Length(group.pairs)] do
      if N[m]=3 then
       if o[j][2]{[1,2]} in group.pairs[m] then
        break;
       fi;
      fi;
     od;
     if not(m in Adjacency(graph,j)) then
      AddEdgeOrbit(graph,[j,m]);
      AddEdgeOrbit(graph,[m,j]);
     fi;
     p:=Concatenation([o[j][1][1]],
          Filtered(o[j][2],q->q^map[o[j][1][1]]=q));
     p:=Set(p);
     for m in [1..Length(group.pairs)] do
      if N[m]=2 then
       if p in group.pairs[m] then
        break;
       fi;
      fi;
     od;
     if not(m in Adjacency(graph,j)) then
      AddEdgeOrbit(graph,[j,m]);
      AddEdgeOrbit(graph,[m,j]);
     fi;
    fi; 
   od;
   group.comps[i]:=ConnectedComponents(graph);
   M:=List(group.comps[i],c->Maximum(List(c,q->N[q])));
   group.shs[i]:=Cartesian(List(M,
                   function(m) if m=2 then 
                                return ["2A","2B"];
                               elif m=3 then
                                return ["3A","3C"];
                               elif m=4 then
                                return ["4A","4B"];
                               elif m=5 then
                                return ["5A"];
                               else
                                return ["6A"];
                               fi; end));
  fi;
 od;
 return group;
end;
