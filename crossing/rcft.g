BindGlobal("DoDiagonalizeMat2",function(arg)
local R,M,transform,divide,swaprow, swapcol, addcol, addrow, multcol, multrow, l, n, start, d,
      typ, ed, posi,posj, a, b, qr, c, i,j,left,right,cleanout, alldivide,basmat,origtran;

  R:=arg[1];
  M:=arg[2];
  transform:=arg[3];
  divide:=arg[4];

  l:=Length(M);
  n:=Length(M[1]);

  basmat:=fail;
  if transform then
    left:=IdentityMat(l,R);
    right:=IdentityMat(n,R);
    if Length(arg)>4 then
      origtran:=arg[5];
      basmat:=IdentityMat(l,R); # for RCF -- transpose of P' in D&F, sec. 12.2
    fi;
  fi;

  swaprow:=function(a,b)
  local r;
    r:=M[a];
    M[a]:=M[b];
    M[b]:=r;
    if transform then
      r:=left[a];
      left[a]:=left[b];
      left[b]:=r;
      if basmat<>fail then
	r:=basmat[a];
	basmat[a]:=basmat[b];
	basmat[b]:=r;
      fi;
    fi;
  end;

  swapcol:=function(a,b)
  local c;
    c:=M{[1..l]}[a];
    M{[1..l]}[a]:=M{[1..l]}[b];
    M{[1..l]}[b]:=c;
    if transform then
      c:=right{[1..n]}[a];
      right{[1..n]}[a]:=right{[1..n]}[b];
      right{[1..n]}[b]:=c;
    fi;
  end;

  addcol:=function(a,b,m)
  local i;
    for i in [1..l] do
      M[i][a]:=M[i][a]+m*M[i][b];
    od;
    if transform then
      for i in [1..n] do
        right[i][a]:=right[i][a]+m*right[i][b];
      od;
    fi;
  end;

  addrow:=function(a,b,m)
    AddCoeffs(M[a],M[b],m);
    if transform then
      AddCoeffs(left[a],left[b],m);
      if basmat<>fail then
        basmat[b]:=basmat[b]-basmat[a]*Value(m,origtran);
      fi;
    fi;
  end;

  multcol:=function(a,m)
  local i;
    for i in [1..l] do
      M[i][a]:=M[i][a]*m;
    od;
    if transform then
      for i in [1..n] do
        right[i][a]:=right[i][a]*m;
      od;
    fi;
  end;

  multrow:=function(a,m)
    MultRowVector(M[a],m);
    if transform then
      MultRowVector(left[a],m);
      if basmat<>fail then
	MultRowVector(basmat[a],1/m);
      fi;
    fi;
  end;

  # clean out row and column
  cleanout:=function()
  local a,i,b,c,qr;
    repeat
      # now do the GCD calculations only in row/column
      for i in [start+1..n] do
        a:=i;
        b:=start;
        if not IsZero(M[start][b]) then
          repeat
            qr:=QuotientRemainder(R,M[start][a],M[start][b]);
            addcol(a,b,-qr[1]);
            c:=a;a:=b;b:=c;
          until IsZero(qr[2]);
          if b=start then
            swapcol(start,i);
          fi;
        fi;

        # normalize
        qr:=StandardAssociateUnit(R,M[start][start]);
        multcol(start,qr);

      od;

      for i in [start+1..l] do
        a:=i;
        b:=start;
        if not IsZero(M[b][start]) then
          repeat
            qr:=QuotientRemainder(R,M[a][start],M[b][start]);
            addrow(a,b,-qr[1]);
            c:=a;a:=b;b:=c;
          until IsZero(qr[2]);
          if b=start then
            swaprow(start,i);
          fi;
        fi;

        qr:=StandardAssociateUnit(R,M[start][start]);
        multrow(start,qr);

      od;
    until ForAll([start+1..n],i->IsZero(M[start][i]));
  end;

  start:=1;
  while start<=Length(M) and start<=n do

    # find element of lowest degree and move it into pivot
    # hope is this will reduce the total number of iterations by making
    # it small in the first place
    d:=infinity;

    for i in [start..l] do
      for j in [start..n] do
        if not IsZero(M[i][j]) then
          ed:=EuclideanDegree(R,M[i][j]);
          if ed<d then
            d:=ed;
            posi:=i;
            posj:=j;
          fi;
        fi;
      od;
    od;

    if d<>infinity then # there is at least one nonzero entry

      if posi<>start then
        swaprow(start,posi);
      fi;
      if posj<>start then
        swapcol(start,posj);
      fi;
      cleanout();

      if divide then
        repeat
          alldivide:=true;
          #make sure the pivot also divides the rest
          for i in [start+1..l] do
            for j in [start+1..n] do
              if Quotient(M[i][j],M[start][start])=fail then
                alldivide:=false;
                # do gcd
                addrow(start,i,One(R));
                cleanout();
              fi;
            od;
          od;
        until alldivide;

      fi;

      # normalize
      qr:=StandardAssociateUnit(R,M[start][start]);
      multcol(start,qr);

    fi;
    start:=start+1;
  od;

  if transform then
   M:=rec(rowtrans:=left,coltrans:=right,normal:=M);
   if basmat<>fail then
     M.basmat:=basmat;
   fi;
   return M;
  else
    return M;
  fi;
end);

DeclareGlobalFunction("RationalCanonicalFormTransform");

InstallGlobalFunction(RationalCanonicalFormTransform,function(mat)
local cr,R,x,com,nf,matt,p,i,j,di,d,v;
  matt:=TransposedMat(mat);
  cr:=DefaultFieldOfMatrix(mat);
  R:=PolynomialRing(cr,1);
  x:=IndeterminatesOfPolynomialRing(R)[1];
  com:=x*mat^0-mat;
  com:=List(com,ShallowCopy);
  nf:=DoDiagonalizeMat2(R,com,true,true,matt);
  di:=DiagonalOfMat(nf.normal);
  p:=[];
  for i in [1..Length(di)] do
    d:=DegreeOfUnivariateLaurentPolynomial(di[i]);
    if d>0 then
      v:=List(nf.basmat[i],x->Value(x,Zero(cr))); # move in base ring
      Add(p,v);
      for j in [1..d-1] do
        v:=v*matt;
	Add(p,v);
      od;
    fi;
  od;
  return TransposedMat(p);
end);
