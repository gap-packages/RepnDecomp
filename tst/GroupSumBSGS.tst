gap> tester := rep -> GroupSumBSGS(Source(rep.rep), g -> Image(rep.rep, g)) = Sum(Source(rep.rep), g -> Image(rep.rep, g));;
gap> TestMany@RepnDecomp(tester, 3);
true