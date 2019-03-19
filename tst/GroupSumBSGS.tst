gap> tester := rep -> GroupSumBSGS(rep.rep) = Sum(Source(rep.rep), g -> Image(rep.rep, g));;
gap> TestMany@RepnDecomp(tester, 3);
true