LoadPackage("RepnDecomp");

BenchMany@RepnDecomp(rep -> REPN_ComputeUsingMyMethod(rep.rep : irreps := rep.irreps), "fast.txt", 100);
