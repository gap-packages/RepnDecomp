LoadPackage("RepnDecomp");

BenchMany@RepnDecomp(rep -> REPN_ComputeUsingSerre(rep.rep : irreps := rep.irreps), "serre.txt", 100);
