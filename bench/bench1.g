LoadPackage("RepnDecomp");

BenchMany@RepnDecomp(rep -> BlockDiagonalRepresentationFast(rep.rep, rep.irreps), "fast.txt", 100);
