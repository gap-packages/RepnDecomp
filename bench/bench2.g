LoadPackage("RepnDecomp");

BenchMany@RepnDecomp(rep -> BlockDiagonalRepresentation(rep.rep), "serre.txt", 100);
