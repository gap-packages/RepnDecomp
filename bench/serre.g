Read("common.g");

start := Runtime();;
x:=BlockDiagonalBasis(rho);;
c:=RepresentationCentralizerBlocks(rho);;
Print("RUNTIME: ", Runtime()-start, "\n");;
