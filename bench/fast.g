Read("common.g");

start := Runtime();;
x:=BlockDiagonalRepresentationFast(rho, irreps);;
Print("RUNTIME: ", Runtime()-start, "\n");;
