Read("common.g");

start := Runtime();;
x:=BlockDiagonalRepresentationFast(rho);;
Print("RUNTIME: ", Runtime()-start, "\n");;
