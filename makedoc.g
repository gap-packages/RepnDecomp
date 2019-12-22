#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# This file is a script which compiles the package manual.
#
if fail = LoadPackage("AutoDoc", "2016.02.16") then
    Error("AutoDoc version 2016.02.16 or newer is required.");
fi;

AutoDoc( rec( extract_examples := true, scaffold := true, autodoc := true, dir := "doc/" ) );
