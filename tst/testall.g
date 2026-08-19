#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "RepnDecomp" );

# Shared deterministic property checks used by the hand-written regression
# tests.  The AutoDoc-generated examples remain independent of these helpers.
ReadPackage( "RepnDecomp", "tst/testutils.g" );

# Several public algorithms use randomized choices internally.  Fixing the
# global source makes failures reproducible while still exercising those paths.
Reset( GlobalMersenneTwister, 27041993 );

TestDirectory(DirectoriesPackageLibrary( "RepnDecomp", "tst" ),
              rec(exitGAP := true,
                  testOptions := rec(compareFunction := "uptowhitespace")));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
