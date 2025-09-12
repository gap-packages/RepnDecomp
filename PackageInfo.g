#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "RepnDecomp",
Subtitle := "Decompose representations of finite groups into irreducibles",
Version := "1.3.1",
Date := "10/09/2025", # dd/mm/yyyy format
License := "GPL-3.0",

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Kaashif",
    LastName := "Hymabaccus",
    WWWHome := "https://kaashif.co.uk",
    Email := "kaashif@kaashif.co.uk",
    Place := "Oxford",
    Institution := "University of Oxford",
  ),
  rec(
    IsAuthor := false,
    IsMaintainer := true,
    FirstNames := "Dmitrii",
    LastName := "Pasechnik",
    WWWHome := "https://pasechnik.info/dima",
    Email := "dima@pasechnik.info",
    Place := "Evanston",
    Institution := "Northwestern University",
  ),
],

SourceRepository := rec( Type := "git",
                         URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
                       ),

IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats  := ".tar.gz .tar.bz2",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "deposited",

AbstractHTML := "The <span class='pkgname'>RepnDecomp</span> package provides functions implementing various algorithms for decomposing linear representations of finite groups.",

PackageDoc := rec(
  BookName  := "RepnDecomp",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Decompose representations of finite groups into irreducibles",
),

Dependencies := rec(
  GAP := ">= 4.10",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.6.1" ],
                           [ "GRAPE", ">= 4.8.1" ] ],
  SuggestedOtherPackages := [ [ "IO", ">= 4.7.0" ] ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [ "representations", "groups", ],

));
