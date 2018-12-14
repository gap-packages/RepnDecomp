#!/bin/sh

gap --nointeract makedoc.g
rm -r public
mv doc public
(cd public; rm -f *.pdf *.lab *.txt *.six *.xml RepnDecomp.*)
