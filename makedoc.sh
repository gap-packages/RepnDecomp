#!/usr/bin/env bash
#
# Copyright (C) 2017-2019 Max Horn
# Copyright (C) 2017-2019 The GAP Team
# 
# This code is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# SPDX-License-Identifier: GPL-2.0-or-later
#
set -ex

GAPROOT=${GAPROOT-$HOME/gap}

# set up a custom GAP root containing only this package, so that
# we can force GAP to load the correct version of this package
# (we already did that in build_pkg.sh, but we do it again here,
# to allow the occasional instance where a package wants to also
# run the tests of others packages, by invoking this script multiple
# times in different directories)
mkdir -p gaproot/pkg/
ln -f -s $PWD gaproot/pkg/

# start GAP with custom GAP root, to ensure correct package version is loaded
GAP="$GAPROOT/bin/gap.sh -l $PWD/gaproot; --quitonbreak"

# Unless explicitly turned off by setting the NO_COVERAGE environment variable,
# we collect coverage data
if [[ -z $NO_COVERAGE ]]; then
    mkdir -p ${COVDIR-coverage}
    GAP="$GAP --cover ${COVDIR-coverage}/$(mktemp XXXXXX).coverage"
fi

# TODO: honor TestFile from PackageInfo record, but make sure that it
# is for the package in the current directory
$GAP -q < makedoc.g
