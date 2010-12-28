#! /bin/bash -e 
# $Id: dos_format.sh 1370 2009-05-31 05:31:50Z bradbell $
# -----------------------------------------------------------------------------
# CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-10 Bradley M. Bell
#
# CppAD is distributed under multiple licenses. This distribution is under
# the terms of the 
#                     Common Public License Version 1.0.
#
# A copy of this license is included in the COPYING file of this distribution.
# Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
# -----------------------------------------------------------------------------
version=`cat configure.ac | grep "^ *AC_INIT(" | 
        sed -e 's/[^,]*, *\([^ ,]*\).*/\1/'`
yyyy_mm_dd=`echo $version | sed -e 's/\(....\)\(..\)/\1-\2-/'`
# ---------------------------------------------------------------------
list="AUTHORS"
svn cat AUTHORS | sed > AUTHORS.$$ \
	-e "s/, [0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\} *,/, $yyyy_mm_dd,/"
# ---------------------------------------------------------------------
list="$list configure"
svn cat configure | sed > configure.$$ \
	-e "s/CppAD [0-9]\{8\}[.0-9]*/CppAD $version/g" \
	-e "s/VERSION='[0-9]\{8\}[.0-9]*'/VERSION='$version'/g" \
	-e "s/configure [0-9]\{8\}[.0-9]*/configure $version/g" \
	-e "s/config.status [0-9]\{8\}[.0-9]*/config.status $version/g" \
	-e "s/\$as_me [0-9]\{8\}[.0-9]*/\$as_me $version/g" \
	-e "s/Generated by GNU Autoconf.*$version/&./"
# ---------------------------------------------------------------------
list="$list configure.ac"
svn cat configure.ac | sed > configure.ac.$$\
     -e "s/(CppAD, [0-9]\{8\}[.0-9]* *,/(CppAD, $version,/"
# ---------------------------------------------------------------------
list="$list cppad/config.h"
# svn_commit.sed will make sure config.h has these values
sed -i.save cppad/config.h \
	-e 's/\(^# *define *CPPAD_BOOSTVECTOR *\) 1 *$/\1 0/' \
	-e 's/\(^# *define *CPPAD_CPPADVECTOR *\) 0 *$/\1 1/' \
	-e 's/\(^# *define *CPPAD_STDVECTOR *\) 1 *$/\1 0/'
#
svn cat cppad/config.h | sed > cppad/config.h.$$ \
	-e "s/CppAD [0-9]\{8\}[.0-9]*/CppAD $version/g" \
	-e "s/VERSION \"[0-9]\{8\}[.0-9]*\"/VERSION \"$version\"/g"
# ---------------------------------------------------------------------
list="$list cppad/configure.hpp"
svn cat cppad/configure.hpp | sed > cppad/configure.hpp.$$ \
	-e "s/CppAD [0-9]\{8\}[.0-9]*/CppAD $version/g"
# ---------------------------------------------------------------------
for file in $list
do
	if diff $file.$$ $file > /dev/null
	then
		printf "%-20s at most the version number is different.\n" $file:
	else
		printf "%-20s changes not counting version number:\n" $file:
		diff $file.$$ $file
	fi
done
mv cppad/config.h.save config.h
