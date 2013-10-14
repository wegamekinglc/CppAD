#! /bin/bash -e
# $Id$
# -----------------------------------------------------------------------------
# CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-13 Bradley M. Bell
#
# CppAD is distributed under multiple licenses. This distribution is under
# the terms of the 
#                     Eclipse Public License Version 1.0.
#
# A copy of this license is included in the COPYING file of this distribution.
# Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
# -----------------------------------------------------------------------------
# $begin get_adolc.sh$$ $newlinech #$$
# $spell
#	tgz
#	Adolc
#	gz
#	CppAD
# $$
#
# $section Download and Install Adolc in Build Directory$$
# $index adolc, download and install$$
# $index download, install adolc$$
# $index install, adolc$$ 
#
# $head Syntax$$
# $code bin/get_adolc.sh$$
#
# $head Purpose$$
# If you are using Unix, this command will download and install 
# $href%https://projects.coin-or.org/ADOL-C%ADOL-C%$$ in the
# CppAD $code build$$ directory.
#
# $children%
#	bin/get_colpack.sh
# %$$
# $head Requirements$$
# You must first use $cref get_colpack.sh$$ to download and install
# $code ColPack$$ (coloring algorithms used for sparse matrix derivatives).
#
# $head Distribution Directory$$
# This command must be executed in the 
# $cref/distribution directory/download/Distribution Directory/$$.
#
# $head External Directory$$
# The Adolc source code is downloaded into the sub-directory
# $code build/external$$ below the distribution directory.
#
# $head Prefix Directory$$
# The Adolc include files are installed in the sub-directory
# $code build/prefix/include/adolc$$ below the distribution directory.
#
# $head Reuse$$
# The files $codei%build/external/ADOL-C-%version%.tgz%$$
# and the directory $codei%build/external/ADOL-C-%version%$$
# will be reused if they exist. Delete this file and directory
# to get a complete rebuild.
#
# $end
# -----------------------------------------------------------------------------
if [ $0 != "bin/get_adolc.sh" ]
then
	echo "bin/get_adolc.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
echo 'Download adolc to build/external and install it to build/prefix'
version='2.4.1'
web_page="http://www.coin-or.org/download/source/ADOL-C"
prefix=`pwd`'/build/prefix'
# --------------------------------------------------------------------------
if [ -e /usr/lib64 ]
then
	libdir='lib64'
else
	libdir='lib'
fi
# -----------------------------------------------------------------------------
if [ ! -d build/external ]
then
	echo_eval mkdir -p build/external
fi
echo_eval cd build/external
# -----------------------------------------------------------------------------
if [ ! -e "ADOL-C-$version.tgz" ]
then
	echo_eval wget --no-check-certificate $web_page/ADOL-C-$version.tgz
fi
# -----------------------------------------------------------------------------
if [ -e "$prefix/include/adolc" ]
then
	echo_eval rm -r "$prefix/include/adolc"
fi
# -----------------------------------------------------------------------------
if [ ! -e ADOL-C-$version ]
then
	echo_eval tar -xzf ADOL-C-$version.tgz
fi
echo_eval cd ADOL-C-$version
# -----------------------------------------------------------------------------
if [ ! -e build ]
then
	echo_eval mkdir build
fi
echo_eval cd build
# -----------------------------------------------------------------------------
flags="--prefix=$prefix --with-colpack=$prefix --libdir=$prefix/$libdir"
#
echo "../configure $flags" 
../configure $flags
#
echo "make install"
make install
# -----------------------------------------------------------------------------
echo "get_adolc: OK"