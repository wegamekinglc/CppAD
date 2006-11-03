# ifndef CPPAD_FUN_CONSTRUCT_INCLUDED
# define CPPAD_FUN_CONSTRUCT_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-06 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin FunConstruct$$
$spell 
	Taylor
	var
	ADvector
	const
$$

$spell
$$

$section Construct an ADFun Object and Stop Recording$$

$index ADFun, construct$$
$index construct, ADFun$$
$index tape, stop recording$$
$index recording, stop tape$$

$head Syntax$$
$syntax%ADFun<%Base%> %f%
%$$
$syntax%ADFun<%Base%> %f%(%x%, %y%)
%$$


$head Purpose$$
The $syntax%AD<%Base%>%$$ object $italic f$$ can 
store an AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$.
It can then be used to calculate derivatives of the corresponding
$xref/glossary/AD Function/AD function/$$
$latex \[
	F : B^n \rightarrow B^m
\] $$
where $latex B$$ is the space corresponding to objects of type $italic Base$$.

$head x$$
If the argument $italic x$$ is present,
it must be the vector argument in the previous call to
$cref/Independent/$$.
Neither its size, or any of its values, are allowed to change
between calling
$syntax%
	Independent(%x%)
%$$
and 
$syntax%
	ADFun<%Base%> %f%(%x%, %y%)
%$$.

$head y$$
The sequence of operations that map $italic x$$
to $italic y$$ are stored in the AD function object $italic f$$.


$head Default Constructor$$
The default constructor 
$syntax%
	ADFun<%Base%> %f%
%$$
creates an 
$syntax%AD<%Base%>%$$ object with no corresponding operation sequence; i.e.,
$syntax%
	%f%.size_var()
%$$
returns the value zero (see $xref/SeqProperty/size_var/size_var/$$).

$head Sequence Constructor$$
The sequence constructor 
$syntax%
	ADFun<%Base%> %f%(%x%, %y%)
%$$
creates the $syntax%AD<%Base%>%$$ object $italic f$$,
stops the recording of AD of $italic Base$$ operations,
stores the corresponding operation sequence in the object $italic f$$.
It then stores the first order Taylor coefficients 
(corresponding to the value of $italic x$$) in $italic f$$.
This is equivalent to the following steps using the default constructor:
$list number$$
Create $italic f$$ with the default constructor
$syntax%
	ADFun<%Base%> %f%;
%$$
$lnext
Stop the tape and storing the operation sequence using
$syntax%
	%f%.Dependent(%y%);
%$$
(see $xref/Dependent/$$).
$lnext
Calculating the first order Taylor coefficients for all 
the variables in the operation sequence using
$syntax%
	%f%.Forward(%p%, %x_p%)
%$$
with $italic p$$ equal to zero and the elements of $italic x_p$$
equal to the corresponding elements of $italic x$$
(see $xref/Forward/$$).
$lend

$head Example$$

$subhead Sequence Constructor$$
The file
$xref/Independent.cpp/$$ 
contains an example and test of the sequence constructor.
It returns true if it succeeds and false otherwise.

$subhead Default Constructor$$
The files
$xref/FunCheck.cpp/$$ 
and
$xref/HesLagrangian.cpp/$$
contain an examples and tests using the default constructor.
They return true if they succeed and false otherwise.

$end
*/


// BEGIN CppAD namespace
namespace CppAD {

template <typename Base>
template <typename VectorAD>
ADFun<Base>::ADFun(const VectorAD &x, const VectorAD &y)
: totalNumVar(0), Taylor(CppADNull), ForJac(CppADNull)
{	size_t i, j, m, n;

	// stop the tape and store the operation sequence
	Dependent(y);

	// allocate memory for one zero order Taylor coefficient
	taylor_per_var= 1;
	TaylorColDim  = 1;
	Taylor        = CppADTrackNewVec(totalNumVar, Taylor);

	// set zero order coefficients corresponding to indpendent variables
	n = ind_taddr.size();
	CppADUsageError(
		n == x.size(),
		"ADFun<Base>: independent variable vector has changed"
	);
	for(j = 0; j < n; j++)
	{	CppADUnknownError( ind_taddr[j] == (j+1) );
		CppADUsageError(
			x[j].taddr == (j+1),
			"ADFun<Base>: independent variable vector has changed"
		);
		Taylor[ ind_taddr[j] ]  = x[j].value;
	}

	// use independent variable values to fill in values for others
	compareChange = ForwardSweep(
		false, 0, totalNumVar, &Rec, TaylorColDim, Taylor
	);
	CppADUnknownError( compareChange == 0 );

	// check the dependent variable values
	m = dep_taddr.size();
	for(i = 0; i < m; i++) CppADUsageError(
		Taylor[dep_taddr[i]] == y[i].value,
		"independent variable not equal its tape evaluation"
		", it may be nan"
	);
}

} // END CppAD namespace

# endif