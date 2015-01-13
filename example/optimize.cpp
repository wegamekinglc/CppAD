/* $Id$ */
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-15 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin optimize.cpp$$

$section ADFun Operation Sequence Optimization: Example and Test$$

$index optimize, operation sequence$$
$index operation, optimize sequence$$
$index sequence, optimize operation$$
$index test, optimize$$
$index example, optimize$$

$code
$verbatim%example/optimize.cpp%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
namespace {
	template <class Float>
	void fun(const Float& x, Float& y, size_t& n_var, size_t& n_opt)
	{	using CppAD::exp;

	 	// Create a variable that is is only used in the comparision operation
		// (was optimized out until 2015-01-12).
		Float a = 1. / x;
		n_var += 1;
		n_opt += 1;

		// Create a variable that is used by the result
		Float b = x * 5.;
		n_var += 1;
		n_opt += 1;

		// only one variable created for this comparison operation
		// but the value depends on which branch is taken.
		Float c;
		if( a < x )  
			c = 2.0 * b;
		else
			c = 3.0 * b;
		n_var += 1;
		n_opt += 1;

		// Create a variable that is optimized out because it
		// will always have the same value as b
		Float d = 5. * x;
		n_var += 1;
		n_opt += 0;

		// Create three variables that will be converted to one
		// cumulative summation. Note that a is not connected to 
		// the result y (in the operation sequence).
		y      = 1. + b + c + d; 
		n_var += 3;
		n_opt += 1;
	}
}

bool optimize(void)
{	bool ok = true;
	using CppAD::AD;

	// domain space vector
	size_t n  = 1;
	CPPAD_TESTVECTOR(AD<double>) ax(n);
	ax[0]      = .5; 

	// declare independent variables and start tape recording
	CppAD::Independent(ax);
	size_t n_var = 1 + n; // one phantom variable at the beginning
	size_t n_opt = 1 + n; // and one for each independent variable

	// range space vector 
	size_t m = 1;
	CPPAD_TESTVECTOR(AD<double>) ay(m);
	fun(ax[0], ay[0], n_var, n_opt);

	// create f: x -> y and stop tape recording
	CppAD::ADFun<double> f(ax, ay);
	ok &= (f.size_var() == n_var);

	// Check zero order forward mode on the original operation sequence
	CPPAD_TESTVECTOR(double) x(n), y(m);
	x[0] = Value(ax[0]);
	y   = f.Forward(0, x);
	ok &= f.CompareChange() == 0;
	ok &= (y[0] == Value(ay[0]));

	// Optimize the operation sequence
	f.optimize();
	ok &= (f.size_var() == n_opt);

	// Check result for a zero order calculation.
	// This has already been checked if NDEBUG is not defined.
	y   = f.Forward(0, x);
	ok &= f.CompareChange() == 0;
	ok &= (y[0] == Value(ay[0]));

	// Check when the result of the comparision would be differnent
	x[0] = 2.0;
	y   = f.Forward(0, x);
	ok &= f.CompareChange() == 1;
	ok &= (y[0] != Value(ay[0]));

	return ok;
}

// END C++
