/* $Id$ */
# ifndef CPPAD_COMP_OP_INCLUDED
# define CPPAD_COMP_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-15 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */


namespace CppAD { // BEGIN_CPPAD_NAMESPACE
/*!
\file comp_op.hpp
Zero order forward mode check how many comparisons changed.
*/

// -------------------------------- <= -----------------------------------
template <class Base>
inline void forward_leqpv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LeqpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LeqpvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base  x = parameter[ arg[0] ];
	Base* y = taylor + arg[1] * cap_order;

	count += GreaterThanZero(x - y[0]);
}
template <class Base>
inline void forward_leqvp_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LeqvpOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LeqvpOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base  y = parameter[ arg[1] ];

	count += GreaterThanZero(x[0] - y);
}
template <class Base>
inline void forward_leqvv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LeqvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LeqvvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base* y = taylor + arg[1] * cap_order;

	count += GreaterThanZero(x[0] - y[0]);
}
// ------------------------------- < -------------------------------------
template <class Base>
inline void forward_ltpv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LtpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LtpvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base  x = parameter[ arg[0] ];
	Base* y = taylor + arg[1] * cap_order;

	count += GreaterThanOrZero(x - y[0]);
}
template <class Base>
inline void forward_ltvp_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LtvpOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LtvpOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base  y = parameter[ arg[1] ];

	count += GreaterThanOrZero(x[0] - y);
}
template <class Base>
inline void forward_ltvv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(LtvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(LtvvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base* y = taylor + arg[1] * cap_order;

	count += GreaterThanOrZero(x[0] - y[0]);
}
// ------------------------------ == -------------------------------------
template <class Base>
inline void forward_eqpv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(EqpvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(EqpvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base  x = parameter[ arg[0] ];
	Base* y = taylor + arg[1] * cap_order;

	count += (x != y[0]);
}
template <class Base>
inline void forward_eqvv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(EqvvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(EqvvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base* y = taylor + arg[1] * cap_order;

	count += (x[0] != y[0]);
}
// -------------------------------- != -----------------------------------
template <class Base>
inline void forward_nepv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(NepvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(NepvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base  x = parameter[ arg[0] ];
	Base* y = taylor + arg[1] * cap_order;

	count += (x == y[0]);
}
template <class Base>
inline void forward_nevv_op_0(
	size_t&       count       ,
	const addr_t* arg         ,
	const Base*   parameter   ,
	size_t        cap_order   ,
	Base*         taylor      )
{
	// check assumptions
	CPPAD_ASSERT_UNKNOWN( NumArg(NevvOp) == 2 );
	CPPAD_ASSERT_UNKNOWN( NumRes(NevvOp) == 0 );

	// Taylor coefficients corresponding to arguments and result
	Base* x = taylor + arg[0] * cap_order;
	Base* y = taylor + arg[1] * cap_order;

	count += (x[0] == y[0]);
}


} // END_CPPAD_NAMESPACE
# endif
