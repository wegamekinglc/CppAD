/* $Id$ */
# ifndef CPPAD_THREAD_TEAM_INCLUDED
# define CPPAD_THREAD_TEAM_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/* 
$begin team_thread.hpp$$
$spell
	pthreads
	const
	cstddef
	bool
	pthread
	initializes
	hpp
	num
	CppAD
	ta
$$
$section Specifications for A Team of AD Threads$$

$index thread, AD team$$
$index AD, thread team$$
$index team, AD threads$$

$head Syntax$$
$codei%include "team_thread.hpp"
%ok%   = team_start(%num_threads%)
%ok%   = team_work(%worker%)
%ok%   = team_stop()
%name% = team_name()
%$$

$head Purpose$$
These routines start, use, and stop a team of threads that can
be used with the CppAD type $code AD<double>$$.
For example,
these could be OpenMP threads, pthreads, or Boost threads to name a few.

$head Restrictions$$
Calls to the routines
$code team_start$$,
$code team_work$$, and
$code team_stop$$, must all be done by the master thread; i.e.,
$cref/thread_num/ta_thread_num/$$ must be zero.
In addition, they must all be done in sequential execution mode; i.e.,
when the master thread is the only thread that is running
($cref/in_parallel/ta_in_parallel/$$ must be false).

$head team_start$$
The argument 
$icode%num_threads% > 0%$$ has type $code size_t$$
and specifies the number of threads in this team.
This initializes both $code AD<double>$$ and $code team_work$$
to be used with $icode num_threads$$.
If $icode%num_threads% > 1%$$,
$icode%num_threads% - 1%$$ new threads are created
and put in a waiting state until $code team_work$$ is called.

$head team_work$$
This routine may be called one or more times
between the call to $code team_start$$ and $code team_stop$$.
The argument $icode worker$$ has type
$codei%bool %worker%(void)%$$.
Each call to $code team_work$$ runs $icode num_threads$$ versions
of $icode worker$$ with the corresponding value of
$cref/thread_num/ta_thread_num/$$ 
between zero and $icode%num_threads% - 1%$$ and
different for each thread,

$head team_stop$$
This routine terminates all the other threads except for
thread number zero; i.e., it terminates the threads corresponding to
$codei%
	%thread_num% = 1 , ... , %num_threads%-1
%$$

$head team_name$$
This routines returns a name that identifies this thread_team.
The return value has prototype
$icode%
	const char* %name%
%$$ 
and is a statically allocated $code '\0'$$ terminated C string.

$head ok$$
The return value $icode ok$$ has type $code bool$$.
It is $code false$$ if an error is detected during the
corresponding call.
Otherwise it is $code true$$.

$head Example Use$$
Example use of these specifications can be found in the file
$cref simple_ad.cpp$$.

$children%
	multi_thread/openmp/team_openmp.cpp%
	multi_thread/bthread/team_bthread.cpp%
	multi_thread/pthread/team_pthread.cpp
%$$
$head Example Implementation$$
Example implementations of these specifications can be found in the files:
$table
$rref team_openmp.cpp$$
$rref team_bthread.cpp$$
$rref team_pthread.cpp$$
$tend

$head Source$$
$codep */
# include <cstddef> // for size_t

extern bool team_start(size_t num_threads);
extern bool team_work(void worker(void));
extern bool team_stop(void);
extern const char* team_name(void);
/* $$
$end
*/

# endif