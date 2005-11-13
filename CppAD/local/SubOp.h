# ifndef CppADSubOpIncluded
# define CppADSubOpIncluded

// BEGIN SHORT COPYRIGHT
/* -----------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-05 Bradley M. Bell

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
------------------------------------------------------------------------ */
// END SHORT COPYRIGHT

/*
$begin ForSubvvOp$$
$spell
	Subpv
	Subvp
	Subvv
	Taylor
	const
	inline
	Op
$$


$index subtract, forward operator$$
$index forward, subtract operator$$
$index operator, subtract forward$$
$index ForSub$$

$section Forward Mode Subtraction Operators$$

$table
$bold Syntax$$ 
$cnext
$syntax%inline void ForSubvvOp(size_t %d%, 
	%Base% *%z%, const %Base% *%x%, const %Base% *%y%)%$$
$rnext $cnext
$syntax%inline void ForSubpvOp(size_t %d%, 
	%Base% *%z%, const %Base% *%p%, const %Base% *%y%)%$$
$rnext $cnext
$syntax%inline void ForSubvpOp(size_t %d%, 
	%Base% *%z%, const %Base% *%x%, const %Base% *%p%)%$$
$tend

$fend 25$$

$head Description$$
Computes the $italic d$$ order Taylor coefficient for $latex Z$$ where
$table
Operation  $cnext Value  $rnext
Subvv       $cnext $latex Z = X - Y$$ $rnext
Subpv       $cnext $latex Z = P - Y$$ $rnext
Subvp       $cnext $latex Z = X - P$$ 
$tend

$head x$$
The vector $italic x$$ has length $latex d+1$$ and contains the
$th d$$ order Taylor coefficient row vector for $italic X$$.

$head y$$
The vector $italic y$$ has length $latex d+1$$ and contains the
$th d$$ order Taylor coefficient row vector for $italic Y$$.

$head p$$
The scalar $syntax%*%p%$$ contains the value of the parameter $italic P$$.

$head z$$
The vector $italic z$$ has length $latex d+1$$.
On input it contains the
$th d-1$$ order Taylor coefficient row vector for $italic Z$$.
On output it contains the
$th d$$ order Taylor coefficient row vector for $italic Z$$; i.e.,
$syntax%%z%[%d%]%$$ is set equal to the $th d$$ Taylor coefficient for
the function $italic Z$$.

$end
------------------------------------------------------------------------------
$begin RevSubvvOp$$
$spell
	Subpv
	Subvp
	Subvv
	Taylor
	const
	inline
	Op
	px
	py
	pz
$$

$mindex RevSubvvOp reverse minus subtract$$
$section Reverse Mode Subtraction Operator$$

$table
$bold Syntax$$ 
$cnext
$syntax%inline void RevSubvvOp(size_t %d%, 
	const %Base% *%pz%, %Base% *%px%, %Base% *%py%)%$$
$cnext
$syntax%inline void RevSubpvOp(size_t %d%, 
	const %Base% *%pz%, %Base% *%py%)%$$
$cnext
$syntax%inline void RevSubvpOp(size_t %d%, 
	const %Base% *%pz%, %Base% *%px%)%$$
$tend

$head Description$$
We are given the partial derivatives for a function
$latex G(z, x, y)$$ and we wish to compute the partial derivatives for
the function
$latex \[
	H(x, y) = G [ Z(x, y) , x , y ]
\]$$
where $latex Z(x, y)$$ is defined as the 
$th d$$ order Taylor coefficient row vector for $italic Z$$ as
a function of the corresponding vectors for 
$italic X$$ and $italic Y$$ where

$table
Operation  $cnext Value  $rnext
Subvv       $cnext $latex Z = X - Y$$ $rnext
Subpv       $cnext $latex Z = P - Y$$ $rnext
Subvp       $cnext $latex Z = X - P$$ 
$tend

Note that $italic Z$$ has been used both the original subtraction 
function and for the corresponding mapping of Taylor coefficients.

$head pz$$
The vector $italic pz$$ has length $latex d+1$$ and 
$syntax%%pz%[%j%]%$$ contains the partial for $italic G$$
with respect to the $th j$$ order Taylor coefficient for $italic Z$$.

$head On Input$$

$subhead px$$
The vector $italic px$$ has length $latex d+1$$ and 
$syntax%%px%[%j%]%$$ contains the partial for $italic G$$
with respect to the $th j$$ order Taylor coefficient for $italic X$$.

$subhead py$$
The vector $italic py$$ has length $latex d+1$$ and 
$syntax%%py%[%j%]%$$ contains the partial for $italic G$$
with respect to the $th j$$ order Taylor coefficient for $italic Y$$.

$head On Output$$

$subhead px$$
If present,
the vector $italic px$$ has length $latex d+1$$ and 
$syntax%%px%[%j%]%$$ contains the partial for $italic H$$
with respect to the $th j$$ order Taylor coefficient for $italic X$$.

$subhead py$$
If present,
the vector $italic py$$ has length $latex d+1$$ and 
$syntax%%py%[%j%]%$$ contains the partial for $italic H$$
with respect to the $th j$$ order Taylor coefficient for $italic Y$$.


$end
------------------------------------------------------------------------------
*/

// BEGIN CppAD namespace
namespace CppAD {

// --------------------------- Subvv -----------------------------------------

template <class Base>
inline void ForSubvvOp(size_t d, 
	Base *z, const Base *x, const Base *y)
{
	z[d] = x[d] - y[d];
}

template <class Base>
inline void RevSubvvOp(size_t d, 
	const Base *pz, Base *px, Base *py)
{
	// number of indices to access
	size_t i = d + 1;

	while(i)
	{	--i;
		px[i] += pz[i];
		py[i] -= pz[i];
	}
}

// --------------------------- Subpv -----------------------------------------

template <class Base>
inline void ForSubpvOp(size_t d, 
	Base *z, const Base *p, const Base *y)
{
	if( d == 0 )
		z[d] = (*p) - y[d];
	else	z[d] = - y[d];

}

template <class Base>
inline void RevSubpvOp(size_t d, 
	const Base *pz, Base *py)
{
	// number of indices to access
	size_t i = d + 1;

	while(i)
	{	--i;
		py[i] -= pz[i];
	}
}

// --------------------------- Subvp -----------------------------------------

template <class Base>
inline void ForSubvpOp(size_t d, 
	Base *z, const Base *x, const Base *p)
{
	if( d == 0 )
		z[d] = x[d] - (*p);
	else	z[d] = x[d];

}

template <class Base>
inline void RevSubvpOp(size_t d, 
	const Base *pz, Base *px)
{
	// number of indices to access
	size_t i = d + 1;

	while(i)
	{	--i;
		px[i] += pz[i];
	}
}

} // END CppAD namespace

# endif
