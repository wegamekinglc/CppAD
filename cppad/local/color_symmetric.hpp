// $Id$
# ifndef CPPAD_COLOR_SYMMETRIC_INCLUDED
# define CPPAD_COLOR_SYMMETRIC_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-15 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*!
\file color_symmetric.hpp
Coloring algorithm for a symmetric sparse matrix.
*/
// --------------------------------------------------------------------------
/*!
Determine which rows of a symmetric sparse matrix can be computed together
(using CppAD algorithm).

\tparam VectorSize
is a simple vector class with elements of type size_t.

\tparam VectorSet
is an unspecified type with the exception that it must support the
operations under pattern and the following operations where
p is a VectorSet object:
\n
<code>VectorSet p</code>
Constructs a new vector of sets object.
\n
<code>p.resize(ns, ne)</code>
resizes p to ns sets with elements between zero ne.
All of the ns sets are initially empty.
\n
<code>p.add_element(s, e)</code>
add element e to set with index s.

\param pattern [in]
Is a representation of the sparsity pattern for the matrix.
Note that color_symmetric does not change the values in pattern,
but it is not const because its iterator facility modifies some of its
internal data.
\n
<code>m = pattern.n_set()</code>
\n
sets m to the number of rows (and columns) in the sparse matrix.
All of the row indices are less than this value. 
\n
<code>n = pattern.end()</code>
\n
sets n to the number of columns in the sparse matrix
(which must be equal to the number of rows).
All of the column indices are less than this value. 
\n
<code>pattern.begin(i)</code>
instructs the iterator facility to start iterating over
columns in the i-th row of the sparsity pattern.
\n
<code>j = pattern.next_element()</code>
Sets j to the next possibly non-zero column 
in the row specified by the previous call to <code>pattern.begin</code>.
If there are no more such columns, the value
<code>pattern.end()</code> is returned.

\param row [in]
is a vector specifying which row indices to compute.

\param col [in]
is a vector, with the same size as row,
that specifies which column indices to compute.
For each  valid index k, the index pair
<code>(row[k], col[k])</code> must be present in the sparsity pattern.
It may be that some entries in the sparsity pattern do not need to be computed;
i.e, do not appear in the set of
<code>(row[k], col[k])</code> entries.

\param color [out]
is a vector with size m.
The input value of its elements does not matter.
Upon return, it is a coloring for the rows of the sparse matrix.
\n
\n
Fix any (i, j) in the sparsity pattern.
Furthermore suppose that there is an i1 with
i1 != i, 
color[i1] == color[i] and 
(i1, j) is in the sparsity pattern.
If follows that for all j1 with
j1 != j and color[j1] == color[j],
(j1, i ) is not in the sparsity pattern.
\n
\n
This routine tries to minimize, with respect to the choice of colors,
the maximum, with respct to k, of <code>color[ row[k] ]</code>,
not counting the indices where row[k] == m.
*/
template <class VectorSet>
void color_symmetric_cppad(
	VectorSet&              pattern   ,
	CppAD::vector<size_t>&  row       ,
	CppAD::vector<size_t>&  col       ,
	CppAD::vector<size_t>&  color     )
{	size_t i, j, k, ell;

	size_t K = row.size();
	size_t m = pattern.n_set();
	CPPAD_ASSERT_UNKNOWN( m == pattern.end() );
	CPPAD_ASSERT_UNKNOWN( color.size() == m );
	CPPAD_ASSERT_UNKNOWN( col.size()   == K );

	// row, column pairs that appear in ( row[k], col[k] )
	CppAD::vector< std::set<size_t> > pair_needed(m);
	std::set<size_t>::iterator itr, itr1;
	for(k = 0;  k < K; k++)
	{	CPPAD_ASSERT_UNKNOWN( pattern.is_element(row[k], col[k]) );
		pair_needed[ row[k] ].insert( col[k] );
		pair_needed[ col[k] ].insert( row[k] );
	}


	// initial coloring
	color.resize(m);
	ell = 0;
	for(i = 0; i < m; i++)
	{	if( pair_needed[i].empty() )
			color[i] = m;
		else
			color[i] = ell++;
	}

	// which colors are forbidden for this row
	CppAD::vector<bool> forbidden(m);

	// must start with row zero so that we remove results computed for it
	for(size_t i1 = 0; i1 < m; i1++) // for each row that appears
	if( color[i1] < m )
	{
		// initial all colors as ok for this row
		// (value of forbidden for ell > initial color[i] does not matter)
		for(ell = 0; ell <= color[i1]; ell++)
			forbidden[ell] = false;

		// -----------------------------------------------------
		// Forbid grouping with rows that would destroy results that are 
		// needed for this row.
		itr1 = pair_needed[i1].begin();
		while( itr1 != pair_needed[i1].end() )
		{	// entry (i1, j1) is needed for this row
			size_t j1 = *itr1;
			// Forbid rows i != i1 that have non-zero sparsity at (i, j1)
			// which is same as non-zero sparsity at (j1, i)
			pattern.begin(j1);
			i = pattern.next_element();
			while( i != pattern.end() )
			{	if( ( i < i1 ) & (color[i] < m) )
					forbidden[ color[i] ] = true;
				i = pattern.next_element();
			}
			itr1++;
		}


		// Forbid grouping with rows that this row would destroy results for
		for( i = 0; i < i1; i++)
		{	itr = pair_needed[i].begin();
			while( itr != pair_needed[i].end() )
			{	j = *itr;
				// (i, j) has is needed for row i
				// Forbid grouping with i1 if (i1, j) has non-zero sparsity
				if( pattern.is_element(i1, j) )
					forbidden[ color[i] ] = true;
				itr++;
			}
		}

		// pick the color with smallest index
		ell = 0;
		while( forbidden[ell] )
		{	ell++;
			CPPAD_ASSERT_UNKNOWN( ell <= color[i] );
		}
		color[i] = ell;

		// no longer need results that are computed by this row
		itr1 = pair_needed[i1].begin();
		while( itr1 != pair_needed[i1].end() )
		{	j = *itr1;
			if( j > i1 )
			{	itr = pair_needed[j].find(i1);
				if( itr != pair_needed[j].end() )
					pair_needed[j].erase(itr);
			}
			itr1++;
		}
	}

	// determine which sparsity entries need to be reflected 
	for(k = 0; k < row.size(); k++)
	{	i   = row[k];
		j   = col[k];
		itr = pair_needed[i].find(j);
		if( itr == pair_needed[i].end() )
		{	row[k] = j;
			col[k] = i;
# ifndef NDEBUG
			itr = pair_needed[j].find(i);
			CPPAD_ASSERT_UNKNOWN( itr != pair_needed[j].end() );
# endif
		}
	}
	return;
}

# if CPPAD_HAS_COLPACK
// --------------------------------------------------------------------------
/*!
Determine which rows of a symmetric sparse matrix can be computed together; 
i.e., do not have non-zero overlapping values, or can be computed by the
reflected entry.

\tparam VectorSize
is a simple vector class with elements of type \c size_t.

\tparam VectorSet
is an unspecified type with the exception that it must support the
operations under \c pattern and the following operations where
\c p is a \c VectorSet object:
\n
<code>VectorSet p</code>
Constructs a new vector of sets object.
\n
<code>p.resize(ns, ne)</code>
resizes \p to \c ns sets with elements between zero \c ne.
All of the \c ns sets are initially empty.
\n
<code>p.add_element(s, e)</code>
add element \c e to set with index \c s.

\param pattern [in]
Is a representation of the sparsity pattern for the matrix.
Note that \c color_general does not change the values in \c pattern,
but it is not \c const because its iterator facility modifies some of its
internal data.
\n
<code>m = pattern.n_set()</code>
\n
sets \c m to the number of rows in the sparse matrix.
All of the row indices are less than this value. 
\n
<code>n = pattern.end()</code>
\n
sets \c n to the number of columns in the sparse matrix.
All of the column indices are less than this value. 
\n
<code>pattern.begin(i)</code>
instructs the iterator facility to start iterating over
columns in the i-th row of the sparsity pattern.
\n
<code>j = pattern.next_element()</code>
Sets \c j to the next possibly non-zero column 
in the row specified by the previous call to <code>pattern.begin</code>.
If there are no more such columns, the value
<code>pattern.end()</code> is returned.

\param row [in/out]
is a vector specifying which row indices to compute.

\param col [in/out]
is a vector, with the same size as \c row,
that specifies which column indices to compute.
\n
\par Input
For each  valid index \c k, the index pair
<code>(row[k], col[k])</code> must be present in the sparsity pattern.
It may be that some entries in the sparsity pattern do not need to be computed;
i.e, do not appear in the set of
<code>(row[k], col[k])</code> entries.
\n
\par Output
On output, some of row and column indices may have been swapped
\code
	std::swap( row[k], col[k] )
\endcode
So the the the color for row[k] can be used to compute entry
(row[k], col[k]).

\param color [out]
is a vector with size \c m.
The input value of its elements does not matter.
Upon return, it is a coloring for the rows of the sparse matrix.
\n
\n
If for come \c i, <code>color[i] == m</code>, then 
the i-th row does not appear in the vector \c row.
Otherwise, <code>color[i] < m</code>.
\n
\n
Suppose two different rows, <code>i != r</code> have the same color and
column index \c j is such that both of the pairs 
<code>(i, j)</code> and <code>(r, j)</code> appear in the sparsity pattern.
It follows that neither of these pairs appear in the set of
<code>(row[k], col[k])</code> entries.
\n
\n
This routine tries to minimize, with respect to the choice of colors,
the maximum, with respct to \c k, of <code>color[ row[k] ]</code>.
*/
template <class VectorSet>
void color_symmetric_colpack(
	VectorSet&              pattern   ,
	CppAD::vector<size_t>&  row       ,
	CppAD::vector<size_t>&  col       ,
	CppAD::vector<size_t>&  color     )
{	size_t i, j, k;	
	size_t m = pattern.n_set();
	CPPAD_ASSERT_UNKNOWN( m == pattern.end() );
	CPPAD_ASSERT_UNKNOWN( row.size() == col.size() );

	// Determine number of non-zero entries in each row
	CppAD::vector<size_t> n_nonzero(m);
	size_t n_nonzero_total = 0;
	for(i = 0; i < m; i++)
	{	n_nonzero[i] = 0;
		pattern.begin(i);
		j = pattern.next_element();
		while( j != pattern.end() )
		{	n_nonzero[i]++;
			j = pattern.next_element();
		}
		n_nonzero_total += n_nonzero[i];
	}

	// Allocate memory and fill in Adolc sparsity pattern
	CppAD::vector<unsigned int*> adolc_pattern(m);
	CppAD::vector<unsigned int>  adolc_memory(m + n_nonzero_total);
	size_t i_memory = 0;
	for(i = 0; i < m; i++)
	{	adolc_pattern[i]    = adolc_memory.data() + i_memory;
		adolc_pattern[i][0] = n_nonzero[i];
		pattern.begin(i);
		j = pattern.next_element();
		k = 1;
		while(j != pattern.end() )
		{	adolc_pattern[i][k++] = j;
			j = pattern.next_element();
		}
		CPPAD_ASSERT_UNKNOWN( k == 1 + n_nonzero[i] );
		i_memory += k;
	}
	CPPAD_ASSERT_UNKNOWN( i_memory == m + n_nonzero_total );

	// Must use an external routine for this part of the calculation because
	// ColPack/ColPackHeaders.h has as 'using namespace std' at global level.
	cppad_colpack_symmetric(color, m, adolc_pattern);

	// determine which sparsity entries need to be reflected 
	size_t i1, i2, j1, j2, k1, k2;
	for(k1 = 0; k1 < row.size(); k1++)
	{	i1 = row[k1];
		j1 = col[k1];
		bool reflect = false;
		for(i2 = 0; i2 < m; i2++) if( (i1 != i2) & (color[i1]==color[i2]) )
		{	for(k2 = 1; k2 <= adolc_pattern[i2][0]; k2++)
			{	j2 = adolc_pattern[i2][k2];	
				reflect |= (j1 == j2);
			}
		}
		if( reflect )
		{	row[k1] = j1;
			col[k1] = i1;
		}
	}
	return;
}
# endif // CPPAD_HAS_COLPACK

# endif
