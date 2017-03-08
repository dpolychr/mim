/**
    MIM
    Copyright (C) 2017 Lorraine A. K. Ayad, Chang Liu and Solon P. Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include "mim.h"

using namespace std;
using namespace seqan;

int find_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims )
{

	fprintf ( stderr, " -Merging exact matches\n" );
	merge( sw, ref, query, q_grams, mims );

	fprintf ( stderr, " -Extending merged matches\n" );
	for( int i=0; i<mims->size(); i++ )
	{
		if( mims->at(i). error < sw . k )
			extend( &mims->at(i).error, (int*) &mims->at(i).startQuery, (int*) &mims->at(i).endQuery, (int*) &mims->at(i).startRef, (int*) &mims->at(i).endRef, ref, query, sw );
	}
	
	fprintf ( stderr, " -Adjusting extended matches\n" );
	for( int j=0; j<mims->size(); j++ )
	{
		adjust(  &mims->at(j).error, (int*) &mims->at(j).startQuery, (int*) &mims->at(j).endQuery, (int*) &mims->at(j).startRef, (int*) &mims->at(j).endRef, ref, query, sw );
	}

return 0;
}

int merge( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims )
{
	for( int i = 0; i<q_grams->size(); i++ )
	{	
		unsigned int current_qgram = i;	
		unsigned int edit_distance = 0;

		unsigned int q_start = q_grams->at(i).occQuery;
		unsigned int q_end = q_start + q_grams->at(i).length ;
		unsigned int r_start = q_grams->at(i).occRef;
		unsigned int r_end = r_start + q_grams->at(i).length ;
		int gap_size_ref = 0;
		int gap_size_query = 0;

		for( int j = i + 1; j<q_grams->size(); j++ )
		{
		
			if( q_grams->at(j).occRef < q_grams->at(current_qgram).occRef )
				break;

			gap_size_ref = 	q_grams->at(j).occRef - ( q_grams->at(current_qgram).occRef + q_grams->at(current_qgram).length ); 
			gap_size_query = q_grams->at(j).occQuery - ( q_grams->at(current_qgram).occQuery + q_grams->at(current_qgram).length );

			if( gap_size_ref <= 0 && gap_size_query >= sw . k  )
			{
				if( edit_distance + gap_size_query <= sw . k )
				{
					edit_distance = edit_distance + gap_size_query;
					r_end = q_grams->at(j).occRef + q_grams->at(j).length; 
					q_end = q_end + gap_size_query + q_grams->at(j).length;
					current_qgram = j;
				}
				else break;
			}
			else if( gap_size_query <= 0 && gap_size_ref >= sw . k ) 
			{
	
				if( edit_distance + gap_size_ref <= sw . k )
				{

					edit_distance = edit_distance + gap_size_ref;
					r_end = r_end + gap_size_ref + q_grams->at(j).length;
					q_end =  q_grams->at(j).occQuery + q_grams->at(j).length; 
					current_qgram = j;
				}
				else break;
			}
			else if( gap_size_query <= 0 && gap_size_ref <= 0 )
			{	
				r_end = q_grams->at(j).occRef + q_grams->at(j).length;
				q_end = q_grams->at(j).occQuery + q_grams->at(j).length;
				current_qgram = j;
			}
			else if ( gap_size_query > 0 && gap_size_ref > 0 )
			{
				if( abs( gap_size_query -  gap_size_ref ) > sw . k )
					break;
				
				unsigned char * m_query = ( unsigned char * ) calloc ( gap_size_query + 1, sizeof ( unsigned char ) );
				unsigned char * m_ref = ( unsigned char * ) calloc ( gap_size_ref + 1, sizeof ( unsigned char ) );
			
				memcpy( &m_query[0], &query[ q_end ], gap_size_query );
				memcpy( &m_ref[0], &ref[ r_end ] , gap_size_ref );

				m_query[ gap_size_query ] = '\0';
				m_ref[ gap_size_ref ] = '\0';

				int matching_qgrams = compute_qgrams( m_ref, m_query );

				if( ( ( strlen( ( char * ) ref ) + 1 - matching_qgrams) / 3 ) - 1 + edit_distance > sw . k )
					break;

				int edit_distance_temp = edit_distance + editDistanceMyers( m_query, m_ref );


				free( m_query );
				free( m_ref );

				if( edit_distance_temp <= sw . k )
				{
					edit_distance = edit_distance_temp;
					r_end = r_end + gap_size_ref + q_grams->at(j).length; 
					q_end = q_end + gap_size_query + q_grams->at(j).length;
					current_qgram = j;
				
				}
				else break;
			}

		}	

		MimOcc occ;
		occ.startRef = r_start;
		occ.endRef = r_end;
		occ.startQuery = q_start;
		occ.endQuery = q_end;
		occ.error = edit_distance;
		mims->push_back(occ);

	}
	return 0;
}


int extend( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw )
{

	int toAddStartQuery = 1;
	int toAddEndQuery = 1;
	int toAddStartRef = 1;
	int toAddEndRef = 1;

	int q_start_temp = *q_start;
	int q_end_temp = *q_end;
	int r_start_temp = *r_start;
	int r_end_temp = *r_end;

	int qS = *q_start;
	int rS = *r_start;	
	int qE = *q_end;
	int rE = *r_end;

	int edit_distance_total_L = 0;
	int edit_distance_total_R = 0;

	int edit_distance_temp = *edit_distance;
	int edit_distance_updated = *edit_distance;

	while( q_start_temp >= 0 || r_start_temp >= 0 || q_end_temp <=  strlen( ( char* ) yInput ) -1 || r_end_temp <=  strlen( ( char* ) xInput ) )
	{

		int edit_distance_R;
		if (  q_end_temp  < strlen( ( char* ) yInput )  &&  r_end_temp  < strlen( ( char* ) xInput )  ) 
		{	
			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
			memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';
		
			edit_distance_R =  editDistanceMyers( m_ref_R, m_query_R );
			free( m_ref_R );
			free( m_query_R );

		}
		else if( qE == strlen( ( char* ) yInput ) && rE != strlen( ( char* ) xInput ) )
		{
			edit_distance_R = edit_distance_total_R + 1;
		}
		else if( rE == strlen( ( char* ) xInput ) && qE != strlen( ( char* ) yInput ) )
		{
			edit_distance_R = edit_distance_total_R + 1;
		}
		else if ( q_end_temp  < strlen( ( char* ) yInput ) && r_end_temp >= strlen( ( char* ) xInput ) )	
		{
			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  strlen( ( char* ) xInput ) - rE );
			memcpy( &m_query_R[0], &yInput[qE], toAddEndQuery  );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';
		
			edit_distance_R =  editDistanceMyers( m_ref_R, m_query_R );
			free( m_ref_R );
			free( m_query_R );
		}
		else if ( q_end_temp  >= strlen( ( char* ) yInput ) - 1 && r_end_temp < strlen( ( char* ) xInput ) - 1 )	
		{
			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef );
			memcpy( &m_query_R[0], &yInput[qE], strlen( ( char* ) yInput ) - qE   );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';
		
			edit_distance_R =  editDistanceMyers( m_ref_R, m_query_R );
			free( m_ref_R );
			free( m_query_R );
		}
		else edit_distance_R = sw . k + 1;

		int edit_distance_L;

		if(  q_start_temp  > 0 &&  r_start_temp > 0 ) 
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
				
			memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
			memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';

		
			edit_distance_L = editDistanceMyers( m_ref_L, m_query_L );
			free( m_ref_L );
			free( m_query_L );

		}
		else if( qS == 0 )
		{
			edit_distance_L = edit_distance_total_L + 1;
		}
		else if( rS == 0 )
		{
			edit_distance_L = edit_distance_total_L + 1;
		}
		else if ( q_start_temp  <= 0 && r_start_temp > 0 )	
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
				
			memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
			memcpy( &m_query_L[0], &yInput [0], qS );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';
			edit_distance_L = editDistanceMyers( m_ref_L, m_query_L );
			free( m_ref_L );
			free( m_query_L );
		}
		else if ( q_start_temp  > 0 && r_start_temp <= 0 )	
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
				
			memcpy( &m_ref_L[0], &xInput [0], rS );
			memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';

		
			edit_distance_L = editDistanceMyers( m_ref_L, m_query_L );
			free( m_ref_L );
			free( m_query_L );
		}
		else edit_distance_L = sw . k + 1;


		if( edit_distance_L + edit_distance_R + edit_distance_temp > sw . k )
		{

			if( edit_distance_temp + edit_distance_total_R + edit_distance_L < edit_distance_temp + edit_distance_R + edit_distance_total_L  && edit_distance_temp + edit_distance_total_R + edit_distance_L <= sw . k ) //extend left
			{
				toAddStartQuery++;
				q_start_temp--;
				toAddStartRef++;
				r_start_temp--;
				edit_distance_total_L = edit_distance_L;
				edit_distance_updated = edit_distance_temp + edit_distance_total_R + edit_distance_L;
			}
			else if ( edit_distance_temp + edit_distance_R + edit_distance_total_L <= edit_distance_temp + edit_distance_total_R + edit_distance_L && edit_distance_temp + edit_distance_R + edit_distance_total_L <= sw . k ) //extend right
			{
				toAddEndQuery++;
				q_end_temp++;
				toAddEndRef++;
				r_end_temp++;
				
				edit_distance_total_R = edit_distance_R;
				edit_distance_updated = edit_distance_temp + edit_distance_R + edit_distance_total_L;
			}
			else
			{
				if( q_start_temp < 0 )
					*q_start = 0;
				else *q_start = q_start_temp;

				if( r_start_temp < 0 )
					*r_start = 0;
				else *r_start = r_start_temp;

				if( q_end_temp > strlen( ( char* ) yInput ) )
					*q_end =  strlen( ( char* ) yInput );
				else *q_end = q_end_temp;

				if( r_end_temp >  strlen( ( char* ) xInput ) )
					*r_end = strlen( ( char* ) xInput );
				else *r_end = r_end_temp;
					
				*edit_distance = edit_distance_updated;

				return 0;
			}
			
		}
		else if( edit_distance_temp +  edit_distance_L + edit_distance_R <= sw . k ) //extend both directions
		{ 
			q_start_temp--;
			toAddStartQuery++;			

			r_start_temp--;
			toAddStartRef++;

			q_end_temp++;
			toAddEndQuery++;
				
			r_end_temp++;
			toAddEndRef++;

			edit_distance_total_L = edit_distance_L;
			edit_distance_total_R = edit_distance_R;

			edit_distance_updated =  edit_distance_temp + edit_distance_L + edit_distance_R;
		}
	}

	*q_start = 0;
	*r_start = 0;
	*q_end =  strlen( ( char* ) yInput );
	*r_end = strlen( ( char* ) xInput );
	*edit_distance = edit_distance_updated;
	
return 0;
}


int adjust( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw )
{

	int rS = *r_start;
	int qS = *q_start;
	int rE = *r_end;
	int qE = *q_end;

	unsigned char * A = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
	unsigned char * B = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

	memcpy( &A[0], &xInput[rS],  rE - rS  );
	memcpy( &B[0], &yInput[qS],  qE - qS);
	A[ rE - rS ] = '\0';
	B[ qE- qS ] = '\0';
		
        *edit_distance = editDistanceMyers( A, B );

	free( A );
	free( B );

	while( *edit_distance < sw . k )
	{
		if( qS == 0 && rS == 0 && qE ==  strlen( ( char* ) yInput ) && rE == strlen( ( char* ) xInput ) )
			break;

		int eD = *edit_distance;
		
		extend( ( unsigned int*) &eD, (int*) &qS, (int*) &qE, (int*) &rS, (int*) &rE,  xInput, yInput, sw );

		*q_start = qS;
		*q_end = qE;
		*r_start = rS;
		*r_end = rE;

		unsigned char * A2 = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
		unsigned char * B2 = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

		memcpy( &A2[0], &xInput[rS],  rE- rS  );
		memcpy( &B2[0], &yInput[qS],  qE - qS );
		A2[ rE - rS ] = '\0';
		B2[ qE- qS ] = '\0';
		
		*edit_distance =  editDistanceMyers( A2, B2 );

		free( A2 );
		free( B2 );
	}

return 0;
}

template <typename TStringSet, typename TIndexSpec>
int q_gram_counting(TStringSet &set, TIndexSpec )
{

	typedef String<char> TString;
	typedef Index<TStringSet, TIndexSpec> TIndex;
	typedef typename Fibre<TIndex, QGramCounts>::Type TCounts;
	typedef typename Fibre<TIndex, QGramCountsDir>::Type TCountsDir;
	typedef typename Value<TCountsDir>::Type TDirValue;
	typedef typename Iterator<TCounts, Standard>::Type TIterCounts;
	typedef typename Iterator<TCountsDir, Standard>::Type TIterCountsDir;

	TIndex index(set);
	indexRequire(index, QGramCounts());

	int seqNum = countSequences(index);
	Matrix<int, 2> distMat;
	setLength(distMat, 0, seqNum);
	setLength(distMat, 1, seqNum);
	resize(distMat, 0);

	int num_q_grams = 0;

	TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
	TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
	TIterCounts itCountsBegin = begin(indexCounts(index), Standard());

	TDirValue bucketBegin = *itCountsDir;
	for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
	{
		TDirValue bucketEnd = *itCountsDir;

		if (bucketBegin != bucketEnd)
		{
			TIterCounts itA = itCountsBegin + bucketBegin;
			TIterCounts itEnd = itCountsBegin + bucketEnd;
			for(; itA != itEnd; ++itA)
				for(TIterCounts itB = itA; itB != itEnd; ++itB)
					distMat((*itA).i1, (*itB).i1)  += _min((*itA).i2, (*itB).i2);
		}
		bucketBegin = bucketEnd;
	}

return distMat(0,1);
}


int compute_qgrams( unsigned char * m_ref, unsigned char * m_query )
{

	typedef String<char> TString;

   	 TString r = m_ref;
   	 TString q = m_query;
	StringSet<DnaString> stringSet;
	reserve(stringSet, 2); //2 is number of sequences

	appendValue(stringSet, r);
	appendValue(stringSet, q);

	int no_q_grams = q_gram_counting(stringSet, IndexQGram<UngappedShape<4>, OpenAddressing>() );

return no_q_grams;
}

/*
Myers Bit-Vector algorithm implemented using SeqAn Library
www.seqan.de
*/
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput )
{
	typedef String<char> TSequence;

	TSequence seq1 = xInput;
	TSequence seq2 = yInput;

	int score = globalAlignmentScore( seq1, seq2, MyersBitVector() )/-1;

	return score;
}
