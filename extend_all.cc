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
#include <fstream>
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

bool order_all(MimOcc a, MimOcc b) 
{ 
	if( a.startRef == b.startRef)
	{
		return( a.startQuery < b.startQuery );
	}
	else return (a.startRef < b.startRef ); 
}

bool uniqueEnt_all(MimOcc a, MimOcc b) 
{
	if( a.startRef == b.startRef && a.endRef == b.endRef && a.startQuery == b.startQuery && a.endQuery == b.endQuery )
	{
  		return ( true );
	}
	else return ( false );
}

int find_all_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims )
{

	fprintf ( stderr, " -Merging exact matches\n" );
	merge( sw, ref, query, q_grams, mims );

	fprintf ( stderr, " -Extending merged matches\n" );

	int mimSize = mims->size();
	for( int i=0; i<mimSize; i++ )
	{
		if( mims->at(i). error < sw . k )
		{
			extend_all( &mims->at(i).error, (int*) &mims->at(i).startQuery, (int*) &mims->at(i).endQuery, (int*) &mims->at(i).startRef, (int*) &mims->at(i).endRef, ref, query, sw, mims );
		}
	}


	mims->erase(mims->begin(), mims->begin() + mimSize );

	fprintf ( stderr, " -Adjusting extended matches\n" );


	adjust_all( ref, query, sw, mims );

	for( int j=0; j<mims->size(); j++ )
	{
		adjust2_all(  &mims->at(j).error, (int*) &mims->at(j).startQuery, (int*) &mims->at(j).endQuery, (int*) &mims->at(j).startRef, (int*) &mims->at(j).endRef, ref, query, sw, mims );
	}


	sort( mims->begin(), mims->end(), order_all );

	unique( mims->begin(), mims->end(), uniqueEnt_all ); 

return 0;
}

int extend_all( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw, vector<MimOcc> * mims )
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

	char operationEnd;
	char operationStart;

	/************************************************ Score for extending right ***************************************************/
	int edit_distance_R = edit_distance_temp;
	int edit_distance_R_prev = edit_distance_temp;

	vector<PrevPos_R> * rightPositions = new vector<PrevPos_R>;
	
	PrevPos_R firstR;
	firstR.prev_R_ref = *r_end;
	firstR.prev_R_query = *q_end;
		
	rightPositions->push_back( firstR );

	while( sw . k  >= edit_distance_R - edit_distance_temp )
	{ 
	
		if( q_start_temp == 0 && r_start_temp == 0 && q_end_temp ==  strlen( ( char* ) yInput ) && r_end_temp == strlen( ( char* ) xInput ) )
		{
			MimOcc newOcc;
			newOcc.startRef = 0;
			newOcc.startQuery = 0;
		
			newOcc.endQuery = q_end_temp ==  strlen( ( char* ) yInput );
			newOcc.endRef =  r_end_temp == strlen( ( char* ) xInput );
			newOcc.error = *edit_distance;

			mims->push_back( newOcc );


			break;
		}	

		if(   q_end_temp == strlen( ( char* ) yInput ) && r_end_temp == strlen( ( char* ) xInput ) )
		{
			edit_distance_R =  edit_distance_R + sw . k;
		}
		if (  q_end_temp  < strlen( ( char* ) yInput )  &&  r_end_temp  < strlen( ( char* ) xInput ) ) 
		{	

			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
			memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';

			int editDist_S = editDistanceMyers( m_ref_R, m_query_R );
			int editDist_I;
			int editDist_D;

			if( toAddEndRef > 1 )
			{
				memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
				memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
				m_ref_R[ toAddEndRef - 1 ] = '\0';
				m_query_R[ toAddEndQuery ] = '\0';

				editDist_I = editDistanceMyers( m_ref_R, m_query_R );

				memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
				memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
				m_ref_R[ toAddEndRef ] = '\0';
				m_query_R[ toAddEndQuery - 1 ] = '\0';

				editDist_D = editDistanceMyers( m_ref_R, m_query_R );

			}
			else
			{
				editDist_I = sw . k + 1;
				editDist_D = sw . k + 1;
			}

			edit_distance_R_prev = edit_distance_R;
			edit_distance_R =  edit_distance_temp + min( editDist_S, min( editDist_I, editDist_D ) );
			
			if( edit_distance_R ==  edit_distance_temp + editDist_S )
				operationEnd = 'S';
			else if( edit_distance_R ==  edit_distance_temp + editDist_I )
				operationEnd = 'I';
			else if( edit_distance_R == edit_distance_temp + editDist_D )
				operationEnd = 'D';


			free( m_ref_R );
			free( m_query_R );

		}
		else if( qE == strlen( ( char* ) yInput ) && rE != strlen( ( char* ) xInput ) )
		{
			edit_distance_R_prev = edit_distance_R;
			edit_distance_R = edit_distance_R + 1;
			operationEnd = 'D';
		}
		else if( rE == strlen( ( char* ) xInput ) && qE != strlen( ( char* ) yInput ) )
		{
			edit_distance_R_prev = edit_distance_R;
			edit_distance_R = edit_distance_R + 1;
			operationEnd = 'I';
		}
		else if ( q_end_temp  < strlen( ( char* ) yInput ) && r_end_temp >= strlen( ( char* ) xInput ) )	
		{
			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  strlen( ( char* ) xInput ) - rE );
			memcpy( &m_query_R[0], &yInput[qE], toAddEndQuery  );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';
	
			edit_distance_R_prev = edit_distance_R;
			edit_distance_R =  edit_distance_temp + editDistanceMyers( m_ref_R, m_query_R );
			free( m_ref_R );
			free( m_query_R );

			operationEnd = 'I';
		}
		else if ( q_end_temp  >= strlen( ( char* ) yInput ) - 1 && r_end_temp < strlen( ( char* ) xInput ) - 1 )	
		{
			unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef );
			memcpy( &m_query_R[0], &yInput[qE], strlen( ( char* ) yInput ) - qE   );
			m_ref_R[ toAddEndRef ] = '\0';
			m_query_R[ toAddEndQuery ] = '\0';
	
			edit_distance_R_prev = edit_distance_R;
			edit_distance_R =  edit_distance_temp + editDistanceMyers( m_ref_R, m_query_R );
			free( m_ref_R );
			free( m_query_R );

			operationEnd = 'D';
		}

		
		if( operationEnd == 'S' )
		{
			q_end_temp++;
			r_end_temp++;
			toAddEndQuery++;
			toAddEndRef++;
		}
		else if( operationEnd == 'I' )
		{		
			toAddEndQuery++;
			q_end_temp++;

		}
		else if( operationEnd == 'D' )
		{	
			toAddEndRef++;
			r_end_temp++;

		}
	
		if( edit_distance_R <= sw . k )
		{
			if( edit_distance_R < edit_distance_R_prev && rightPositions->size() != 0  )
			{
				while ( rightPositions->size() - 1 > edit_distance_R - edit_distance_temp )
				{
					rightPositions->pop_back();
				}	
			}	
			else if ( edit_distance_R_prev == edit_distance_R && rightPositions->size() != 0 )
				rightPositions->pop_back();


			firstR.prev_R_ref = r_end_temp;
			firstR.prev_R_query = q_end_temp;

			rightPositions->push_back( firstR );
		}
	}


	/*********************************************** score for extending left *************************************************/
	int edit_distance_L = edit_distance_temp;
	int edit_distance_L_prev = edit_distance_temp;

	vector<PrevPos_L>  * leftPositions = new vector<PrevPos_L>;

	PrevPos_L firstL;
	firstL.prev_L_ref = *r_start;
	firstL.prev_L_query = *q_start;
	
	leftPositions->push_back( firstL );

	while( sw . k  >= edit_distance_L - edit_distance_temp )
	{ 
		if( q_start_temp == 0 &&  r_start_temp == 0 )
		{
			edit_distance_L = edit_distance_L + sw . k;
		}
		if(  q_start_temp  > 0 &&  r_start_temp > 0   )  
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );

			memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
			memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';

			int editDist_S = editDistanceMyers( m_ref_L, m_query_L );
			int editDist_I;
			int editDist_D;

			if( toAddStartRef > 1 )
			{
				memcpy( &m_ref_L[0], &xInput [ rS - toAddStartRef + 1  ], toAddStartRef  );
				memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
				m_ref_L[ toAddStartRef - 1 ] = '\0';
				m_query_L[ toAddStartQuery ] = '\0';
				
				editDist_I = editDistanceMyers( m_ref_L, m_query_L );

				memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
				memcpy( &m_query_L[0], &yInput [ qS - toAddStartQuery + 1 ], toAddStartQuery   );
				m_ref_L[ toAddStartRef ] = '\0';
				m_query_L[ toAddStartQuery - 1 ] = '\0';

				editDist_D = editDistanceMyers( m_ref_L, m_query_L );
			

			}
			else
			{
				editDist_I = sw . k + 1;
				editDist_D = sw . k + 1;
			}

			edit_distance_L_prev = edit_distance_L;
			edit_distance_L = edit_distance_temp + min( editDist_S, min( editDist_I, editDist_D ) );

			if( edit_distance_L ==  edit_distance_temp + editDist_S )
				operationStart = 'S';
			else if( edit_distance_L == edit_distance_temp +  editDist_I )
				operationStart = 'I';
			else if( edit_distance_L ==  edit_distance_temp + editDist_D )
				operationStart = 'D';

			free( m_ref_L );

			free( m_query_L );

		}
		else if( qS == 0 && rS != 0 )
		{
			edit_distance_L_prev = edit_distance_L;
			edit_distance_L = edit_distance_L + 1;
			operationStart = 'D';
		}
		else if( rS == 0 && qS != 0 )
		{
			edit_distance_L_prev = edit_distance_L;
			edit_distance_L = edit_distance_L + 1;
			operationStart = 'I';
		}
		else if ( q_start_temp  <= 0 && r_start_temp > 0 )	
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
			
			memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
			memcpy( &m_query_L[0], &yInput [0], qS );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';

			edit_distance_L_prev = edit_distance_L;
			edit_distance_L =  edit_distance_temp + editDistanceMyers( m_ref_L, m_query_L );
			free( m_ref_L );
			free( m_query_L );

			operationStart = 'D';
		}
		else if ( q_start_temp  > 0 && r_start_temp <= 0 )	
		{
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
			
			memcpy( &m_ref_L[0], &xInput [0], rS );
			memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
			m_ref_L[ toAddStartRef ] = '\0';
			m_query_L[ toAddStartQuery ] = '\0';

			edit_distance_L_prev = edit_distance_L;
			edit_distance_L = edit_distance_temp + editDistanceMyers( m_ref_L, m_query_L );
			free( m_ref_L );
			free( m_query_L );

			operationStart = 'I';
		}
		

		if( operationStart == 'S' )
		{
			q_start_temp--;
			r_start_temp--;
			toAddStartQuery++;
			toAddStartRef++;
		}
		else if( operationStart == 'I' )
		{
			toAddStartQuery++;
			q_start_temp--;
		}

		else if( operationStart == 'D' )
		{
			toAddStartRef++;
			r_start_temp--;
		}


		if( edit_distance_L <= sw . k)
		{
			if( edit_distance_L_prev > edit_distance_L && leftPositions->size() != 0 )
			{
				while ( leftPositions->size() >  edit_distance_L - edit_distance_temp )
				{
					leftPositions->erase( leftPositions->begin() );
				}		
			}
			else if ( edit_distance_L_prev == edit_distance_L && leftPositions->size() != 0 )
			{
				leftPositions->erase( leftPositions->begin() );
			}
			firstL.prev_L_ref = r_start_temp;
			firstL.prev_L_query = q_start_temp;

			leftPositions->insert( leftPositions->begin(),firstL );
		}
	}

	/*********************************************** computing extension *************************************************/

	int k =  min( leftPositions->size(), rightPositions->size() ); 
	while(  k > 0 )
	{
		MimOcc newOcc;
		newOcc.startRef = leftPositions->at(0).prev_L_ref;
		newOcc.startQuery = leftPositions->at(0).prev_L_query;

		newOcc.endQuery = rightPositions->at( sw .k - ( leftPositions->size() - 1 + edit_distance_temp) ).prev_R_query;
		newOcc.endRef =  rightPositions->at(sw .k - ( leftPositions->size() - 1 + edit_distance_temp) ).prev_R_ref;
		newOcc.error = sw . k;

		mims->push_back( newOcc );

		leftPositions->erase( leftPositions->begin() );
		k--;
	}

return 0;
}

int adjust_all( unsigned char * xInput, unsigned char * yInput, TSwitch sw, vector<MimOcc> * mims )
{


	vector<MimOcc> * newmims = new vector<MimOcc>;
	int mimsSize = mims->size();
	for(int i = 0; i<mims->size(); i++ )
	{
		int rS = mims->at(i).startRef;
		int qS = mims->at(i).startQuery;
		int rE = mims->at(i).endRef;
		int qE = mims->at(i).endQuery;

		unsigned char * A = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
		unsigned char * B = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

		memcpy( &A[0], &xInput[rS],  rE - rS  );
		memcpy( &B[0], &yInput[qS],  qE - qS);
		A[ rE - rS ] = '\0';
		B[ qE- qS ] = '\0';
		
		
		int editDistance = editDistanceMyers( A, B );
		free( A );
		free( B );
	
		if( editDistance < sw . k )
		{

			int eD = editDistance;

			extend_all( ( unsigned int*) &eD, (int*) &qS, (int*) &qE, (int*) &rS, (int*) &rE,  xInput, yInput, sw, mims );
				
		}
		else newmims->push_back( mims->at( i ) );
	}

	mims->clear();

	for(int i=0; i<newmims->size(); i++)
		mims->push_back( newmims->at(i) ) ;

	delete( newmims );
	

	
return 0;
}

int adjust2_all( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw, vector<MimOcc> * mims )
{

	int rS = *r_start;
	int qS = *q_start;
	int rE = *r_end;
	int qE = *q_end;

	/* can still extend if distance = sw . k */


	while( *edit_distance <= sw . k )
	{
		/* extending right again */

		int editDist_Sr;
		int editDist_Ir;
		int editDist_Dr;
		unsigned char * m_ref_R = ( unsigned char * ) calloc (  rE - rS + 2, sizeof ( unsigned char ) );
		unsigned char * m_query_R = ( unsigned char * ) calloc ( qE - qS + 2, sizeof ( unsigned char ) );

		if( qE <  strlen( ( char* ) yInput ) && rE < strlen( ( char* ) xInput ) )
		{ 
			memcpy( &m_ref_R[0], &xInput[rS],  rE - rS + 1  );
			memcpy( &m_query_R[0], &yInput[qS],  qE - qS + 1  );
			m_ref_R[ rE - rS + 1 ] = '\0';
			m_query_R[ qE - qS + 1 ] = '\0';

			editDist_Sr = editDistanceMyers( m_ref_R, m_query_R );
		}
		else editDist_Sr = sw . k + 2;

		if( qE <  strlen( ( char* ) yInput ) )
		{
			memcpy( &m_ref_R[0], &xInput[rS],  rE - rS   );
			memcpy( &m_query_R[0], &yInput[qS],  qE - qS + 1  );
			m_ref_R[ rE - rS ] = '\0';
			m_query_R[ qE - qS + 1 ] = '\0';

			editDist_Ir = editDistanceMyers( m_ref_R, m_query_R );
		}
		else editDist_Ir = sw . k + 2;

		if( rE <  strlen( ( char* ) yInput ) )
		{
			memcpy( &m_ref_R[0], &xInput[rS],  rE - rS  + 1 );
			memcpy( &m_query_R[0], &yInput[qS],  qE - qS );
			m_ref_R[ rE - rS  + 1] = '\0';
			m_query_R[ qE - qS ] = '\0';

			editDist_Dr = editDistanceMyers( m_ref_R, m_query_R );
		}
		else editDist_Dr = sw . k + 2;

		free(m_ref_R);
		free(m_query_R);

		/* extending left again */
		
		int editDist_Sl;
		int editDist_Il;
		int editDist_Dl;

		unsigned char * m_ref_L = ( unsigned char * ) calloc (  rE - rS + 2, sizeof ( unsigned char ) );
		unsigned char *	m_query_L = ( unsigned char * ) calloc ( qE - qS + 2, sizeof ( unsigned char ) );;

		if( qS > 0 && rS > 0 )
		{
	
			memcpy( &m_ref_L[0], &xInput [ rS - 1  ], rE - rS + 1 );
			memcpy( &m_query_L[0], &yInput[ qS - 1 ], qE - qS + 1 );
			m_ref_L[ rE - rS + 1 ] = '\0';
			m_query_L[ qE - qS + 1] = '\0';

			editDist_Sl = editDistanceMyers( m_ref_L, m_query_L );

		}	
		else editDist_Sl = sw . k + 2;

		if( qS > 0 )
		{
			memcpy( &m_ref_L[0], &xInput [rS ], rE - rS );
			memcpy( &m_query_L[0], &yInput [ qS - 1 ], qE - qS + 1  );
			m_ref_L[ rE - rS  ] = '\0';
			m_query_L[ qE - qS + 1 ] = '\0';

			editDist_Il = editDistanceMyers( m_ref_L, m_query_L );
		}
		else editDist_Il = sw . k + 2;
		
		if( rS > 0 )
		{
			memcpy( &m_ref_L[0], &xInput [rS - 1], rE - rS + 1 );
			memcpy( &m_query_L[0], &yInput [ qS ], qE - qS   );
			m_ref_L[ rE - rS + 1  ] = '\0';
			m_query_L[ qE - qS ] = '\0';

			editDist_Dl = editDistanceMyers( m_ref_L, m_query_L );
		}
		else editDist_Dl = sw . k + 2;

		free(m_ref_L);
		free(m_query_L);

		int left1 = min( editDist_Sl, min( editDist_Il, editDist_Dl ) );
		int right1 = min( editDist_Sr, min( editDist_Ir, editDist_Dr ) );
		int minDist =  min( left1, right1 );
		
		if( minDist == editDist_Sr )
		{
			rE = rE + 1;
			qE = qE + 1;
		}
		else if( minDist == editDist_Ir )
		{
			rE = rE;
			qE = qE + 1;
		}
		else if( minDist == editDist_Dr )
		{
			rE = rE + 1;
			qE = qE;
		}
		else if( minDist == editDist_Sl )
		{
			rS = rS - 1;
			qS = qS - 1;	

		}
		else if( minDist == editDist_Il )
		{
			rS = rS;
			qS = qS - 1;

		}
		else if( minDist == editDist_Dl )
		{
			rS = rS - 1;
			qS = qS;
		}


		if( minDist <= sw . k )
		{
			*edit_distance = minDist;
			*q_start = qS;
			*q_end = qE;
			*r_start = rS;
			*r_end = rE;
		}
		else break;
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
