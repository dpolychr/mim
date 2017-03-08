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
#include "mim.h"

int main(int argc, char **argv)
{
	struct TSwitch  sw;

	FILE *          ref_fd;                	 	// the reference file descriptor
	FILE *          query_fd;              	 	// the query file descriptor
	FILE *          out_fd;                 	// the output file descriptor
	//FILE *		ref;
	//FILE * 		query;
        char *          ref_filename;           	// the reference file name
	char *          query_filename;         	// the query file name
        char *          output_filename;        	// the output file name

        unsigned char ** ref    = NULL;         	// the reference sequence in memory
	unsigned char ** query    = NULL;       	// the query sequence in memory
        unsigned char ** seq_id_ref = NULL;     	// the sequences id in memory
	unsigned char ** seq_id_query = NULL;     	// the sequences id in memory

	char *          alphabet;              		// the alphabet
	unsigned int    i, j;
	unsigned int    q, l;
	unsigned int    total_length = 0;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	omp_set_num_threads( sw.T );

	/* Check the arguments */
        if ( i < 3 )
        {
                usage ();
                return ( 1 );
        }
        else
        {

                if      ( ! strcmp ( "DNA", sw . alphabet ) )   { alphabet = ( char * ) DNA;  sw . matrix = 0; }
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  { alphabet = ( char * ) PROT; sw . matrix = 1; }
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences!\n" );
                        return ( 1 );
                }

		if ( ! strcmp ( "PROT", sw . alphabet ) && sw . v == 1 )
		{ 
			
                        fprintf ( stderr, " Error: Can only work out reverse compliment matches for DNA alphabet!\n" );
                        return ( 1 );
                }

		if ( sw . k >= sw . l )
		{
			fprintf ( stderr, " Error: The number of errors must be smaller than the length of the match.\n" );
			return ( 1 );	
		}

                ref_filename         = sw . ref_filename;
		query_filename       = sw . query_filename;
		if ( ref_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for reference sequence!\n" );
			return ( 1 );
		}
		if ( query_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for query sequence!\n" );
			return ( 1 );
		}
		if ( sw . M != 0 && sw . M != 1 )
		{
			fprintf ( stderr, " Error: Invalid value for M, choose 0 for all MIM or 1 for all longest increasing MIM!\n" );
			return ( 1 );
		}
		if ( sw . c <= 0 )
		{
			fprintf ( stderr, " Error: c must be greater than 0!\n" );
			return ( 1 );
		}
		output_filename         = sw . output_filename;

        }

	
	/* Read the FASTA reference in memory */
	fprintf ( stderr, " Reading the reference file: %s\n", ref_filename );
	if ( ! ( ref_fd = fopen ( ref_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", ref_filename );
		return ( 1 );
	}

	char c;
        unsigned int num_seqs = 0;           // the total number of sequences considered
	unsigned int max_alloc_seq_id = 0;
	unsigned int max_alloc_seq = 0;
	c = fgetc( ref_fd );

	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", ref_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs >= max_alloc_seq_id )
			{
				seq_id_ref = ( unsigned char ** ) realloc ( seq_id_ref,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len = 0;
			unsigned int seq_id_len = 0;

			seq_id_ref[ num_seqs ] = NULL;

			while ( ( c = fgetc( ref_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id_len )
				{
					seq_id_ref[ num_seqs ] = ( unsigned char * ) realloc ( seq_id_ref[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len += ALLOC_SIZE;
				}
				seq_id_ref[ num_seqs ][ seq_id_len++ ] = c;
			}
			seq_id_ref[ num_seqs ][ seq_id_len ] = '\0';
			
		}

		if ( num_seqs >= max_alloc_seq )
		{
			ref = ( unsigned char ** ) realloc ( ref,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		unsigned int seq_len = 0;
		unsigned int max_alloc_seq_len = 0;

		ref[ num_seqs ] = NULL;

		while ( ( c = fgetc( ref_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", ref_filename );
				c = fgetc( ref_fd );
				break;
			}
			if( c == '\n' || c == ' ' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq_len )
			{
				ref[ num_seqs ] = ( unsigned char * ) realloc ( ref[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				ref[ num_seqs ][ seq_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", ref_filename, c );
				return ( 1 );
			}

		}

		if( seq_len != 0 )
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				ref[ num_seqs ] = ( unsigned char * ) realloc ( ref[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len += ALLOC_SIZE;
			}
			ref[ num_seqs ][ seq_len ] = '\0';

			total_length += seq_len;
			num_seqs++;
		}
		
	} while( c != EOF );

	if ( num_seqs > 1 )
	{
        	fprintf( stderr, " Warning: %d sequences were read from file %s. Only the first sequence will be processed\n", num_seqs, ref_filename );
	}

	ref[ num_seqs ] = NULL;

	if ( fclose ( ref_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	/* Complete reading reference */


		
	/* Read the FASTA query in memory */
	fprintf ( stderr, " Reading the query file: %s\n", query_filename );
	if ( ! ( query_fd = fopen ( query_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", query_filename );
		return ( 1 );
	}

	char cq;
        unsigned int num_seqs_q = 0;           // the total number of sequences considered
	unsigned int max_alloc_seq_id_q = 0;
	unsigned int max_alloc_seq_q = 0;
	cq = fgetc( query_fd );

	do
	{
		if ( cq != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", query_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs_q >= max_alloc_seq_id_q )
			{
				seq_id_query = ( unsigned char ** ) realloc ( seq_id_query,   ( max_alloc_seq_id_q + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id_q += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len_q = 0;
			unsigned int seq_id_len_q = 0;

			seq_id_query[ num_seqs_q ] = NULL;

			while ( ( cq = fgetc( query_fd ) ) != EOF && cq != '\n' )
			{
				if ( seq_id_len_q >= max_alloc_seq_id_len_q )
				{
					seq_id_query[ num_seqs_q ] = ( unsigned char * ) realloc ( seq_id_query[ num_seqs_q ],   ( max_alloc_seq_id_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len_q += ALLOC_SIZE;
				}
				seq_id_query[ num_seqs_q ][ seq_id_len_q++ ] = cq;
			}
			seq_id_query[ num_seqs_q ][ seq_id_len_q ] = '\0';
			
		}

		if ( num_seqs_q >= max_alloc_seq_q )
		{
			query = ( unsigned char ** ) realloc ( query,   ( max_alloc_seq_q + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq_q += ALLOC_SIZE;
		}

		unsigned int seq_len_q = 0;
		unsigned int max_alloc_seq_len_q = 0;

		query[ num_seqs_q ] = NULL;

		while ( ( cq = fgetc( query_fd ) ) != EOF && cq != '>' )
		{
			if( seq_len_q == 0 && cq == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", query_filename );
				cq = fgetc( query_fd );
				break;
			}
			if( cq == '\n' || cq == ' ' ) continue;

			cq = toupper( cq );

			if ( seq_len_q >= max_alloc_seq_len_q )
			{
				query[ num_seqs_q ] = ( unsigned char * ) realloc ( query[ num_seqs_q ],   ( max_alloc_seq_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len_q += ALLOC_SIZE;
			}

			if( strchr ( alphabet, cq ) )
			{
				query[ num_seqs_q ][ seq_len_q++ ] = cq;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", query_filename, cq );
				return ( 1 );
			}

		}

		if( seq_len_q != 0 )
		{
			if ( seq_len_q >= max_alloc_seq_len_q )
			{
				query[ num_seqs_q ] = ( unsigned char * ) realloc ( query[ num_seqs_q ],   ( max_alloc_seq_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len_q += ALLOC_SIZE;
			}
			query[ num_seqs_q ][ seq_len_q ] = '\0';

			total_length += seq_len_q;
			num_seqs_q++;
		}
		
	} while( cq != EOF );

	if ( num_seqs_q > 1 )
	{
        	fprintf( stderr, " Warning: %d sequences were read from file %s. Only the first sequence will be processed\n", num_seqs_q, query_filename );
	}

	query[ num_seqs_q ] = NULL;

	if ( fclose ( query_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	/* Complete reading query */

	double start = gettime();

	if( sw . k < 0 )
	{
		fprintf( stderr, " Error: The error size cannot be negative!\n");
		return ( 1 );	
	}

	unsigned int q_gram_size = sw . l / ( sw . k + 1 );

	fprintf ( stderr, " Finding all maximal inexact matches \n" );

	vector<QGramOcc> * q_grams = new vector<QGramOcc>;
	vector<MimOcc> * mims = new vector<MimOcc>;

	if( sw . v == 1 )
	{
		unsigned char * rc_seq = ( unsigned char * ) calloc ( ( strlen( ( char* ) query[1] ) + 1 ) , sizeof( unsigned char ) );
			
		rev_compliment( query[0], rc_seq , strlen( ( char* ) query[0] ) - 1 );
		rc_seq[  strlen( ( char* ) query[0] ) ] = '\0';

		find_maximal_exact_matches( q_gram_size , ref[0], rc_seq , q_grams, ref_filename, query_filename );
		find_maximal_inexact_matches( sw , ref[0], rc_seq, q_grams, mims );

		free( rc_seq );


	}
	else 
	{
		find_maximal_exact_matches( q_gram_size , ref[0], query[0] , q_grams,  ref_filename, query_filename );

		find_maximal_inexact_matches( sw , ref[0], query[0], q_grams, mims );
	}

	delete( q_grams );

	fprintf ( stderr, " Preparing the output\n" );

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	if( sw . M == 0 )
	{
		for ( int i = 0; i < mims->size(); i++ )
		{
			if ( mims->at(i).endQuery - mims->at(i).startQuery >= sw . l || mims->at(i).endRef - mims->at(i).startRef >= sw . l )
			{
				fprintf( out_fd, "%i%s%i%s%i%s%i%s%i\n", mims->at(i).startRef, " ", mims->at(i).endRef, " ", mims->at(i).startQuery, " ", mims->at(i).endQuery," ",  mims->at(i).error );
			}		
		}

		delete( mims );
	}
	else
	{
		vector< vector <MimOcc>* > * lims = new vector< vector<MimOcc> *>();

		vector<MimOcc> * mims_del = new vector<MimOcc>;

		for(int i=0; i<mims->size(); i++ )
		{
			if ( mims->at(i).endQuery - mims->at(i).startQuery >= sw . l && mims->at(i).endRef - mims->at(i).startRef >= sw . l )
			{	
				mims_del->push_back( mims->at(i) );
			}
		}
		
		delete( mims );

	
		longest_increasing_matches( lims, mims_del );

		for ( int i = 0; i < lims->size(); i++ )
		{
			if( lims->at(i)->size() >= sw . c )
			{
				for(int j=0; j< lims->at(i)->size(); j++ )
				{ 
					fprintf( out_fd, "%i%s%i%s%i%s%i%s%i\n", lims->at(i)->at(j).startRef, " ", lims->at(i)->at(j).endRef, " ", lims->at(i)->at(j).startQuery, " ", lims->at(i)->at(j).endQuery," ",  lims->at(i)->at(j).error );
				}
			fprintf(out_fd, "\n");
			}	
			else continue;	
		}

		delete( mims_del );
	}

		
	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}


	double end = gettime();

        fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs.\n", num_seqs, ( end - start ) );
	
	/* Deallocate */
	
	for ( i = 0; i < num_seqs; i ++ )
	{
		free ( ref[i] );
		free( seq_id_ref[i] );
	}
	for ( i = 0; i < num_seqs_q; i ++ )
	{
		free ( query[i] );
		free( seq_id_query[i] );
	}		
	free ( ref );
	free ( seq_id_ref );
	free ( query );
	free ( seq_id_query );
        free ( sw . ref_filename );
	free ( sw . query_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );
	
return 1;
}

unsigned int rev_compliment( unsigned char * str, unsigned char * str2, int len )
{
	int i=0;
	while ( len >= 0 )
	{
		if ( str[len] == 'A' )
			str2[i++] = 'T';
		else if( str[len] == 'C')
			str2[i++] = 'G';
		else if( str[len] == 'G')
			str2[i++] = 'C';
		else if( str[len] == 'T')
			str2[i++] = 'A';
		else if( str[len] == 'N')
			str2[i++] = 'N';
		len--;
	}
return 1;
}
