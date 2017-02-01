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

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
	FILE *		ref;
	FILE * 		query;
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name

        unsigned char ** seq    = NULL;         // the sequence(s) in memory
        unsigned char ** seq_id = NULL;         // the sequence(s) id in memory

	char *          alphabet;               // the alphabet
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

		if ( sw . q < 2 )
		{
			fprintf ( stderr, " Error: The q-gram length is too small.\n" );
			return ( 1 );	
		}

		if ( sw . k > sw . l )
		{
			fprintf ( stderr, " Error: The number of errors cannot be larger than the length of the match.\n" );
			return ( 1 );	
		}

                input_filename       = sw . input_filename;
		if ( input_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for input!\n" );
			return ( 1 );
		}
		output_filename         = sw . output_filename;

        }

	double start = gettime();


	
	/* Read the (Multi)FASTA file in memory */
	fprintf ( stderr, " Reading the input file: %s\n", input_filename );
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
        unsigned int num_seqs = 0;           // the total number of sequences considered
	unsigned int max_alloc_seq_id = 0;
	unsigned int max_alloc_seq = 0;
	c = fgetc( in_fd );

	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs >= max_alloc_seq_id )
			{
				seq_id = ( unsigned char ** ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len = 0;
			unsigned int seq_id_len = 0;

			seq_id[ num_seqs ] = NULL;

			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id_len )
				{
					seq_id[ num_seqs ] = ( unsigned char * ) realloc ( seq_id[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len += ALLOC_SIZE;
				}
				seq_id[ num_seqs ][ seq_id_len++ ] = c;
			}
			seq_id[ num_seqs ][ seq_id_len ] = '\0';
			
		}

		if ( num_seqs >= max_alloc_seq )
		{
			seq = ( unsigned char ** ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		unsigned int seq_len = 0;
		unsigned int max_alloc_seq_len = 0;

		seq[ num_seqs ] = NULL;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' || c == ' ' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ num_seqs ][ seq_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
				return ( 1 );
			}

		}

		if( seq_len != 0 )
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len += ALLOC_SIZE;
			}
			seq[ num_seqs ][ seq_len ] = '\0';

			total_length += seq_len;
			num_seqs++;
		}
		
	} while( c != EOF );

	seq[ num_seqs ] = NULL;

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if( sw . k < 0 )
	{
		fprintf( stderr, " Error: The error size cannot be negative!\n");
		return ( 1 );	
	}

	ref = fopen ( "ref.fasta", "w");
	fprintf( ref, ">%s\n%s\n", seq_id[0], seq[0] );
	fclose ( ref );
	query = fopen ( "query.fasta", "w");
	fprintf( query, ">%s\n%s\n", seq_id[1], seq[1] );
	fclose ( query );

	fprintf ( stderr, " Finding all maximal inexact matches \n" );

	vector<QGramOcc> * q_grams = new vector<QGramOcc>;

	find_maximal_exact_matches( sw . q , seq[0], seq[1] , q_grams );

	if( q_grams->size() == 0 )
	{
		fprintf( stderr, " Error: The q-gram size is too large!\n");
		return ( 1 );	
	}

	vector<MimOcc> * mims = new vector<MimOcc>;

	find_maximal_inexact_matches( sw , seq[0], seq[1], q_grams, mims );

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
	}
	else
	{
		vector< vector <MimOcc>* > * lims = new vector< vector<MimOcc> *>();

		longest_increasing_matches( lims, mims );

		for ( int i = 0; i < lims->size(); i++ )
		{
			for(int j=0; j< lims->at(i)->size(); j++ )
			{ 
				if ( lims->at(i)->at(j).endQuery - lims->at(i)->at(j).startQuery >= sw . l || lims->at(i)->at(j).endRef - lims->at(i)->at(j).startRef >= sw . l )
				{
					fprintf( out_fd, "%i%s%i%s%i%s%i%s%i\n", lims->at(i)->at(j).startRef, " ", lims->at(i)->at(j).endRef, " ", lims->at(i)->at(j).startQuery, " ", lims->at(i)->at(j).endQuery," ",  lims->at(i)->at(j).error );
				}
			}	
		fprintf(out_fd, "\n");	
		}
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
		free ( seq[i] );
		free( seq_id[i] );
	}	
	free ( seq );
	free ( seq_id );
        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );
	delete( q_grams );
	delete( mims );

	return ( 0 );
}
