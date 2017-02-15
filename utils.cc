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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>
#include "mim.h"

static struct option long_options[] =
 {
   { "alphabet",                	required_argument, NULL, 'a' },
   { "seqs-file",              		required_argument, NULL, 'i' },
   { "output-file",             	required_argument, NULL, 'o' },
   { "min-seq-length",          	required_argument, NULL, 'l' },
   { "max-error-size",          	required_argument, NULL, 'k' },
   { "threads", 			optional_argument, NULL, 'T' },
   { "longest-inc-matches", 	        optional_argument, NULL, 'M' },
   { "rev-compliment",                  optional_argument, NULL, 'r' },
   { "min-cluster-size",                optional_argument, NULL, 'c' },
   { "help",                    	no_argument,       NULL, 'h' },
   { NULL,                      	0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> l                              = 10;
   sw -> k				= 1;
   sw -> M				= 0;
   sw -> c				= 5;
   sw -> T                              = 1;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:l:k:T:r:M:c:h", long_options, &oi ) ) != -1 ) 
    {

      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
          break;

          case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           break;

	case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           break;

	case 'M':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> M = val;
           break;

	case 'c':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> c = val;
           break;

	case 'r':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> r = val;
           break;
	
	 case 'T':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> T = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 3 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }

/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: MIM <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet                 <str>		'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.\n" );
   fprintf ( stdout, "  -i, --input-file               <str>		MultiFASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file              <str>		Output filename with maximal inexact matches.\n" );    
   fprintf ( stdout, "  -l, --min-seq-length           <int>		Minimum length of match.\n" );   
   fprintf ( stdout, "  -k, --max-error-size           <int>		Maximum error size between matches.\n" );
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -M, --longest-inc-matches      <int>		Choose 1 to return all longest increasing maximal inexact matches\n"
                     "                                                or 0 to return all maximal inexact matches. Default: 0.\n" );
   fprintf ( stdout, "  -c, --min-cluster-size         <int>		Minimum number of MIM in each cluster when M=1. Default: 5.\n");
   fprintf ( stdout, "  -r, --rev-compliment           <int>		Choose 1 to compute reverse compliment matches and 0 otherwise. Default: 0.\n");
   fprintf ( stdout, " Number of threads:\n" ); 
   fprintf ( stdout, "  -T, --threads                  <int>		Number of threads to use. Default: 1. \n" );
 }

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

