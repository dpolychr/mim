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

#include <vector>
#include <string>
#define INITIAL_SC		-100000
#define ALLOC_SIZE               104857
#define NA			'N'
#define GAP 			'-'
#define ERR                      24
#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*"   //Proteins alphabet
#define DNA                     "ATGCSWRYKMBVHDN"            //IUPAC alphabet
#define INS			1
#define DEL			1
#define SUB			1
#define MAT			0

#define NUC_SCORING_MATRIX_SIZE 15
#define PRO_SCORING_MATRIX_SIZE 24

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)  
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))
#define nuc_delta(a,b) EDNAFULL_matrix[ EDNA[(int)(a)] ][ EDNA[(int)(b)] ]
#define pro_delta(a,b) EBLOSUM62_matrix[ BLOSUM[(int)(a)] ][ BLOSUM[(int)(b)] ]

using namespace std;

struct TSwitch
 {
   //char               * alphabet;
   char               * ref_filename;
   char               * query_filename;
   char               * output_filename;
   unsigned int         matrix;
   int 			T, M, z, x, n, m, a;
   unsigned int         l, k, c, v;
 };

struct QGramOcc
 {
   unsigned int  occRef;
   unsigned int	 occQuery;
   unsigned int  length;
 };

struct MimOcc
 {
   unsigned int  startRef;
   unsigned int	 endRef;
   unsigned int  startQuery;
   unsigned int	 endQuery;
   unsigned int  error;
 };

struct PrevPos_L
 {
  unsigned int prev_L_ref;
  unsigned int prev_L_query;
 };

struct PrevPos_R
 {
  unsigned int prev_R_ref;
  unsigned int prev_R_query;
 };

typedef int32_t INT;


int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
double gettime ( void );
void usage ( void );
int compute_qgrams( unsigned char * m_ref, unsigned char * m_query );
int adjust_all( unsigned char * xInput, unsigned char * yInput, TSwitch sw , vector<MimOcc> * mims );
int find_maximal_exact_matches( unsigned int l, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams );
int find_all_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mnms );
int extend_all( unsigned int * edit_distance,  int * q_start, int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw, vector<MimOcc> * mims);
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput );
int merge( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims );
int longest_increasing_matches( vector<vector< MimOcc > *> * lims, vector<MimOcc> * mims );
unsigned int rev_compliment( unsigned char * str, unsigned char * str2, int iLen );
int adjust2_all(unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw , vector<MimOcc> * mims );
int adjust( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw  );
int find_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mnms );
int extend( unsigned int * edit_distance,  int * q_start, int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw );

