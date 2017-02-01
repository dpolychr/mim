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
#include <algorithm>
#include <omp.h>
#include "mim.h"

using namespace std;

bool ordering(MimOcc a, MimOcc b) 
{ 
	if( a.startQuery == b.startQuery )
	{
		return( a.startRef < b.startRef );
	}
	else return (a.startQuery < b.startQuery); 
}

bool includesAll( vector<MimOcc> * a,  vector<MimOcc> * b )
{
	if( a->size() == b->size() )
		return false;

	if ( a->size() > b->size())
	{
		swap(a, b); 
	}

	bool val = false;
	for(int i = 0; i<a->size(); i++)
	{
		val = false;
		for(int j = 0; j<b->size(); j++ )		
		{
			if( a->at(i).startQuery == b->at(j).startQuery && a->at(i).startRef == b->at(j).startRef )
			{
				val = true;
			}
   				
		}

		if( val == false )
			return false;
		else continue;
	}
	return true;
}

int longest_increasing_matches( vector<vector< MimOcc > *> * lims, vector<MimOcc> * mims )
{

	sort(mims->begin(), mims->end(), ordering);
	
	vector<vector<MimOcc> *> * table = new vector <vector<MimOcc> *>;

	bool used = false;
  	for (int k = 0; k <mims->size(); k++)
 	{
   		vector< vector<MimOcc> *> * tempTable = new vector<vector<MimOcc> *>;
   		for (int j=0; j<table->size(); j++)
   		{
    			if ( table->at(j)->at( table->at(j)->size() - 1).startRef <= mims->at(k).startRef )
    			{
				used = true;
     				vector<MimOcc> * temp = new vector<MimOcc>();
				for(int i=0; i<table->at(j)->size(); i++)
				{
     					temp->push_back(table->at(j)->at(i));
				}
     				temp->push_back(mims->at(k));
     				tempTable->push_back(temp);
    			}
			else used = false;
   		}
   			
		for(int j=0; j<tempTable->size(); j++)
		{
			table->push_back( tempTable->at(j) );
		}

   		vector <MimOcc> * temp = new vector <MimOcc >();

		if( used == false )
		{
   			temp->push_back(mims->at(k));
			table->push_back( temp );
		}
  	}

	vector < MimOcc > * toDelete = new vector < MimOcc >();

  	for (int i = 0; i<table->size(); i++)
  	{
		for( int j=0; j< table->size(); j++ )
		{
			if( table->at(j) == NULL )
				continue;
			
			if( table->at(j) !=NULL && table->at(i) != NULL && includesAll( table->at(i), table->at(j) ) == true )
			{	
				if( table->at(i)-> size() > table->at(j) -> size() )
				{	
					table->at( j ) = NULL;
				}
				else table->at( i ) = NULL;
			}	
		}
	}

	for (int i=0; i<table->size(); i++)
  	{
   		if ( table->at(i) != NULL )
  	 	{
    			lims->push_back( table->at(i) );
   		}
		
  	}

return 0;
}


			

