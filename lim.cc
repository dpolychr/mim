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
	if( a.startRef == b.startRef)
	{
		return( a.startQuery < b.startQuery );
	}
	else return (a.startRef < b.startRef ); 
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

	int j = 0;

	bool used = false;
  	for (int k = 0; k <mims->size(); k++)
 	{
   		vector< vector<MimOcc> *> * tempTable = new vector<vector<MimOcc> *>;

   		for (int j=0; j<table->size(); j++)
   		{
    			if ( table->at(j)->at( table->at(j)->size() - 1).startRef < mims->at(k).startRef && table->at(j)->at( table->at(j)->size() - 1).endQuery < mims->at(k).startQuery && table->at(j)->at( table->at(j)->size() - 1).endRef < mims->at(k).startRef )
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

   		delete( tempTable );
		vector <MimOcc> * temp = new vector <MimOcc >();
		
		if( used == false )
		{
   			temp->push_back(mims->at(k));
			table->push_back( temp );
		}
  	}

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

