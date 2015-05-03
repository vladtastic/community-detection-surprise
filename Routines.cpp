#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Routines.h"


int init_expedge( vector<struct exp>* expedge, vector< struct adj>* adjstack, int nv, int ne, int nc ){
	
	if( expedge == NULL ){

		printf("Error: expedge pointer is NULL!\n");
		printf("Source: init_expedge \n");
		return 1;

	}


	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL!\n");
		printf("Source: init_expedge \n");
		return 1;

	}

	
	if( nv < 1 ){

		printf("Error: Invalid number of vertices!\n");
		printf("Source: init_expedge \n");
		return 1;

	}


	if( nc < 1 ){

		printf("Error: Invalid number of clusters! \n");
		printf("Source: init_expedge \n");
		return 1;

	}



	struct exp tempexp;
	double sum = 0.0;

	for( int i = 0; i < nv; i++){
		
		for( int j = 0; j < (int) adjstack->size(); j++){

			if( adjstack->at(j).key1 == i ){

				sum += adjstack->at(j).value;

			}

		}

		tempexp.key1 = i;
		tempexp.value = sum;
		
		expedge->push_back( tempexp );

		sum = 0.0;

	}

	return 0;

}


int locate( vector< struct adj >* adjstack, int keya, int keyb, int index ){


	if( adjstack == NULL ){

		printf("Error: adjstack pointer is null!\n");
		printf("Source: locate \n");
		return 1;

	}


	if( keya < 0 ){

		printf("Error passing key a!\n");
		printf("Source: locate \n");
		return 1;

	}


	if( keyb < 0 ){

		printf("Error passing key b!\n");
		printf("Source: locate \n");
		return 1;

	}


	if( index < 0 ){

		printf("Error passing index!\n");
		printf("Source: locate \n");
		return 1;

	}


	int pos = -1;

	for( int i = index; i < (int) adjstack->size(); i++ ){
			
		if( adjstack->at(i).key1 == keya ){

			if( adjstack->at(i).key2 == keyb ){

				pos = i;

			}

		}

	}

	return pos;

}


int delete_entry( vector< struct adj >* adjstack, int pos ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is null! \n");
		printf("Source: delete \n");
		return 1;
	
	}

	if( pos < 0 ){

		printf("Error passing pos!\n");
		printf("Source: delete \n");
		return 1;

	}


	struct adj tempadj;

	int stack_length = adjstack->size();

	tempadj.key1 = adjstack->at( stack_length - 1 ).key1;
	tempadj.key2 = adjstack->at( stack_length - 1 ).key2;
	tempadj.value = adjstack->at( stack_length - 1 ).value;

	adjstack->at( stack_length - 1).key1 = adjstack->at( pos ).key1;
	adjstack->at( stack_length - 1).key2 = adjstack->at( pos ).key2;
	adjstack->at( stack_length - 1).value = adjstack->at( pos ).value;

	adjstack->at( pos ).key1 = tempadj.key1;
	adjstack->at( pos ).key2 = tempadj.key2;
	adjstack->at( pos ).value = tempadj.value;

	adjstack->pop_back();
			
	return 0;

}



int aggregate_duplicates( vector< struct adj >* adjstack ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is null!\n");
		printf("Source: aggregate_duplicates\n");
		return 1;

	}


	stable_sort( adjstack->begin(), adjstack->end(), UDLess );

	for( unsigned int i = 0; i < adjstack->size(); i++ ){

		int location = locate( adjstack, adjstack->at(i).key1, adjstack->at(i).key2, i+1 );

		while( location != -1 ){

			// A duplicate exists
			// Add the value at location to the adj at i
			// Delete the entry at location
			// Check for another duplicate

			adjstack->at(i).value = adjstack->at(i).value + adjstack->at(location).value;
			adjstack->erase( adjstack->begin() + location );
			location = locate( adjstack, adjstack->at(i).key1, adjstack->at(i).key2, i+1 );


		}

	}

	return 0;

}


int refresh_expedge( vector< struct exp >* expedge, vector< struct adj >* adjstack, int* nc ){

	if( expedge == NULL ){

		printf("Error: expedge pointer is NULL!\n");
		printf("Source: refresh_expedge \n");
		return 1;

	}


	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL!\n");
		printf("Source: refresh_expedge \n");
		return 1;

	}


	if( *nc < 1 ){

		printf("Error: number of clusters is less than one!\n");
		printf("Source: refresh_expedge \n");
		return 1;

	}



	int stack_size = adjstack->size();

	expedge->clear();

	struct exp tempexp;
	double sum = 0.0;

	for( int i = 0; i < *nc; i++ ){

		for( int j = 0; j < stack_size; j++ ){

			if( adjstack->at(j).key1 == i ){

				sum += adjstack->at(j).value;

			}

		}

		tempexp.key1 = i;
		tempexp.value = sum;
		expedge->push_back( tempexp );
		sum = 0.0;

	}


	return 0;

}


int merge_subcluster_supercluster( vector< struct adj >* adjstack, int keya, int keyb ){

	if( adjstack == NULL ){
	
		printf("Error: Null pointer to adjstack!\n");
		printf("Source: merge_subcluster_supercluster \n");
		return 1;

	}

	if( keya < 0 ){

		printf("Error: Invalid value for key a \n");
		printf("Source: merge_subcluster_supercluster \n");
		return 1;

	}

	if( keyb < 0 ){

		printf("Error: Invalid value for key b \n");
		printf("Source: merge_subcluster_supercluster \n");
		return 1;

	}


	double scratch = 0.0;

	int super_cluster_pos = locate( adjstack, keya, keyb, 0 );
	int sub_cluster_pos = locate( adjstack, keyb, keya, 0 );

	scratch = adjstack->at( sub_cluster_pos ).value;

	adjstack->at( super_cluster_pos ).value += scratch;

	scratch = 0.0;

	adjstack->erase( adjstack->begin() + sub_cluster_pos );

	return 0;

}


int join_subcluster_entries_to_supercluster( vector< struct adj >* adjstack, int keya, int keyb ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL! \n ");
		printf("Source: join_subcluster_entries_to_supercluster \n");
		return 1;

	}

	if( keya < 0 ){

		printf("Error: Invalid value of key a!\n");
		printf("Source: join_subcluster_entries_to_supercluster \n");
		return 1;
	
	}

	if( keyb < 0 ){

		printf("Error: Invalid value of key b!\n");
		printf("Source: join_subcluster_entries_to_supercluster \n");
		return 1;

	}
			

	for( int i = 0; i < (int) adjstack->size(); i++ ){

		if( adjstack->at(i).key1 == keyb ){

			adjstack->at(i).key1 = keya;

		}

	 }


	return 0;

}


int relabel_incident_subcluster( vector<struct adj>* adjstack, int keya, int keyb ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL !\n");
		printf("Source: relabel_incident_subcluster \n");
		return 1;

	}


	if( keya < 0 ){

		printf("Error: Invalid value for key a!\n");
		printf("Source: relabel_incident_subcluster\n");
		return 1;

	}

	if( keyb < 0 ){

		printf("Error: Invalid value for key b!\n");
		printf("Source: relabel_incident_subcluster\n");
		return 1;

	}


		
		
	for( int i = 0; i < (int) adjstack->size(); i++ ){

		if( adjstack->at(i).key2 == keyb ){

			adjstack->at(i).key2 = keya;

		}

	}

	return 0;

}


int update_join_list( vector<struct adj>* joinlist , int keya, int keyb ){

	if( joinlist == NULL ){

		printf("Error: joinlist pointer is NULL!\n");
		printf("Source: update_join_list \n");
		return 1;

	}


	if( keya < 0 ){

		printf("Invalid value for key a!\n");
		printf("Source: update_join_list \n");
		return 1;

	}


	if( keyb < 0 ){

		printf("Invalid value for key b!\n");
		printf("Source; update_join_list \n");
		return 1;

	}


	struct adj tempadj;

	tempadj.key1 = keya;
	tempadj.key2 = keya;
	tempadj.value = keyb;

	joinlist->push_back( tempadj );

	return 0;

}


int merge_clusters( vector< struct adj >* adjstack, vector< struct adj>* joinlist, vector<struct exp>* expedge, vector<struct cnfo>* cluster_info, int** edgelist, int keya, int keyb, int* nc, int* nj, int* ne, int* nv ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	if( joinlist == NULL ){

		printf("Error: joinlist pointer is NULL!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	if( cluster_info == NULL ){

		printf("Error: cluster_info is NULL!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	if( keya < 0 ){

		printf("Error: invalid value for key a!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	if( keyb < 0 ){

		printf("Error: invalid value for key b!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	if( nc < 0 ){

		printf("Error: invalid value for number of clusters!\n");
		printf("Source: merge_clusters \n");
		return 1;

	}


	int err = 0;


	// PART ZERO

		err = update_cluster_info( adjstack, cluster_info, edgelist, keya, keyb );

		if( err != 0){

			printf("Error updating cluster info list !\n");
			printf("Source: merge_clusters \n");
			return 1;

		}


	// PART ONE


		err =  merge_subcluster_supercluster( adjstack, keya, keyb );

		if( err != 0 ){

			printf("Error merging subcluster and supercluster \n");
			printf("Source: merge_clusters \n");
			return 1;

		}



	// PART TWO

		err =  join_subcluster_entries_to_supercluster( adjstack, keya, keyb );


		if( err != 0 ){

			printf("Error joining subcluster entries to supercluster \n");
			printf("Source: merge_clusters \n");
			return 1;

		}


	// PART THREE

		err = relabel_incident_subcluster( adjstack, keya, keyb );

		if( err != 0 ){

			printf("Error relabeling incident subclusters \n");
			printf("Source: merge_clusters \n");
			return 1;

		}


	// PART FOUR

		err = aggregate_duplicates( adjstack );
		if( err != 0 ){

			printf("Error aggregating duplicates \n");
			printf("Source: merge_clusters \n");
			return 1;

		}


	// PART FIVE

		*nj = *nj + 1;
		*nc = *nc - 1;

	// PART SIX

		err = update_join_list( joinlist, keya, keyb );

		if( err != 0 ){

			printf("Error updating join list!\n");
			printf("Source: merge_clusters \n");
			return 1;

		}

	// PART SEVEN

		err = refresh_expedge( expedge, adjstack, nc );

		if( err != 0 ){

			printf("Error refreshing expedge!\n");
			printf("Source: merge_clusters \n");
			return 1;

		}

	// PART EIGHT

		double surprise_score = calculate_S( cluster_info, *ne, *nv );

		printf("\tS: %lf", surprise_score );


	return 0;

}

int write_join_list( vector< struct adj>* joinlist, char* filename, int* nj ){

	if( joinlist == NULL ){

		printf("Error: joinlist pointer is NULL!\n");
		printf("Source: write_join_list \n");
		return 1;

	}

	if( filename == NULL ){

		printf("Error: Filename to write joinlist to is NULL \n");
		printf("Source: write_join_list \n");
		return 1;

	}


	if( nj == NULL ){

		printf("Error: number of Joins pointer is NULL!\n");
		printf("Source: write_join_list \n");
		return 1;

	}



	FILE *fp;

	fp = fopen( filename, "w" );

	if( fp == NULL ){

		printf("Error: Unable to open file for join list writing!\n");
		printf("Source: write_join_list \n");
		return 1;

	}



	fprintf(fp, "Logging Checkpoint after %d Joins\n", *nj );

	for( int i = 0; i < (int) joinlist->size(); i++){

		fprintf(fp, "Cluster %d & Cluster %d merged to Cluster %d\n", joinlist->at(i).key2, (int) joinlist->at(i).value, joinlist->at(i).key1 );

	}

	fclose(fp);


	return 0;

}


void print_adjstack( vector<struct adj>* adjstack ){

		int size = adjstack->size();

		for( int i = 0; i < size; i++ ){

			printf(" %d \t %d \t %lf\n", adjstack->at(i).key1, adjstack->at(i).key2, adjstack->at(i).value );

		}

}


void print_expedge( vector<struct exp>* expedge ){

	int size = expedge->size();

	for( int i = 0; i < size; i++ ){

		printf(" %d\t%lf\t\n", expedge->at(i).key1, expedge->at(i).value );

	}

}

bool UDLess( struct adj a1, struct adj a2 ){

	//	Is a1 < a2?

	if( a1.key1 < a2.key1 ){

		return true;

	}

	if( a1.key1 == a2.key1 ){

		if( a1.key2 < a2.key2 ){

			return true;

		}


		else{

			return false;

		}

	}

	return false;

}


int update_cluster_info( vector< struct adj >* adjstack, vector< struct cnfo >* cluster_info, int** edgelist, int keya, int keyb ){

	if( adjstack == NULL ){

		printf("Error: adjstack pointer is NULL!\n");
		printf("Source: update_cluster_info \n");
		return 1;

	}

	if(  cluster_info == NULL ){

		printf("Error: cluster_info pointer is NULL\n");
		printf("Source: update_cluster_info \n");
		return 1;

	}

	if( keya < 0 ){

		printf("Error: Invalid value for keya!\n");
		printf("Source: update_cluster_info \n");
		return 1;

	}

	if( keyb < 0 ){

		printf("Error: Invalid value for keyb!\n");
		printf("Source: update_cluster_info \n");
		return 1;

	}




	int supercluster_location = locate_cluster( cluster_info, keya );

	if( supercluster_location == -1 ){

		printf("Error: Supercluster could not be found in the cluster information list \n");
		printf("Source: update_cluster_info  \n");
		return 1;

	}


	int subcluster_location = locate_cluster( cluster_info, keyb );

	if( subcluster_location == -1 ){

		printf("Error: Subcluster could not be found in the cluster information list \n");
		printf("Source: update_cluster_info \n");
		return 1;

	}


	//	ADD SUBCLUSTER DATA TO SUPERCLUSTER DATA

		cluster_info->at( supercluster_location ).value1  += cluster_info->at( subcluster_location ).value1;


		// UPDATE INTRACLUSTER EDGES
		// DOES AN ELEMENT OF THE SUPERCLUSTER HAVE AN EDGE TO AN ELEMENT OF THE SUBCLUSTER

		vector<int> tmpverts;
		int intraclusterEdges = 0;


		for( unsigned int i = 0; i < cluster_info->at( supercluster_location ).value2.size(); i++ ){

			tmpverts.push_back( cluster_info->at( supercluster_location ).value2.at(i) );

		}

		for( unsigned int j = 0; j < cluster_info->at( subcluster_location ).value2.size(); j++ ){

			tmpverts.push_back( cluster_info->at( subcluster_location ).value2.at(j) );

		}


		for( unsigned int i = 0; i < ( tmpverts.size() - 1); i++ ){

			for( unsigned int j = i+1; j < tmpverts.size(); j++ ){


				if( edgelist[ tmpverts.at(i) -  1][ tmpverts.at(j) - 1 ] == 1 ){

					intraclusterEdges++;

				}

			}

		}


		cluster_info->at( supercluster_location ).value3 = intraclusterEdges;
		intraclusterEdges = 0;


		// PUSH MEMBER LIST OF SUBCLUSTER ONTO SUPERCLUSTER
		for( unsigned int i = 0; i < cluster_info->at( subcluster_location ).value2.size(); i++ ){

			int vertex_member_id = cluster_info->at( subcluster_location ).value2.at(i);

			cluster_info->at( supercluster_location ).value2.push_back( vertex_member_id );

		}


		cluster_info->erase( cluster_info->begin() + subcluster_location );
		tmpverts.clear();

/*

		for( unsigned int i = 0; i < cluster_info->at( supercluster_location ).value2.size(); i++ ){
			for( unsigned int j = 0; j < cluster_info->at( subcluster_location ).value2.size(); j++ ){

				int flag = locate( adjstack, cluster_info->at( supercluster_location ).value2.at(i), cluster_info->at( subcluster_location ).value2.at(j), 0 );

				printf("Flag is %d for Vertex %d and Vertex %d\n", flag, cluster_info->at( supercluster_location ).value2.at(i), cluster_info->at( subcluster_location ).value2.at(j) );

				if( flag != -1 ){

					cluster_info->at( supercluster_location ).value3 = cluster_info->at( supercluster_location ).value3 + 1;

				}

			}

		}

		cluster_info->at( supercluster_location ).value3 += cluster_info->at( subcluster_location ).value3;

		cluster_info->erase( cluster_info->begin() + subcluster_location );
*/
		return 0;

}

int locate_cluster( vector< struct cnfo >* cluster_info, int keya ){

	if( cluster_info == NULL ){

		printf("Error: Null pointer to cluster info list!\n");
		printf("Source: locate_cluster \n");
		return 1;

	}


	if( keya <  0 ){

		printf("Error: Invalid cluster id!\n");
		printf("Source: locate_cluster \n");
		return 1;

	}


	for( unsigned int i = 0; i < cluster_info->size(); i++ ){

		if( cluster_info->at(i).key1 == keya ){

			return i;

		}

	}

	return -1;

}


void print_cluster_info( vector< struct cnfo >* cluster_info ){

	int csize = cluster_info->size();

	for( int i = 0; i < csize; i++ ){

		int member_list_size = cluster_info->at(i).value2.size();

		printf("\n%d\t", cluster_info->at(i).key1 );
		printf("%d\t{", cluster_info->at(i).value1 );

		for( int j = 0; j < member_list_size; j++ ){

			printf("%d ", cluster_info->at(i).value2.at(j) );

		}

		printf("\b}\t%d\n", cluster_info->at(i).value3 );

	}

}


void print_edge_List( int** edgeList, int dimA, int dimB ){

	for( int i = 0; i < dimA; i++ ){

		for( int j = 0; j < dimB; j++ ){

			printf("%d ", edgeList[i][j]  );

		}

		printf("\n");

	}

}


double calculate_S( vector< struct cnfo >* cluster_info, int ne, int nv ){

	unsigned int F = ( nv * ( nv - 1 ) ) / 2;

	unsigned int n = ne;

	unsigned int M = 0;

	unsigned int p = 0;

	for( unsigned int i = 0; i < cluster_info->size(); i++ ){

		unsigned int temp = cluster_info->at(i).value1;

		M = M + ( ( temp ) * ( temp - 1 ) / 2 );

		p = p + cluster_info->at(i).value3;

	}


//	printf("%d\t%d\n", M, n);
	unsigned int upper_limit = min( M , n );
//	printf("%d\n", upper_limit);


	double surprise = 0.0;


//	printf(" Stats:\n ");
//	printf("\tM:%d\tp:%d\tF-M:%d\tn:%d\n", M, p, F-M, n);


	int i = p;

	do{

		surprise += gsl_ran_hypergeometric_pdf( i, M, ( F - M ), n );
//		surprise += gsl_cdf_hypergeometric_P( i, M, (F-M), n );
//		printf("%lf\n", surprise );

		 p = p + 1;

	}

	while( p < upper_limit );

//	printf("\tcalcS: %lf \n", surprise );

//	surprise = gsl_sf_log( surprise );

	return surprise;

}
