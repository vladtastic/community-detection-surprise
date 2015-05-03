#ifndef ROUTINES_H
#define ROUTINES_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

struct exp{

	int key1;
	double value;

};


struct adj{

	int key1;
	int key2;
	double value;
};


struct cnfo{

	int key1;
	int value1;
	vector< int > value2;
	int value3;

};

using namespace std;

int init_expedge( vector<struct exp>* expedge, vector< struct adj>* adjstack, int nv, int ne, int nc );
int locate( vector< struct adj >* adjstack, int keya, int keyb, int index );
int delete_entry( vector< struct adj >* adjstack, int pos );
int aggregate_duplicates( vector< struct adj >* adjstack, int keyb, int nc );
int refresh_expedge( vector< struct exp >* expedge, vector< struct adj >* adjstack, int* nc );
int merge_subcluster_supercluster( vector< struct adj >* adjstack, int keya, int keyb );
int join_subcluster_entries_to_supercluster( vector< struct adj >* adjstack, int keya, int keyb );
int relabel_incident_subcluster( vector<struct adj>* adjstack, int keya, int keyb );
int update_join_list( vector<struct adj>* joinlist , int keya, int keyb );
int merge_clusters( vector< struct adj >* adjstack, vector< struct adj>* joinlist, vector< struct exp >* expedge, vector< struct cnfo >* cluster_info, int** edgelist, int keya, int keyb, int* nc, int* nj, int* ne, int* nv );
int write_join_list( vector< struct adj>* joinlist, char* filename, int* nj );
void print_adjstack( vector<struct adj>* adjstack );
bool UDLess( struct adj a1, struct adj a2 );
void print_expedge( vector<struct exp>* expedge );
int locate_cluster( vector< struct cnfo >* cluster_info, int keya );
int update_cluster_info( vector< struct adj >* adjstack, vector< struct cnfo >* cluster_info, int** edgelist, int keya, int keyb );
void print_cluster_info( vector< struct cnfo >* cluster_info );
void print_edge_List( int** edgeList, int dimA, int dimB );
double calculate_S( vector< struct cnfo >* cluster_info, int ne, int nv );

#endif
