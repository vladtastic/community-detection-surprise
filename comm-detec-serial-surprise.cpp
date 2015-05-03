#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <sys/time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Routines.h"


int main( int argc, char* argv[] ){



	int nv, ne, nc = 0;
	int numJoins = 0;
	int k;
	int cluster_candidate_a, cluster_candidate_b = -1;
	int temp1, temp2 = 0;
	int err = 0;
	struct timeval t0, t1;
	double read_io_time, cluster_time, write_io_time = 0.0;
	double Q = 0.0;

	double deltaQ, maxDeltaQ = 0.0;

	char* JOIN_OUT = "joinlist.txt";
	char* TIMING_OUT = "timing.txt";

	vector<struct adj>* adjstack = (vector<struct adj>*) malloc(sizeof(vector<struct adj>));
	vector<struct exp>* expedge = (vector<struct exp>*) malloc(sizeof(vector<struct exp>));
	vector<struct adj>* joinlist = (vector<struct adj>*) malloc(sizeof(vector<struct adj>));
	vector<struct cnfo>*  cluster_info = (vector<struct cnfo>*) malloc(sizeof(vector<struct cnfo>));

	srand( time( NULL ) );


	if( adjstack == NULL ){

		printf("Error allocating adjstack!\n");
		return 1;

	}

	if( expedge == NULL ){

		printf("Error allocating expedge!\n");
		return 1;

	}


	if( joinlist == NULL ){

		printf("Error allocating joinlist!\n");
		return 1;

	}


	//	CHECK CLI ARGS

	if( argc < 2 ){

		printf(" Usage: exec edgelist.txt \n" );
		return 1;

	}

	if( argv[1] == NULL ){

		printf("Error: Invalid edge file!\n");
		return 1;

	}


	//	READ I/O TIMER START
		gettimeofday( &t0, NULL );

	//	READIN NV, NE

	FILE *fp;
	fp = fopen( argv[1], "r" );

	if( fp == NULL ){

		printf("Error opening edge list!\n");
		return 1;

	}

	err = fscanf( fp, "%d\t%d\n", &nv, &ne );

	if( nv == 0 || ne == 0 ){

		printf("Error: Number of Edges and/or Number of Vertices cannot be zero!\n");
		return 1;

	}


	nc = nv;


	//	INIT ADJSTACK
	//		THIS I/O IS SERIAL


	int temp_key1, temp_key2 = 0;
	double temp_value = (1.0 / ( 2.0 * ne ));
	struct adj tempadj;

	int id = 0;
	struct cnfo temp_cnfo;
	int numAdjacent = 0;


	//	ASSUME EDGE LIST IS SORTED (ASCENDING ORDER)
	for( int i = 0; i < ne; i++ ){


		err = fscanf(fp, " %d %d \n ", &temp_key1, &temp_key2);

		tempadj.key1 = temp_key1;
		tempadj.key2 = temp_key2;
		tempadj.value = temp_value;


		adjstack->push_back( tempadj );


		tempadj.key1 = temp_key2;
		tempadj.key2 = temp_key1;
		tempadj.value = temp_value;

		adjstack->push_back( tempadj );


	}


	fclose(fp);





	fp = fopen( argv[1] , "r" );

	int** edgelist = (int**) malloc( nv * sizeof(int*) );

	for( int i = 0; i < nv; i++ ){

		edgelist[i] = (int*) malloc( nv * sizeof(int) );

	}

	for( int i = 0; i < nv; i++ ){

		for( int j = 0; j < nv; j++ ){

			edgelist[i][j] = 0;

		}

	}


	err = fscanf( fp, " %*d %*d " );
	int vert1, vert2 = 0;

	for( int i = 0; i < ne; i++ ){

		err = fscanf( fp," %d %d ", &vert1, &vert2 );

		edgelist[vert1-1][vert2-1] = 1;
		edgelist[vert2-1][vert1-1] = 1;

	}

	fclose(fp);



	//	INITIALIZE CLUSTER INFO LIST

	fp = fopen( argv[1], "r" );

	vector< int > parsed_clusters;

	// Discard the headers
	err = fscanf(fp, " %*d %*d \n" );
	int clusterA, clusterB = 0;
	int flag = -1;
	for( int i = 0; i < ne; i++ ){

		err = fscanf( fp, " %d  %d \n ", &clusterA, &clusterB );

		for( unsigned int j = 0; j < parsed_clusters.size(); j++ ){

			if( parsed_clusters.at(j) == clusterA ){
				flag = 1;
			}
		}

		switch( flag ){

			case( -1 ):

				temp_cnfo.key1 = clusterA;
				temp_cnfo.value1 = 1;
				temp_cnfo.value2.push_back( temp_cnfo.key1 );
				temp_cnfo.value3 = 0;
				cluster_info->push_back( temp_cnfo );
				temp_cnfo.value2.clear();

				parsed_clusters.push_back( clusterA );

				break;
		}


		flag = -1;


		for( unsigned int j = 0; j < parsed_clusters.size(); j++ ){

			if( parsed_clusters.at(j) == clusterB ){
				flag = 1;
			}
		}

		switch( flag ){

			case( -1 ):

				temp_cnfo.key1 = clusterB;
				temp_cnfo.value1 = 1;
				temp_cnfo.value2.push_back( temp_cnfo.key1 );
				temp_cnfo.value3 = 0;
				cluster_info->push_back( temp_cnfo );
				temp_cnfo.value2.clear();

				parsed_clusters.push_back( clusterB );

				break;
		}

		flag = -1;

	}


	fclose(fp);


	//	READ I/O TIMER END

		gettimeofday( &t1, NULL );
		read_io_time = ( t1.tv_sec - t0.tv_sec ) + ( ( t1.tv_usec - t0.tv_usec ) / 1e6 );


	//	INIT EXPEDGE

	err = init_expedge( expedge, adjstack, nv, ne, nc );



	//	CLUSTER TIMER START

		gettimeofday( &t0, NULL );


	//	MAIN ROUTINE

	for( int i = 0; i < (ne - 1); i++){

		maxDeltaQ = -1;

		if( i < (ne/2) ){

			k = 1;

		}

		else{

			k = 2;

		}

		for( int j = 1; j <= k; j++){

			cluster_candidate_a = adjstack->at( rand() % adjstack->size()  ).key1;

			for( int m = 0; m < (int) adjstack->size(); m++ ){

				if( adjstack->at(m).key1 == cluster_candidate_a ){

					if( adjstack->at(m).key2 != cluster_candidate_a ){

						cluster_candidate_b = adjstack->at(m).key2;

						for( int n = 0; n < (int) expedge->size(); n++ ){

							if( expedge->at(n).key1 == cluster_candidate_a ){

								temp1 = expedge->at(n).value;

							}


							if( expedge->at(n).key1 == cluster_candidate_b ){

								temp2 = expedge->at(n).value;

							}

						}

						deltaQ = 2 * ( adjstack->at(m).value - ( temp1 * temp2 ) ) ;

						if( deltaQ > maxDeltaQ ){

							maxDeltaQ = deltaQ;

						}

					}
				}

			}

			if( cluster_candidate_b != -1 ){

				err = merge_clusters( adjstack, joinlist, expedge, cluster_info, edgelist, cluster_candidate_a, cluster_candidate_b, &nc, &numJoins, &ne, &nv );

				if( err != 0 ){

					printf("Error: Could not join clusters!\n");
					printf("Source: main \n");
					return 1;

				}

			}

			cluster_candidate_a = -1;
			cluster_candidate_b = -1;

		}

		Q += maxDeltaQ;

		printf("\tQ: %lf\n", Q);


	}




	//	CLUSTER TIMER END
		gettimeofday( &t1, NULL );
		cluster_time = ( t1.tv_sec - t0.tv_sec ) + ( ( t1.tv_usec - t0.tv_usec ) / 1e6 );


	printf("Clustering Completed.\n");
	printf("Writing Join List.\n");


	//	WRITE I/O TIMER START
		gettimeofday( &t0, NULL );

	err = write_join_list( joinlist, JOIN_OUT, &numJoins );


	//	WRITE I/O TIMER END

		gettimeofday( &t1, NULL );
		write_io_time = (t1.tv_sec - t0.tv_sec ) + ( ( t1.tv_usec - t0.tv_usec ) / 1e6 );

	if( err != 0 ){

		printf("Error: Unable to write join list!\n");
		printf("Source: main \n");
		return 1;

	}


	printf("Writing Timing File.\n");

	FILE *fp2;
	fp2 = fopen( TIMING_OUT, "w" );

	fprintf( fp2, "Number of Vertices: %d\n", nv );
	fprintf( fp2, "Number of Edges: %d\n\n", ne );
	fprintf( fp2, "Read I/O Time: %lf\n", read_io_time );
	fprintf( fp2, "Cluster Time: %lf\n", cluster_time );
	fprintf( fp2, "Write I/O Time: %lf\n", write_io_time );

	fclose(fp2);

	free( adjstack );
	free( expedge );
	free( joinlist );

	return 0;

}
