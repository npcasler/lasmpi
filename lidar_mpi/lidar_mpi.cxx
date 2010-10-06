/*******************************************************************/
/* LIDAR Processing with MPI - Shared memory work and sequential   */
/* dumping of results at the end of parallel run into un-ordered   */
/* LAS-file. Since LAS file is un-ordered data can be written as it*/
/* arrives from processing nodes. If striped read is available     */
/* each node can read in and maintain a copy of the data in memory */
/*******************************************************************/

#include <liblas/capi/liblas.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

#define PIPE_MSG 0  /* type of message containing a point to 
                       be checked */
#define END_MSG 1  /* type of message indicating no more data will
                      be coming */

#define inner_ring 5.0 /* Points to ignore in IDW2 */
#define outer_ring 10.0 /* Points to include in IDW2 */

int NNodes,  /* number of nodes in computation*/
    N,  /* Lidar points processed */
    Me,  /* my node number */
    ToCheck;  /* current number to check for passing on to next node;
                 stylistically this might be nicer as a local in
		 Node*(), but I have placed it here to dramatize
		 the fact that the globals are NOT shared among
		 the nodes, thus no problem */

double T1,T2;  /* start and finish times */

typedef double* point_data_t;
vector<point_data_t> las_points;

/*******************************/
/* Read all data into a vector */
/*******************************/
void Read_LAS_file(const char* las_name)
{	
	LASReaderH reader = LASReader_Create(las_name);
	LASPointH laspoint;
	double* coords;
	printf("Attempting file read on Node %d\n",Me);
	
	do
	{
		laspoint = LASReader_GetNextPoint(reader);
		if(LASPoint_GetClassification(laspoint)==2)
		{
			coords = new double[3];
			coords[0] = LASPoint_GetX(laspoint);
			coords[1] = LASPoint_GetY(laspoint);
			coords[2] = LASPoint_GetZ(laspoint);
			las_points.push_back(coords);
		}
		//las_points.push_back(LASPoint_Copy(laspoint));
	}while(laspoint != NULL);
	
	LASReader_Destroy(reader);
	printf("Finished reading %d points on Node %d\n",las_points.size(),Me);
}


/******************************************/
/* Release point data buffers allocated   */
/* on the heap                            */
/******************************************/
void Destroy_LAS_Vector()
{
	vector<point_data_t>::iterator pointVecIterator;
	printf("Start destroying Point vector on Node %d\n",Me);
	for(pointVecIterator = las_points.begin(); 
        pointVecIterator != las_points.end();
        pointVecIterator++)
	{
		delete *pointVecIterator;
		//LASPoint_Destroy(*pointVecIterator);
	}
	printf("Finish destroying Point vector on Node %d\n",Me);

}


/***********************************************/
/* Initialize the MPI blocks with arguments    */
/***********************************************/
void Init(int argc, char ** argv)
{
	MPI_Init(&argc,&argv);

	/* Get the name of the LIDAR File to be processed */
	char* las_name = argv[1];	
	/* puts the number of nodes in NNodes */
	MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
	/* puts the node number of this node in Me */
	MPI_Comm_rank(MPI_COMM_WORLD,&Me);

	Read_LAS_file(las_name);
}

/**************************************/
/* Tear down MPI and Data vectors     */
/**************************************/
void Finalize()
{
	Destroy_LAS_Vector();
	MPI_Finalize();
}

/****************************************************/
/* NodeArb processes every element in the vector      */
/* to calculate distance from received element      */
/* and generate a prediction of height using points */
/* withing 5 units                                  */
/****************************************************/
void NodeArb()
{
	int data_size = las_points.size();
	double test_point[3];
	int nodes_in_play = NNodes-2;
	MPI_Status status;  /* see below */ 
	int Error;
	//Receive points till end and calculate distance to N-th other points
	//for(int i=0;i<=data_size-nodes_in_play;i+=nodes_in_play)
	while(1)
	{
		Error = MPI_Recv(test_point,3,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		//Error in RECV break
		if(Error != MPI_SUCCESS) break;
		//End of points break;
		if(status.MPI_TAG == END_MSG) break;
	
		double wt_sum = 0.0;
		double z_wt_sum = 0.0;

		for(int i=0;i<data_size;i++)
		{
			double xdiff = las_points[i][0]-test_point[0];
			double ydiff = las_points[i][1]-test_point[1];

			double dist_sqr = xdiff*xdiff + ydiff*ydiff;

			double dist = sqrt(dist_sqr);
			
			//Accumulate points in a ring
			if(dist>inner_ring && dist<outer_ring)
			{
				double z = las_points[i][2];
				wt_sum += 1.0/dist_sqr;
				z_wt_sum += z*1.0/dist_sqr;
			}
		}

		double zdiff = test_point[2] - (z_wt_sum/wt_sum);

		//Send off result to serializing node
		Error = MPI_Send(&zdiff,1,MPI_DOUBLE,NNodes-1,PIPE_MSG,MPI_COMM_WORLD);
	}
}

/******************************************************/
/* Last Node in the system performs reduction		  */
/* and serialization								  */
/******************************************************/
void NodeN()
{
	int data_size = las_points.size();
	double zdiff;
	MPI_Status Status;  /* see below */
	FILE* f = fopen("Out.txt","w");
	for(int i=0;i<data_size;i++)
	{
		MPI_Recv(&zdiff,1,MPI_DOUBLE,MPI_ANY_SOURCE,PIPE_MSG,MPI_COMM_WORLD,&Status);
		if(fabs(zdiff)>1.0)
		{
			fprintf(f,"Height Disparity: %f\n",zdiff);
		}
	}
	fclose(f);
}

/*****************************************************/
/* Node0 read points sequentially and messages other */
/* nodes to perfrom IDW2 based checks                */
/*****************************************************/
void Node0()
{
	int data_size = las_points.size();
	double test_point[3];

	//MPI_Status status;  /* see below */ 
	int Error;
	//Read points from the vector and send to all other nodes
	int nodes_in_play = NNodes-2;
	for(int i=0;i<data_size;i++)
	{
		int dest = (i%nodes_in_play)+1;
		//printf("Point %d of %d going to %d\n",i+1,data_size,dest);
		Error = MPI_Send(las_points[i],3,MPI_DOUBLE,dest,PIPE_MSG,MPI_COMM_WORLD);
		if(Error != MPI_SUCCESS) break;
	}
	//Done with all points send end-msg
	for(int dest=1;dest<NNodes-1;dest++)
	{
		printf("Stopping %d\n",dest);
		Error = MPI_Send(test_point,3,MPI_DOUBLE,dest,END_MSG,MPI_COMM_WORLD);
	}
}

void NodeJobs()
{
	if(Me==0)
	{
		Node0();
	}
	else if(Me==NNodes-1)
	{
		NodeN();
	}
	else
	{
		NodeArb();
	}

}

int main(int argc, char ** argv) {

	Init(argc,argv);
	NodeJobs();
	Finalize();
}