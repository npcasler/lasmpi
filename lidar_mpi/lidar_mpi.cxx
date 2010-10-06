/********************************************************************
                   GNU LESSER GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.


  This version of the GNU Lesser General Public License incorporates
the terms and conditions of version 3 of the GNU General Public
License, supplemented by the additional permissions listed below.

  0. Additional Definitions.

  As used herein, "this License" refers to version 3 of the GNU Lesser
General Public License, and the "GNU GPL" refers to version 3 of the GNU
General Public License.

  "The Library" refers to a covered work governed by this License,
other than an Application or a Combined Work as defined below.

  An "Application" is any work that makes use of an interface provided
by the Library, but which is not otherwise based on the Library.
Defining a subclass of a class defined by the Library is deemed a mode
of using an interface provided by the Library.

  A "Combined Work" is a work produced by combining or linking an
Application with the Library.  The particular version of the Library
with which the Combined Work was made is also called the "Linked
Version".

  The "Minimal Corresponding Source" for a Combined Work means the
Corresponding Source for the Combined Work, excluding any source code
for portions of the Combined Work that, considered in isolation, are
based on the Application, and not on the Linked Version.

  The "Corresponding Application Code" for a Combined Work means the
object code and/or source code for the Application, including any data
and utility programs needed for reproducing the Combined Work from the
Application, but excluding the System Libraries of the Combined Work.

  1. Exception to Section 3 of the GNU GPL.

  You may convey a covered work under sections 3 and 4 of this License
without being bound by section 3 of the GNU GPL.

  2. Conveying Modified Versions.

  If you modify a copy of the Library, and, in your modifications, a
facility refers to a function or data to be supplied by an Application
that uses the facility (other than as an argument passed when the
facility is invoked), then you may convey a copy of the modified
version:

   a) under this License, provided that you make a good faith effort to
   ensure that, in the event an Application does not supply the
   function or data, the facility still operates, and performs
   whatever part of its purpose remains meaningful, or

   b) under the GNU GPL, with none of the additional permissions of
   this License applicable to that copy.

  3. Object Code Incorporating Material from Library Header Files.

  The object code form of an Application may incorporate material from
a header file that is part of the Library.  You may convey such object
code under terms of your choice, provided that, if the incorporated
material is not limited to numerical parameters, data structure
layouts and accessors, or small macros, inline functions and templates
(ten or fewer lines in length), you do both of the following:

   a) Give prominent notice with each copy of the object code that the
   Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the object code with a copy of the GNU GPL and this license
   document.

  4. Combined Works.

  You may convey a Combined Work under terms of your choice that,
taken together, effectively do not restrict modification of the
portions of the Library contained in the Combined Work and reverse
engineering for debugging such modifications, if you also do each of
the following:

   a) Give prominent notice with each copy of the Combined Work that
   the Library is used in it and that the Library and its use are
   covered by this License.

   b) Accompany the Combined Work with a copy of the GNU GPL and this license
   document.

   c) For a Combined Work that displays copyright notices during
   execution, include the copyright notice for the Library among
   these notices, as well as a reference directing the user to the
   copies of the GNU GPL and this license document.

   d) Do one of the following:

       0) Convey the Minimal Corresponding Source under the terms of this
       License, and the Corresponding Application Code in a form
       suitable for, and under terms that permit, the user to
       recombine or relink the Application with a modified version of
       the Linked Version to produce a modified Combined Work, in the
       manner specified by section 6 of the GNU GPL for conveying
       Corresponding Source.

       1) Use a suitable shared library mechanism for linking with the
       Library.  A suitable mechanism is one that (a) uses at run time
       a copy of the Library already present on the user's computer
       system, and (b) will operate properly with a modified version
       of the Library that is interface-compatible with the Linked
       Version.

   e) Provide Installation Information, but only if you would otherwise
   be required to provide such information under section 6 of the
   GNU GPL, and only to the extent that such information is
   necessary to install and execute a modified version of the
   Combined Work produced by recombining or relinking the
   Application with a modified version of the Linked Version. (If
   you use option 4d0, the Installation Information must accompany
   the Minimal Corresponding Source and Corresponding Application
   Code. If you use option 4d1, you must provide the Installation
   Information in the manner specified by section 6 of the GNU GPL
   for conveying Corresponding Source.)

  5. Combined Libraries.

  You may place library facilities that are a work based on the
Library side by side in a single library together with other library
facilities that are not Applications and are not covered by this
License, and convey such a combined library under terms of your
choice, if you do both of the following:

   a) Accompany the combined library with a copy of the same work based
   on the Library, uncombined with any other library facilities,
   conveyed under the terms of this License.

   b) Give prominent notice with the combined library that part of it
   is a work based on the Library, and explaining where to find the
   accompanying uncombined form of the same work.

  6. Revised Versions of the GNU Lesser General Public License.

  The Free Software Foundation may publish revised and/or new versions
of the GNU Lesser General Public License from time to time. Such new
versions will be similar in spirit to the present version, but may
differ in detail to address new problems or concerns.

  Each version is given a distinguishing version number. If the
Library as you received it specifies that a certain numbered version
of the GNU Lesser General Public License "or any later version"
applies to it, you have the option of following the terms and
conditions either of that published version or of any later version
published by the Free Software Foundation. If the Library as you
received it does not specify a version number of the GNU Lesser
General Public License, you may choose any version of the GNU Lesser
General Public License ever published by the Free Software Foundation.

  If the Library as you received it specifies that a proxy can decide
whether future versions of the GNU Lesser General Public License shall
apply, that proxy's public statement of acceptance of any version is
permanent authorization for you to choose that version for the
Library.
***********************************************************************/



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