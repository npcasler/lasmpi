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

/**********************************************************************/
/* This program aims to reclassify short-vegetation using a mask and  */
/* a dilated buffer to mark points in the mask for reclassification   */
/* then re-added to the ground class if they are within a threshold   */
/* of an IDW interpolation they are added back to the ground class    */
/**********************************************************************/

#include <liblas/capi/liblas.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <limits>
#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()


using namespace std;

#define PIPE_MSG 0  /* type of message containing a point to 
                       be checked */
#define END_MSG 1  /* type of message indicating no more data will
                      be coming */

typedef LASPointH point_data_t;
vector<point_data_t> las_points;


char* las_out_file; /* LAS Output file parameter */
LASWriterH las_writer; /* Handle to LASWriter Valid only on NodeN */
 
int NNodes,  /* number of nodes in computation*/
    N,  /* Lidar points processed */
    Me;  /* my node number */

double add_back_threshold;


/* GDAL data buffers */
void* maskArr;
void* bufArr;
int xsize,ysize;
double transform[6];

/****************************************/
/* Initiate GDAL files for reading		*/
/* Depending on architecture 1 or all	*/
/* nodes have access to them, possibly	*/
/* read whole file into a buffer to		*/
/* reduce disk costs					*/
/****************************************/
void Read_GDAL_files(const char* mask_file_name, const char* buffer_file_name)
{
	GDALAllRegister();

	GDALDataset* maskDataset = (GDALDataset*) GDALOpen( mask_file_name, GA_ReadOnly );
	GDALDataset* bufferDataset = (GDALDataset*) GDALOpen( buffer_file_name, GA_ReadOnly );

	if(maskDataset!=NULL && bufferDataset!=NULL)
	{
		xsize = maskDataset->GetRasterXSize();
		ysize = maskDataset->GetRasterYSize();

		maskDataset->GetGeoTransform(transform);

		CPLAssert( maskDataset->GetRasterBand(1)->GetRasterDataType() == GDT_Byte );
		CPLAssert( bufferDataset->GetRasterBand(1)->GetRasterDataType() == GDT_Byte );

		maskArr = CPLMalloc(xsize*ysize);
		bufArr = CPLMalloc(xsize*ysize);

		maskDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, xsize, ysize, 
			maskArr, xsize, ysize, GDT_Byte, 
			0, 0 );

		maskDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, xsize, ysize, 
			bufArr, xsize, ysize, GDT_Byte, 
			0, 0 );

	}
}

/********************************************/
/* Release GDAL buffer allocated during Init*/
/********************************************/

void Destroy_GDAL_Array()
{
	CPLFree((void*)maskArr);
	CPLFree((void*)bufArr);
}

/*******************************/
/* Read all data into a vector */
/*******************************/
void Read_LAS_file(const char* las_name)
{	
	LASReaderH reader = LASReader_Create(las_name);
	LASPointH laspoint;
	
	/* Set up LASWriter on NodeN */
	if (Me==NNodes-1)
	{
		LASHeaderH header = LASReader_GetHeader(reader);
		las_writer = LASWriter_Create(las_out_file,header,LAS_MODE_WRITE);
	}

	do
	{
		laspoint = LASReader_GetNextPoint(reader);
		if(LASPoint_GetClassification(laspoint)==2)
		{
			las_points.push_back(LASPoint_Copy(laspoint));
		}
		else
		{
			if (Me==NNodes-1)
			{
				LASWriter_WritePoint(las_writer,laspoint);
			}
		}

	}while(laspoint != NULL);
	printf("Finished reading %d points on Node %l\n",las_points.size(),Me);
	LASReader_Destroy(reader);
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
		LASPoint_Destroy(*pointVecIterator);
	}
	printf("Finish destroying Point vector on Node %d\n",Me);
	/* Clean up the node containing the LASWriter */
	if(Me==NNodes-1)
	{
		LASWriter_Destroy(las_writer);
	}
}


/***********************************************/
/* Initialize the MPI blocks with arguments    */
/***********************************************/
void Init(int argc, char ** argv)
{
	MPI_Init(&argc,&argv);

	/* Get the name of the LIDAR File to be processed */
	char* las_name = argv[1];
	/* Get the name of the output las file */
	las_out_file = argv[2];
	/* Get the name of the mask raster file */
	char* mask_file_name = argv[3];
	/* Get the name of the buffer raster file */
	char* buffer_file_name = argv[4];
	/* Get some points back */
	add_back_threshold = atof(argv[5]);

	/* puts the number of nodes in NNodes */
	MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
	/* puts the node number of this node in Me */
	MPI_Comm_rank(MPI_COMM_WORLD,&Me);

	/* Buffer up GDAL based mask */
	Read_GDAL_files(mask_file_name,buffer_file_name);

	/* Las data is only needed on intermediate nodes */
	Read_LAS_file(las_name);
}

/**************************************/
/* Tear down MPI and Data vectors     */
/**************************************/
void Finalize()
{
	/* Las data is only needed on intermediate nodes */
	Destroy_LAS_Vector();
	MPI_Finalize();
}


/******************************************/
/* Node0 sends ground points off to other */
/* Nodes for consideration against masks, */
/* IDW estimates etc					  */
/******************************************/
void Node0()
{
	unsigned int data_size = las_points.size();
	int nodes_in_play = NNodes-2;

	//MPI_Status status;  /* see below */ 
	int Error;

	int dummy[2];

	for(unsigned int i = 0; i< data_size; i++)
	{
		/* Should really do something more useful here */
		/* Like look up GDAL array for whether it is in */
		/* Mask or Buffer or outside both but w/e */
		int dest = i%nodes_in_play + 1;
		dummy[0] = i;
		LASPointH test_point = las_points[i];
		
		double xloc = LASPoint_GetX(test_point);
		double yloc = LASPoint_GetY(test_point);

		int xoff = (int)ceil((xloc-transform[0])/transform[1]);
		int yoff = (int)ceil((yloc-transform[3])/transform[5]);

		/******************************************/
		/* What to do with LASPoint			      */
		/* 0 - pass it on to written		      */
		/* 1 - we are in mask do something !!     */
		/* 2 - we are in buffer save me and write */
		/******************************************/
		int decision = 0;
		/* Ensure xoff/yoff is in range */
		if(xoff<xsize && yoff <ysize && xoff>0 && yoff>0)
		{
			GByte mask = ((GByte*)maskArr)[xoff+yoff*xsize];
			GByte buffer = ((GByte*)bufArr)[xoff+yoff*xsize];

			if(mask == 255) decision = 1;
			else if(buffer == 255) decision = 2;
		}

		dummy[2] = decision;

		Error = MPI_Send(dummy,2,MPI_INTEGER,dest,PIPE_MSG,MPI_COMM_WORLD);
	}

	//Done with all points send end-msg
	for(int dest=1;dest<NNodes-1;dest++)
	{
		printf("Stopping %d\n",dest);
		Error = MPI_Send(dummy,2,MPI_INTEGER,dest,END_MSG,MPI_COMM_WORLD);
	}
}

/******************************************/
/* NodeArb Receives LASPoints and decides */
/* Whether to send them off to be written */
/* or to reclassify them after mask+IDW	  */
/* consideration						  */
/******************************************/
void NodeArb()
{
	int data_size = las_points.size();
	int nodes_in_play = NNodes-2;
	MPI_Status status;  /* see below */ 
	int Error;

	int test_data[2];

	vector<int>  mask_points;
	vector<int> buffer_points;

	while(1)
	{
		Error = MPI_Recv(test_data,2,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		//Error in RECV break
		if(Error != MPI_SUCCESS) break;
		//End of points break;
		if(status.MPI_TAG == END_MSG) break;

		if (test_data[1] == 0)
		{
			test_data[1] = 2;
			Error = MPI_Send(test_data,2,MPI_INTEGER,NNodes-1,PIPE_MSG,MPI_COMM_WORLD);
		}
		else if(test_data[1] == 1)
		{
			mask_points.push_back(test_data[0]);
		}
		else if(test_data[1] == 2)
		{
			test_data[1] = 2;
			buffer_points.push_back(test_data[0]);
			Error = MPI_Send(test_data,2,MPI_INTEGER,NNodes-1,PIPE_MSG,MPI_COMM_WORLD);
		}
	}

	for(unsigned int i=0;i<mask_points.size();i++)
	{
		LASPointH test_point = las_points[mask_points[i]];
		test_data[0] = i;
		test_data[1] = 3;
		Error = MPI_Send(test_data,2,MPI_INTEGER,NNodes-1,PIPE_MSG,MPI_COMM_WORLD);
	}
}

/*******************************************/
/* NodeN keeps receiving LASPoints and	   */
/* writes them to the open LASWriterHandle */
/* no thinking needed except ensuring data */
/* size may be							   */
/*******************************************/
void NodeN()
{
	unsigned int data_size = las_points.size();
	int test_data[2];
	MPI_Status Status;
	for(unsigned int i=0;i<data_size-1;i++)
	{
		MPI_Recv(test_data,2,MPI_INTEGER,MPI_ANY_SOURCE,PIPE_MSG,MPI_COMM_WORLD,&Status);
		int test_index = test_data[0];
		int test_class = test_data[1];

		LASPointH test_point = las_points[test_index];

		LASPoint_SetClassification(test_point,test_class);
		LASWriter_WritePoint(las_writer,test_point);
	}
}

/****************************************/
/* Job locator based on Node number     */
/* basically a switch case with if-else */
/****************************************/
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


/*********************************************/
/* Very simplistic main method could do more */
/*********************************************/
int main(int argc, char ** argv) {

	Init(argc,argv);
	NodeJobs();
	Finalize();
}