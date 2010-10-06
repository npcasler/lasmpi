
/* this include file is mandatory */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* MPI sample program; NOT INTENDED TO BE EFFICIENT as a prime
   finder, either in algorithm or implementation 

   MPI (Message Passing Interface) is a popular package using
   the "message passing" paradigm for communicating between
   processors in parallel applications; as the name implies,
   processors communicate by passing messages using "send" and 
   "receive" functions

   finds and reports the number of primes less than or equal to N

   uses a pipeline approach:  node 0 looks at all the odd numbers
   (i.e. has already done filtering out of multiples of 2) and
   filters out those that are multiples of 3, passing the rest
   to node 1; node 1 filters out the multiples of 5, passing
   the rest to node 2; in this simple example, we just have node
   2 filter out all the rest and then report the number of primes 

   note that we should NOT have a node run through all numbers
   before passing them on to the next node, since we would then
   have no parallelism at all; on the other hand, passing on just
   one number at a time isn't efficient either, due to the high
   overhead of sending a message if it is a network (tens of
   microseconds until the first bit reaches the wire, due to
   software delay); thus efficiency would be greatly improved if 
   each node saved up a chunk of numbers before passing them to 
   the next node */

#define MAX_N 100000
#define PIPE_MSG 0  /* type of message containing a number to 
                       be checked */
#define END_MSG 1  /* type of message indicating no more data will
                      be coming */

int NNodes,  /* number of nodes in computation*/
    N,  /* find all primes from 2 to N */
    Me,  /* my node number */
    ToCheck;  /* current number to check for passing on to next node;
                 stylistically this might be nicer as a local in
		 Node*(), but I have placed it here to dramatize
		 the fact that the globals are NOT shared among
		 the nodes, thus no problem */

double T1,T2;  /* start and finish times */

Init(Argc,Argv)
   int Argc;  char **Argv;

{  int DebugWait;  

   /* mandatory to begin any MPI program */
   MPI_Init(&Argc,&Argv);

   /* do not access command-line arguments until after MPI_Init()
      is called */
   /* get command-line arguments */
   N = atoi(Argv[1]); 
   DebugWait = atoi(Argv[2]);

   /* this loop is here to synchronize all nodes for debugging;
      if DebugWait is specified as 1 on the command line, all nodes 
      wait here until the debugging programmer starts gdb at all 
      nodes and within gdb sets DebugWait to 0 to then proceed; 
      see http://heather.cs.ucdavis.edu/~matloff/pardebug.html */
   while (DebugWait) ;

   /* puts the number of nodes in NNodes */
   MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
   /* puts the node number of this node in Me */
   MPI_Comm_rank(MPI_COMM_WORLD,&Me); 

   /* OK, get started; first record current time in T1 */
   if (Me == 2) T1 = MPI_Wtime();
}

Node0()

{  int I,Dummy,
       Error;  /* not checked in this example */
   for (I = 1; I <= N/2; I++)  {
      ToCheck = 2 * I + 1;
      if (ToCheck > N) break;
      /* MPI_Send  --  send a message
               parameters:
               pointer to place where message is to be drawn from
               number of items in message
               item type
               destination node
               message type ("tag") programmer-defined
               node group number (in this case all nodes) */
      if (ToCheck % 3 > 0)
         Error = MPI_Send(&ToCheck,1,MPI_INT,1,PIPE_MSG,MPI_COMM_WORLD);
   }
   Error = MPI_Send(&Dummy,1,MPI_INT,1,END_MSG,MPI_COMM_WORLD);
}

Node1()

{  int ToCheck,  /* current number to check from Node 0 */
       Error,  /* not checked in this example */
       Dummy;
   MPI_Status Status;  /* see below */ 

   while (1)  {
      /* MPI_Recv  --  receive a message
               parameters:
               pointer to place to store message
               number of items in message (see notes on
                  this at the end of this file)
               item type
               accept message from which node(s) 
               message type ("tag"), programmer-defined (in this
	          case any type)
               node group number (in this case all nodes)
               status (see notes on this at the end of this file) */
      Error = MPI_Recv(&ToCheck,1,MPI_INT,0,MPI_ANY_TAG,
                  MPI_COMM_WORLD,&Status);
      if (Status.MPI_TAG == END_MSG) break;
      if (ToCheck % 5 > 0) 
         Error = MPI_Send(&ToCheck,1,MPI_INT,2,PIPE_MSG,MPI_COMM_WORLD);
   }
   /* now send our end-of-data signal, which is conveyed in the
      message type, not the message (we have a dummy message just
      as a placeholder */
   Error = MPI_Send(&Dummy,1,MPI_INT,2,END_MSG,MPI_COMM_WORLD);
}

Node2()

{  int ToCheck,  /* current number to check from Node 0 */
       Error,  /* not checked in this example */
       PrimeCount,I,IsComposite;
   MPI_Status Status;  /* see below */ 

   PrimeCount = 3;  /* must account for the primes 2, 3 and 5, which
                       won't be detected below */
   while (1)  {
      Error = MPI_Recv(&ToCheck,1,MPI_INT,1,MPI_ANY_TAG,
                  MPI_COMM_WORLD,&Status);
      if (Status.MPI_TAG == END_MSG) break;
      IsComposite = 0;
      for (I = 7; I*I <= ToCheck; I += 2)
         if (ToCheck % I == 0)  {
	    IsComposite = 1;
	    break;
	 }
      if (!IsComposite) PrimeCount++;  
   }
   /* check the time again, and subtract to find run time */
   T2 = MPI_Wtime();
   printf("elapsed time = %f\n",(float)(T2-T1));
   /* print results */
   printf("number of primes = %d\n",PrimeCount);
}

main(argc,argv)
   int argc; char **argv;

{  Init(argc,argv);
   /* note:  instead of having a switch statement, we could write
      three different programs, each running on a different node */
   switch (Me)  {
      case 0:  Node0();
	  break;
      case 1:  Node1();
          break;
      case 2:  Node2();
   };
   /* mandatory for all MPI programs */
   MPI_Finalize();
}

/* explanation of "number of items" and "status" arguments at the end 
   of MPI_Recv():

   when receiving a message you must anticipate the longest possible
   message, but the actual received message may be much shorter than
   this; you can call the MPI_Get_count() function on the status
   argument to find out how many items were actually received

   the status argument will be a pointer to a struct, containing the
   node number, message type and error status of the received
   message

   say our last parameter is Status; then Status.MPI_SOURCE
   will contain the number of the sending node, and 
   Status.MPI_TAG will contain the message type; these are
   important if used MPI_ANY_SOURCE or MPI_ANY_TAG in our
   node or tag fields but still have to know who sent the
   message or what kind it is */


