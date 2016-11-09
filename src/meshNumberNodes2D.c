#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

// create a global numbering of nodes
// [  requires (EToE,EToF) and halo information ]
void meshNumberNodes2D(mesh2D *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // number of local and halo nodes
  iint NlocalNodes = mesh->Nelements*mesh->Np;
  iint NhaloNodes  = mesh->totalHaloPairs*mesh->Np;

  // number of local nodes on each process
  iint *allNlocalNodes = (iint*) calloc(size, sizeof(iint));
  MPI_Allgather(&NlocalNodes, 1, MPI_IINT,
		allNlocalNodes, 1, MPI_IINT,
		MPI_COMM_WORLD);

  iint *cumulativeNlocalNodes = (iint*) calloc(size, sizeof(iint));
  for(iint r=1;r<size;++r){
    cumulativeNlocalNodes[r] = cumulativeNlocalNodes[r-1] + allNlocalNodes[r-1];
  }
  
  // build initial index array
  iint *globalNumbering
    = (iint*) calloc( (NlocalNodes+NhaloNodes), sizeof(iint));

  // start using accumulated count of nodes on lower rank processes
  iint cnt = cumulativeNlocalNodes[rank];
  for(iint e=0;e<mesh->Nelements;++e)
    for(iint n=0;n<mesh->Np;++n)
      globalNumbering[e*mesh->Np+n] = cnt++;
  
  // halo send buffer for globalNumbering
  iint *sendBuffer = (iint*) calloc(NhaloNodes, sizeof(iint));
  iint maxCnt;
  
  do{  // keep comparing face node numbers

    // swap halos between partitions
    meshHaloExchange2D(mesh,mesh->Np*sizeof(iint),
		       globalNumbering,sendBuffer,globalNumbering+NlocalNodes);

    // loop over face nodes
    cnt = 0;
    for(int n=0;n<mesh->Nelements*mesh->Nfaces*mesh->Nfp;++n){
      iint idM = mesh->vmapM[n];
      iint idP = mesh->vmapP[n];
      iint gnM = globalNumbering[idM];
      iint gnP = globalNumbering[idP];
      // use maximum of current globalNumbering for either trace
      if(gnM!=gnP){
	globalNumbering[idM] = mymax(gnM,gnP);
	globalNumbering[idP] = mymax(gnM,gnP);
	++cnt;
      }
    }

    // find maximum number of changes made by any process
    MPI_Allreduce(&cnt, &maxCnt, 1, MPI_IINT, MPI_MAX, MPI_COMM_WORLD);
    if(rank==0) printf("number of changes = %d\n", maxCnt);

  }while(maxCnt>0); // until no face node numbers change

  // at this point we can set up a communication pattern
  // for finite element gather scatters
  // [ can use global numbers to determine process
  //   that owns each global node. May be lopsided ]
  mesh->globalNumbering = globalNumbering;
  
}
