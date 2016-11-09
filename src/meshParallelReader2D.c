
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "mpi.h"

#include "mesh2D.h"

/* 
   purpose: read gmsh triangle mesh 
*/
mesh2D* meshParallelReader2D(char *fileName){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp = fopen(fileName, "r");
  int n;

  mesh2D *mesh = (mesh2D*) calloc(1, sizeof(mesh2D));

  mesh->Nverts = 3; // number of vertices per element
  mesh->Nfaces = 3;

  if(fp==NULL){
    printf("meshReader2D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat *VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat *VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));

  /* load nodes */
  for(n=0;n<mesh->Nnodes;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat, VX+n, VY+n);
  }
  
  /* look for section with Element node data */
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* find # of triangles */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Ntriangles = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==2) ++Ntriangles;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Ntriangles/size;
  int remainder = Ntriangles - chunk*size;

  int NtrianglesLocal = chunk + (rank<remainder);

  /* where do these elements start ? */
  int start = rank*chunk + mymin(rank, remainder); 
  int end = start + NtrianglesLocal-1;
  
  /* allocate space for Element node index data */

  mesh->EToV 
    = (iint*) calloc(NtrianglesLocal*mesh->Nverts, 
		     sizeof(iint));

  /* scan through file looking for triangle elements */
  int cnt=0;
  Ntriangles = 0;
  for(n=0;n<mesh->Nelements;++n){
    iint elementType, v1, v2, v3;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType==2){  // triangle
      if(start<=Ntriangles && Ntriangles<=end){
	sscanf(buf, "%*d%*d%*d%*d%*d %d%d%d", 
	       &v1, &v2, &v3);
	/* read vertex triplet for trianngle */
	mesh->EToV[cnt*mesh->Nverts+0] = v1-1;
	mesh->EToV[cnt*mesh->Nverts+1] = v2-1;
	mesh->EToV[cnt*mesh->Nverts+2] = v3-1;
	++cnt;
      }
      ++Ntriangles;
    }
  }
  fclose(fp);

  /* record number of found triangles */
  mesh->Nelements = NtrianglesLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(n=0;n<mesh->Nverts;++n){
      mesh->EX[e*mesh->Nverts+n] = VX[mesh->EToV[e*mesh->Nverts+n]];
      mesh->EY[e*mesh->Nverts+n] = VY[mesh->EToV[e*mesh->Nverts+n]];
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);

  return mesh;

}
  
