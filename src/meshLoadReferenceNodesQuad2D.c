#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, "nodes/quadrilateralN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Ncheck;
  sscanf(buf, "%d", &Ncheck);
  if(Ncheck != N) printf("bu55er - wrong data file\n");
  mesh->N = N;
  mesh->Nfp = N+1;

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Npcheck;
  sscanf(buf, "%d", &Npcheck);
  mesh->Np = Npcheck;

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->r = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->r+n, mesh->s+n);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Dr = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Dr+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Ds = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Ds+n);
  }
  fgets(buf, BUFSIZ, fp); // read EOL

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->faceNodes = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));
  for(int f=0;f<mesh->Nfaces;++f){
    for(int n=0;n<mesh->Nfp;++n){
      fscanf(fp, iintFormat, mesh->faceNodes+n + f*mesh->Nfp);
    }
  }
  fgets(buf, BUFSIZ, fp); // read EOL

    
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->LIFT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Nfaces*mesh->Nfp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->LIFT+n);
  }
  fgets(buf, BUFSIZ, fp);

  /* 1D collocation differentiation matrix on GLL nodes */
  fgets(buf, BUFSIZ, fp); // read comment
  printf("comment = %s\n", buf);
  mesh->D = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    for(int m=0;m<mesh->N+1;++m){
      fscanf(fp, dfloatFormat, mesh->D+m+n*(mesh->N+1));
      printf("%f ",  mesh->D[m+n*(mesh->N+1)]);
    }
    printf("\n");
  }
  fgets(buf, BUFSIZ, fp); // read comment

  /* 1D GLL node coordinates */
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->gllz = (dfloat*) calloc(mesh->N+1, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    fscanf(fp, dfloatFormat, mesh->gllz+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fclose(fp);
}
