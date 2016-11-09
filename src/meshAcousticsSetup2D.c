#include "mpi.h"
#include <math.h>
#include "mesh2D.h"

void meshAcousticsSetup2D(mesh2D *mesh){

  mesh->Nfields = 4;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  mesh->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

  // fix this later (initial conditions)
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];

      acousticsCavitySolution2D(x, y, t,
				mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);

      cnt += mesh->Nfields;
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  // set time step
  dfloat hmin = 1e9;
  dfloat hmax = 1e-9;
  for(iint e=0;e<mesh->Nelements;++e){

    for(iint f=0;f<(mesh->Nfaces*mesh->Nfp);++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e*mesh->Nfp + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/4,   sJ/J = 2L/A = 2L/(h*L) = 2/h
      // h = 2/(sJ/J)

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
  }
  dfloat cfl = .1; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh->N+1.)*(mesh->N+1.)*mesh->Lambda2);

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  //
  mesh->finalTime = 0.1;
  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  mesh->dt = mesh->finalTime/mesh->NtimeSteps;

  // errorStep
  mesh->errorStep = 100;

  printf("h = %g\n", hmax);
  printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  dfloat w = 2.0/(mesh->N*(mesh->N+1));

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 1);
  mesh->device.setup(deviceConfig);

  // build D transpose
  dfloat *DT = (dfloat*) calloc((mesh->N+1)*(mesh->N+1), sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    for(int m=0;m<mesh->N+1;++m){
      DT[n+m*(mesh->N+1)] = mesh->D[n*(mesh->N+1)+m];
    }
  }

  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_resq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->resq);

  mesh->o_D =
    mesh->device.malloc((mesh->N+1)*(mesh->N+1)*sizeof(dfloat), mesh->D);

  mesh->o_DT =
    mesh->device.malloc((mesh->N+1)*(mesh->N+1)*sizeof(dfloat), DT);

  mesh->o_vgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Nvgeo*sizeof(dfloat),
			mesh->vgeo);

  mesh->o_sgeo =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nfp*mesh->Nsgeo*sizeof(dfloat),
			mesh->sgeo);

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(iint),
			mesh->vmapP);

  //// need to double check
  //dfloat *q4 = (dfloat*) calloc(4*mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*sizeof(dfloat), sizeof(dfloat));

  //dfloat *rhsq4 = (dfloat*) calloc(4*mesh->Np*mesh->Nelements*sizeof(dfloat), sizeof(dfloat));

  //dfloat *sgeo4 = (dfloat*) calloc(4*mesh->Np*mesh->Nelements*mesh->Nfaces*sizeof(dfloat), sizeof(dfloat));

  //// need to double check.
  //for(iint ii=0;ii<mesh->Np*(mesh->totalHaloPairs+mesh->Nelements);++ii){
  //  q4[ii*4+0]=mesh->q[ii];
  //  q4[ii*4+1]=mesh->q[ii+1];
  //  q4[ii*4+2]=mesh->q[ii+2];
  //  q4[ii*4+3]=mesh->q[ii+3];
  //}

  //for(iint ii=0;ii<mesh->Np*mesh->Nelements;++ii){
  //  rhsq4[ii*4+0]=mesh->rhsq[ii];
  //  rhsq4[ii*4+1]=mesh->rhsq[ii+1];
  //  rhsq4[ii*4+2]=mesh->rhsq[ii+2];
  //  rhsq4[ii*4+3]=mesh->rhsq[ii+3];
  //}

  //mesh->o_q4 =
  //  mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*4*sizeof(dfloat), q4);

  //mesh->o_rhsq4 =
  //  mesh->device.malloc(mesh->Np*mesh->Nelements*4*sizeof(dfloat), rhsq4);

  occa::kernelInfo kernelInfo;

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_N", mesh->N);
  kernelInfo.addDefine("p_w", w);
  kernelInfo.addDefine("p_Np", mesh->Np);
  kernelInfo.addDefine("p_Nfp", mesh->Nfp);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfp", mesh->Nfp*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  kernelInfo.addDefine("p_Lambda2", 0.5f);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
  }

  if(sizeof(iint)==4){
    kernelInfo.addDefine("iint","int");
  }
  if(sizeof(iint)==8){
    kernelInfo.addDefine("iint","long long int");
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }

  mesh->acousticsVolumeKernel =
    mesh->device.buildKernelFromSource("src/meshAcousticsVolume2D.okl",
				       "meshAcousticsVolume2D_o0",
				       kernelInfo);

  mesh->acousticsSurfaceKernel =
    mesh->device.buildKernelFromSource("src/meshAcousticsSurface2D.okl",
				       "meshAcousticsSurface2D_s0",
				       kernelInfo);

  mesh->acousticsUpdateKernel =
    mesh->device.buildKernelFromSource("src/meshAcousticsUpdate2D.okl",
				       "meshAcousticsUpdate2D",
				       kernelInfo);

  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource("src/meshHaloExtract2D.okl",
				       "meshHaloExtract2D",
				       kernelInfo);

  int Ntests = 10;


#define maxNkernels 100

  int NvolumeKernels = 6;
  occa::kernel *meshAcousticsVolumeKernels = new occa::kernel[maxNkernels];
  char kernelNames[maxNkernels][BUFSIZ];
  double bestElapsed = 1e9;

  for(int ker=0;ker<NvolumeKernels;++ker){
    sprintf(kernelNames[ker], "meshAcousticsVolume2D_o%d", ker);

    meshAcousticsVolumeKernels[ker] =
      mesh->device.buildKernelFromSource("src/meshAcousticsVolume2D.okl", kernelNames[ker], kernelInfo);

    mesh->device.finish();
    occa::tic(kernelNames[ker]);
    for(int test=0;test<Ntests;++test)
      meshAcousticsVolumeKernels[ker](mesh->Nelements,
				      mesh->o_vgeo,
				      mesh->o_DT,
				      mesh->o_q,
				      mesh->o_rhsq);
    mesh->device.finish();
    double elapsed = occa::toc(kernelNames[ker]);
    if(elapsed<bestElapsed){
      mesh->acousticsVolumeKernel = meshAcousticsVolumeKernels[0];
      printf("promoting kernel: %d (time %g)\n", ker, elapsed);
      bestElapsed = elapsed;
    }
    else{
      printf("not promoting kernel: %d (time %g)\n", ker, elapsed);
    }
  }

  int NsurfaceKernels = 7;
  char surfaceKernelNames[maxNkernels][BUFSIZ];
  occa::kernel *meshAcousticsSurfaceKernels = new occa::kernel[maxNkernels];
  bestElapsed = 1e9;

  for(int ker=0;ker<NsurfaceKernels;++ker){
    sprintf(surfaceKernelNames[ker], "meshAcousticsSurface2D_s%d", ker);

    meshAcousticsSurfaceKernels[ker] =
      mesh->device.buildKernelFromSource("src/meshAcousticsSurface2D.okl",
					 surfaceKernelNames[ker], kernelInfo);

    mesh->device.finish();
    occa::tic(surfaceKernelNames[ker]);
    for(int test=0;test<Ntests;++test)
      meshAcousticsSurfaceKernels[ker](mesh->Nelements,
				       mesh->o_sgeo,
				       mesh->o_vmapM,
				       mesh->o_vmapP,
				       mesh->o_q,
      				       mesh->o_rhsq);
    mesh->device.finish();
    double elapsed = occa::toc(surfaceKernelNames[ker]);
    if(elapsed<bestElapsed){
      mesh->acousticsSurfaceKernel = meshAcousticsSurfaceKernels[0];
      printf("promoting kernel: %d (time %g)\n", ker, elapsed);
      bestElapsed = elapsed;
    }
    else{
      printf("not promoting kernel: %d (time %g)\n", ker, elapsed);
    }
  }
}
