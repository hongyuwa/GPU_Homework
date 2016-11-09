#include "mesh2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void meshAcousticsVolume2D(mesh2D *mesh){

  // for all elements
  for(iint e=0;e<mesh->Nelements;++e){

    // compute gradient at each node
    for(iint j=0;j<mesh->N+1;++j){
      for(iint i=0;i<mesh->N+1;++i){

	// local node index
	iint n = i + (mesh->N+1)*j;

	// prefetch geometric factors (quadrilateral)
	iint gid = mesh->Np*mesh->Nvgeo*e + n;
        dfloat drdx = mesh->vgeo[gid + mesh->Np*RXID];
        dfloat drdy = mesh->vgeo[gid + mesh->Np*RYID];
        dfloat dsdx = mesh->vgeo[gid + mesh->Np*SXID];
        dfloat dsdy = mesh->vgeo[gid + mesh->Np*SYID];

        // compute 'r' and 's' derivatives of (u,v,p) at node n
        dfloat dudr = 0, duds = 0;
        dfloat dvdr = 0, dvds = 0;
        dfloat dpdr = 0, dpds = 0;

        for(iint m=0;m<mesh->N+1;++m){
	  dudr +=mesh->D[i*(mesh->N+1) + m]*mesh->q[0 + m*mesh->Nfields + j*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];

	  duds +=mesh->D[j*(mesh->N+1) + m]*mesh->q[0 + i*mesh->Nfields + m*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];

	  dvdr +=mesh->D[i*(mesh->N+1) + m]*mesh->q[1 + m*mesh->Nfields + j*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];

	  dvds +=mesh->D[j*(mesh->N+1) + m]*mesh->q[1 + i*mesh->Nfields + m*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];

	  dpdr +=mesh->D[i*(mesh->N+1) + m]*mesh->q[2 + m*mesh->Nfields + j*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];

	  dpds +=mesh->D[j*(mesh->N+1) + m]*mesh->q[2 + i*mesh->Nfields + m*(mesh->N+1)*mesh->Nfields + e*mesh->Np*mesh->Nfields];
	}

	// chain rule
        dfloat dudx = drdx*dudr + dsdx*duds;
        dfloat dvdy = drdy*dvdr + dsdy*dvds;
        dfloat dpdx = drdx*dpdr + dsdx*dpds;
        dfloat dpdy = drdy*dpdr + dsdy*dpds;

        // indices for writing the RHS terms
        iint id = mesh->Nfields*(e*mesh->Np + n);

        // store acoustics rhs contributions from collocation differentiation
        mesh->rhsq[id+0] = -dpdx;
        mesh->rhsq[id+1] = -dpdy;
        mesh->rhsq[id+2] = -dudx-dvdy;
      }
    }
  }
}
