#include <stdlib.h>
#include "mesh2D.h"

// function to compute surface contributions
// for nodal DG acoustics right hand side
void meshAcousticsSurface2D(mesh2D *mesh){

  // temporary storage for flux terms
  dfloat fluxu = 0.0;
  dfloat fluxv = 0.0;
  dfloat fluxp = 0.0;

  // end weights
  dfloat w = 2.0/(mesh->N*(mesh->N+1));

  // for all elements
  for(iint e=0;e<mesh->Nelements;++e){
    // for all face nodes of all elements
    for(iint f=0;f<mesh->Nfaces;++f){ // for each face

      for(iint i=0;i<mesh->Nfp;++i){ // for each node on face
        // find index
        iint base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        dfloat nx = mesh->sgeo[base+0];
        dfloat ny = mesh->sgeo[base+1];
        dfloat sJ = mesh->sgeo[base+2];
        dfloat invJ = mesh->sgeo[base+3];

        // indices of negative and positive traces of face node
        iint id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp +i;
        iint idM = mesh->Nfields*mesh->vmapM[id];
        iint idP = mesh->Nfields*mesh->vmapP[id];

        // load negative trace node values of q
        dfloat uM = mesh->q[idM+0];
        dfloat vM = mesh->q[idM+1];
        dfloat pM = mesh->q[idM+2];

        // load positive trace node values of q
        dfloat uP = mesh->q[idP+0]; // fix BCs later
        dfloat vP = mesh->q[idP+1];
        dfloat pP = mesh->q[idP+2];

        // boundary conditions
        if(idM==idP || idP<0){
	  // assert Neumann for pressure and no penetration for velocity
	  uP = -uM;
	  vP = -vM;
	  pP =  pM;
        }

        // compute (q^* - q^-)
        dfloat duS = 0.5*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
        dfloat dvS = 0.5*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
        dfloat dpS = 0.5*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM));

        // evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
        fluxu = 1.0/w*invJ*sJ*(-nx*dpS);
        fluxv = 1.0/w*invJ*sJ*(-ny*dpS);
        fluxp = 1.0/w*invJ*sJ*(-nx*duS-ny*dvS);

        mesh->rhsq[idM]   += fluxu;
        mesh->rhsq[idM+1] += fluxv;
        mesh->rhsq[idM+2] += fluxp;
      }
    }
  }
}

