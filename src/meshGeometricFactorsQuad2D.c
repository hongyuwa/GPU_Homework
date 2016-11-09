
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshGeometricFactorsQuad2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 5;
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, 
				sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    int id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0];
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0];
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    
    for(iint n=0;n<mesh->Np;++n){

      /* local node coordinates */
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* Jacobian matrix */
      dfloat xr = 0.25*( (1-sn)*(xe2-xe1) + (1+sn)*(xe3-xe4) );
      dfloat xs = 0.25*( (1-rn)*(xe4-xe1) + (1+rn)*(xe3-xe2) );
      dfloat yr = 0.25*( (1-sn)*(ye2-ye1) + (1+sn)*(ye3-ye4) );
      dfloat ys = 0.25*( (1-rn)*(ye4-ye1) + (1+rn)*(ye3-ye2) );

      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr*ys - xs*yr;
      dfloat rx =  ys/J;
      dfloat ry = -xs/J;
      dfloat sx = -yr/J;
      dfloat sy =  xr/J;
      
      /* store geometric factors */
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] = rx;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] = ry;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] = sx;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] = sy;
      mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]  = J;
    }
  }
}
