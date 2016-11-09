#include <stdlib.h>
#include "mesh2D.h"

void meshGatherScatterQuad2D(mesh2D *mesh){

  for(iint n=0;n<mesh->NgatherNodes;++n){
    const iint Ngather = mesh->gatherCounts[n];
    cosnt iint offset  = mesh->gatherOffsets[n];

    dfloat sumqn = 0;
    for(iint m=0;m<Ngather;++m){
      const iint id = mesh->gatherIds[offset+m];
      sumqn += mesh->AqL[id];
    }
    
    for(iint m=0;m<Ngather;++m){
      const iint id = mesh->gatherIds[offset+m];
      mesh->ZtZAqL[id] = sumqn;
    }
  }
}
