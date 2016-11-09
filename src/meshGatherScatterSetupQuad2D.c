#include <stdlib.h>
#include "mesh2D.h"

typedef struct {
  iint localId;
  iint globalId;
}gatherNode_t;

int compareGatherNodes(const void *a, const void *b){
  const iint *nodea = (iint*) a;
  const iint *nodeb = (iint*) b;

  if(nodea->globalId <nodeb->globalId) return -1;
  if(nodea->globalId >nodeb->globalId) return +1;
  return 0;
}

// serial -
// in parallel - add a isend/irecv step before doing this then finish up after wait
void meshGatherScatterSetupQuad2D(mesh2D *mesh){

  // make sure global numbering is set up
  meshNumberNodes2D(mesh);

  iint Nnodes = mesh->Np*mesh->Nelements;

  // set up local-global number pairs
  gatherNode_t *gatherNodes = (gatherNode_t*) calloc(Nnodes, sizeof(gatherNode_t));
  int lnode = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      gatherNodes[lnode].localId  = lnode;
      gatherNodes[lnode].globalId = mesh->globalNumbering[lnode];
      ++lnode;
    }
  }
  
  // sort based on gatherNode index
  qsort(gatherNodes, Nnodes, sizeof(gatherNode_t), compareGatherNodes);

  // extract: ordered list of local node indices for gather
  mesh->gatherNodeIds = (iint*) calloc(Nnodes, sizeof(iint));
  for(n=0;n<Nnodes;++n)   
    mesh->gatherNodeIds[n] = gatherNodes[n].localId;

  // count unique global nodes
  mesh->NgatherNodes = 1;
  for(n=1;n<Nnodes;++n)
    if(gatherNodes[n].globalId!=gatherNodes[n-1].globalId)
      ++(mesh->NgatherNodes);
  
  // count multiplicity of gather nodes and find offsets
  mesh->gatherCounts  = (iint*) calloc(mesh->NgatherNodes, sizeof(iint)); 
  mesh->gatherOffsets = (iint*) calloc(mesh->NgatherNodes, sizeof(iint)); 
  iint gnode = 0;
  for(n=1;n<Nnodes;++n){ 
    if(gatherNodes[n].globalId!=gatherNodes[n-1].globalId){
      ++gnode;
      mesh->gatherOffsets[gnode] = n;
      mesh->gatherCounts[gnode-1] = n-mesh->gatherOffsets[gnode-1]; 
    }
  }
  mesh->gatherCounts[gnode] = Nnodes-gatherOffsets[gnode]; // ??

  free(gatherNodes);
}
