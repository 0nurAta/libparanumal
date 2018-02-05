#include "ellipticTet3D.h"

void ellipticParallelGatherScatterTet3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  // use gather map for gather and scatter
  occaTimerTic(mesh->device,"meshParallelGatherScatter3D");
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq);
  occaTimerToc(mesh->device,"meshParallelGatherScatter3D");
  
}
