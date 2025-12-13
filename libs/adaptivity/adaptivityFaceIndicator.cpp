/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "adaptivity.hpp"

namespace libp {

void adaptivity_t::LongestEdgeNEW(memory<dlong>& FaceFlag,
                                  memory<dlong>& RefFlag)
{
  printf("Finding Longest Edge!\n");
  dlong const MAX_REFINEMENT_LEVEL = 1;
  const dlong Nverts = mesh.Nverts; // 3
  const dlong Nfaces = mesh.Nfaces; // 3

  for (dlong e = 0; e < mesh.Nelements; ++e) {
    if (RefFlag[e] == 1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL){

    const dlong id  = e * mesh.Nverts;
    const dlong idf = e * mesh.Nfaces;

    const dfloat x0 = mesh.EX[id+0], y0 = mesh.EY[id+0];
    const dfloat x1 = mesh.EX[id+1], y1 = mesh.EY[id+1];
    const dfloat x2 = mesh.EX[id+2], y2 = mesh.EY[id+2];

    const dfloat mag0 = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)); // edge (0,1)
    const dfloat mag1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // edge (1,2)
    const dfloat mag2 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0)); // edge (2,0)

      
      dlong Face_id;
      if (mag0 >= mag1 && mag0 >= mag2) {
        Face_id = idf+0;
      } else if (mag1 >= mag2) {
        Face_id = idf+1;
      } else {
        Face_id = idf+2;
      }
      FaceFlag[Face_id] = 1;        

      
  }
  }
}


void adaptivity_t::LongestEdge(memory<dlong>& FaceFlag, memory<dlong>& RefFlag){ 
 dlong const MAX_REFINEMENT_LEVEL = 3;
  printf("Finding Longest Edge!\n");

  for (int e = 0; e < mesh.Nelements; ++e)
  {
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {      
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces; 
              
      // Find vertex locations of elements to be refined
      const dfloat x0 = mesh.EX[id+0]; const dfloat x1 = mesh.EX[id+1]; const dfloat x2 = mesh.EX[id+2];    
      const dfloat y0 = mesh.EY[id+0]; const dfloat y1 = mesh.EY[id+1]; const dfloat y2 = mesh.EY[id+2];
      
      // Find Longest Edge
      const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
      const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
      
      dlong Face_id;
      if (mag0 >= mag1 && mag0 >= mag2) {
        Face_id = idf+0;
      } else if (mag1 >= mag2) {
        Face_id = idf+1;
      } else {
        Face_id = idf+2;
      }
      FaceFlag[Face_id] = 1;        

      }          
  }           
}

} // namespace libp
