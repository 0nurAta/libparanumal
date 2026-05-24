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

void adaptivity_t::LongestEdge(memory<dlong>& FaceFlag, 
                               memory<dlong>& RefFlag, 
                               dlong level){

  for (int e = 0; e < mesh.Nelements; ++e)
  {
      if (RefFlag[e]==1 && EToRefLevel[e]<level)
      {
      const dlong id = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces; 
              
      // Find vertex locations of elements to be refined
      const dfloat x0 = mesh.EX[id+0]; 
      const dfloat x1 = mesh.EX[id+1]; 
      const dfloat x2 = mesh.EX[id+2];

      const dfloat y0 = mesh.EY[id+0]; 
      const dfloat y1 = mesh.EY[id+1]; 
      const dfloat y2 = mesh.EY[id+2];
      
      // Find Longest Edge
      const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
      const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
      
      dlong Face_id;
      if (mag0 >= mag1 && mag0 >= mag2) {Face_id = idf+0;} 
      else if (mag1 >= mag2) {Face_id = idf+1;} 
      else {Face_id = idf+2;}

      FaceFlag[Face_id] = 1;        
      }          
  }           
}

void adaptivity_t::LongestEdge2(memory<dlong>& FaceFlag, 
                               memory<dlong>& RefFlag, 
                               dlong level, dlong Nrefine,
                               memory<hlong>& new_v_id,
                               hlong* new_vertex){

  dlong counter = 0;
  memory<hlong> elemList(Nrefine*3,-1);
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      if (RefFlag[e]==1 && EToRefLevel[e]<level)
      {
      const dlong id = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces; 
              
      // Find vertex locations of elements to be refined
      const dfloat x0 = mesh.EX[id+0]; 
      const dfloat x1 = mesh.EX[id+1]; 
      const dfloat x2 = mesh.EX[id+2];

      const dfloat y0 = mesh.EY[id+0]; 
      const dfloat y1 = mesh.EY[id+1]; 
      const dfloat y2 = mesh.EY[id+2];
      
      // Find Longest Edge
      const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
      const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
      
      dlong Face_id;
      hlong vid0; hlong vid1;
      if (mag0 >= mag1 && mag0 >= mag2) {Face_id = idf+0; 
      vid0= mesh.EToV[idf+0]; vid1= mesh.EToV[idf+1];} 
      else if (mag1 >= mag2) {Face_id = idf+1; vid0= mesh.EToV[idf+1]; vid1= mesh.EToV[idf+2];} 
      else {Face_id = idf+2; vid0= mesh.EToV[idf+2]; vid1= mesh.EToV[idf+0];}

      if (vid0>vid1)
      {
       hlong c;
       c = vid0;
       vid0=vid1;
       vid1=c;

      }
         printf("e=%d,vid0=%d,vid1=%d\n", e,vid0,vid1);
      FaceFlag[Face_id] = 1;   
          
      elemList[counter*3+0] = e;
      elemList[counter*3+1] = vid0;
      elemList[counter*3+2] = vid1;
      counter ++; 
      }          
  } 

  // sorting algorithm 
    bool swp;
    for (int i = 0; i < counter-1; ++i)
  {
    swp = false;
    for (int j = 0; j < counter-i-1; ++j)
    {
      if (elemList[j*3+1] > elemList[(j+1)*3+1] ||
   (elemList[j*3+1] == elemList[(j+1)*3+1] &&
    elemList[j*3+2] > elemList[(j+1)*3+2])) {
        std::swap(elemList[j*3+0],elemList[(j+1)*3+0]);
        std::swap(elemList[j*3+1],elemList[(j+1)*3+1]);
        std::swap(elemList[j*3+2],elemList[(j+1)*3+2]);
        swp = true;
      }
    }
       if (!swp)
            break;
  }           

  //memory<dlong> new_v_id(counter*2,-1);
  dlong k = 0;
for (int i = 0; i < counter; ++i)
{
  if (i > 0) {
    bool sameEdge =
      (elemList[i*3+1] == elemList[(i-1)*3+1]) &&
      (elemList[i*3+2] == elemList[(i-1)*3+2]);

    if (!sameEdge) {k++;}
  }

  new_v_id[i*2+0] = elemList[i*3+0];
  new_v_id[i*2+1] = mesh.Nnodes + k;
}
*new_vertex = (counter > 0) ? (k + 1) : 0;
  printf("new_vertex=%d\n",*new_vertex );

 /*for (int i = 0; i < Nrefine; ++i)
  {

   printf("e=%d,vid0=%d,vid1=%d\n", elemList[i*3+0],elemList[i*3+1],elemList[i*3+2]);
   
  }*/
}


void adaptivity_t::NewestVertex(memory<dlong>& FaceFlag, memory<dlong>& RefFlag){ 
 dlong const MAX_REFINEMENT_LEVEL = 2;

  for (int e = 0; e < mesh.Nelements; ++e)
  {   

            if (RefFlag[e]==1 && EToRefLevel[e]<1)
      {
      const dlong id = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces; 
              
      // Find vertex locations of elements to be refined
      const dfloat x0 = mesh.EX[id+0]; 
      const dfloat x1 = mesh.EX[id+1]; 
      const dfloat x2 = mesh.EX[id+2];

      const dfloat y0 = mesh.EY[id+0]; 
      const dfloat y1 = mesh.EY[id+1]; 
      const dfloat y2 = mesh.EY[id+2];

      // Find Longest Edge
      const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
      const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
      
      dlong Face_id;
      if (mag0 >= mag1 && mag0 >= mag2) {Face_id = idf+0;} 
      else if (mag1 >= mag2) {Face_id = idf+1;} 
      else {Face_id = idf+2;}
      
      FaceFlag[Face_id] = 1;          
      } 
      else if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {      
      const dlong idf = e*mesh.Nfaces;  
      FaceFlag[idf+0] = 1;        
      }          
  }           
}

} // namespace libp
