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

}

void adaptivity_t::RGB_flag(memory<dlong>& FaceFlag, 
                            memory<dlong>& RefFlag, 
                            dlong level, dlong Nrefine,
                            memory<hlong>& new_v_id,
                            memory<dlong>& ConfRed,
                            memory<dlong>& ConfGreen,
                            memory<hlong>& EtoNewV,
                            hlong* new_vertex)
{
  // --------------------------------------------------------------------------
  // Pass 0: clear temporary counters
  // ConfRed[e] counts how many sides of element e are marked by red neighbors
  // --------------------------------------------------------------------------
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    ConfRed[e] = 0;
  }

  // --------------------------------------------------------------------------
  // Pass 1: from currently refined red elements, count neighbor-side requests
  // If an element is touched from all 3 sides, promote it to red refinement
  // --------------------------------------------------------------------------
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e] == 1 && EToRefLevel[e] < level)
    {
      const dlong idf = e*mesh.Nfaces;

      dlong n0 = mesh.EToE[idf+0];
      dlong n1 = mesh.EToE[idf+1];
      dlong n2 = mesh.EToE[idf+2];

      if (n0 > -1) ConfRed[n0]++;
      if (n1 > -1) ConfRed[n1]++;
      if (n2 > -1) ConfRed[n2]++;
      //if (EToRefLevel[e]>=1 && (RedFlag[n0]!=1 || RedFlag[n1]!=1 || RedFlag[n2]!=1))
      //{RefFlag[e]=0;}
    }
  }

  // promote fully surrounded elements to red refinement
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (ConfRed[e] >= 3 && RefFlag[e] != 1 && EToRefLevel[e] < level)
    {
      RefFlag[e] = 1;
      printf("RGB promoted e=%d by 3-side conform marking\n", e);
    }
  }

  

  // --------------------------------------------------------------------------
  // Count final number of refined elements after promotion
  // --------------------------------------------------------------------------
  dlong finalRefineCount = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e] == 1 && EToRefLevel[e] < level)
    {
      finalRefineCount++;
    }
  }

  // --------------------------------------------------------------------------
  // Build final edge request list:
  // one row per edge request = [elem, local_edge, vmin, vmax]
  // --------------------------------------------------------------------------
  dlong counter = 0;
  memory<hlong> edgeList(3*finalRefineCount*4, -1);

  for (int i = 0; i < mesh.Nelements*mesh.Nfaces; ++i)
  {
    FaceFlag[i] = 0;
  }

  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e] == 1 && EToRefLevel[e] < level)
    {
      const dlong id  = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces;

      hlong v0 = mesh.EToV[id+0];
      hlong v1 = mesh.EToV[id+1];
      hlong v2 = mesh.EToV[id+2];

      hlong a, b;

      // local edge 0 = (v0,v1)
      a = (v0 < v1) ? v0 : v1;
      b = (v0 > v1) ? v0 : v1;
      edgeList[counter*4+0] = e;
      edgeList[counter*4+1] = 0;
      edgeList[counter*4+2] = a;
      edgeList[counter*4+3] = b;
      FaceFlag[idf+0] = 1;
      counter++;

      // local edge 1 = (v1,v2)
      a = (v1 < v2) ? v1 : v2;
      b = (v1 > v2) ? v1 : v2;
      edgeList[counter*4+0] = e;
      edgeList[counter*4+1] = 1;
      edgeList[counter*4+2] = a;
      edgeList[counter*4+3] = b;
      FaceFlag[idf+1] = 1;
      counter++;

      // local edge 2 = (v2,v0)
      a = (v2 < v0) ? v2 : v0;
      b = (v2 > v0) ? v2 : v0;
      edgeList[counter*4+0] = e;
      edgeList[counter*4+1] = 2;
      edgeList[counter*4+2] = a;
      edgeList[counter*4+3] = b;
      FaceFlag[idf+2] = 1;
      counter++;
    }
  }

  // --------------------------------------------------------------------------
  // Initialize per-element midpoint-id output:
  // [elem, mid01, mid12, mid20]
  // --------------------------------------------------------------------------
  dlong elemCounter = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e] == 1 && EToRefLevel[e] < level)
    {
      new_v_id[elemCounter*4+0] = e;
      new_v_id[elemCounter*4+1] = -1;
      new_v_id[elemCounter*4+2] = -1;
      new_v_id[elemCounter*4+3] = -1;
      elemCounter++;
    }
  }

  // --------------------------------------------------------------------------
  // Sort by canonical edge (vmin, vmax)
  // --------------------------------------------------------------------------
  bool swp;
  for (int i = 0; i < counter-1; ++i)
  {
    swp = false;
    for (int j = 0; j < counter-i-1; ++j)
    {
      if (edgeList[j*4+2] > edgeList[(j+1)*4+2] ||
         (edgeList[j*4+2] == edgeList[(j+1)*4+2] &&
          edgeList[j*4+3] > edgeList[(j+1)*4+3]))
      {
        std::swap(edgeList[j*4+0], edgeList[(j+1)*4+0]);
        std::swap(edgeList[j*4+1], edgeList[(j+1)*4+1]);
        std::swap(edgeList[j*4+2], edgeList[(j+1)*4+2]);
        std::swap(edgeList[j*4+3], edgeList[(j+1)*4+3]);
        swp = true;
      }
    }
    if (!swp) break;
  }

  // --------------------------------------------------------------------------
  // Assign one midpoint id per unique edge group
  // If any row in the group already has a stored midpoint in EtoNewV,
  // use that id for the whole group.
  // EtoNewV layout:
  // [elem, mid01, mid12, mid20]
  // --------------------------------------------------------------------------
  dlong newCount = 0;

  for (int i = 0; i < counter; )
  {
    dlong j = i + 1;

    for (; j < counter; ++j)
    {
      bool sameEdge =
        (edgeList[j*4+2] == edgeList[i*4+2]) &&
        (edgeList[j*4+3] == edgeList[i*4+3]);

      if (!sameEdge) break;
    }

    hlong group_mid = -1;

    // First check whether any element in this edge group already has a midpoint id
    for (int p = i; p < j; ++p)
    {
      hlong elem       = edgeList[p*4+0];
      hlong local_edge = edgeList[p*4+1];

      if (EtoNewV[4*elem + 1 + local_edge] != -1)
      {
        group_mid = EtoNewV[4*elem + 1 + local_edge];
        printf("mid id is = %d \n",group_mid);
        break;
      }
    }

    // If not found, create a fresh midpoint id for the whole group
    if (group_mid == -1)
    {
      group_mid = mesh.Nnodes + newCount;
      newCount++;
    }

    // Assign the same midpoint id to every row in this edge group
    for (int p = i; p < j; ++p)
    {
      hlong elem       = edgeList[p*4+0];
      hlong local_edge = edgeList[p*4+1];

      for (int r = 0; r < elemCounter; ++r)
      {
        if (new_v_id[r*4+0] == elem)
        {
          new_v_id[r*4 + 1 + local_edge] = group_mid;
          break;
        }
      }
    }

    i = j;
  }

  *new_vertex = newCount;

  printf("RGB finalRefineCount=%d new_vertex=%lld\n",
         finalRefineCount, (long long)*new_vertex);

  for (int i = 0; i < elemCounter; ++i)
  {
    printf("e=%lld, v01=%lld, v12=%lld, v20=%lld\n",
           (long long)new_v_id[i*4+0],
           (long long)new_v_id[i*4+1],
           (long long)new_v_id[i*4+2],
           (long long)new_v_id[i*4+3]);
  }
}

void adaptivity_t::LongestEdgeConform(memory<dlong>& FaceFlag, 
                                      memory<dlong>& RefFlag,
                                      memory<dlong>& ConfFlag, 
                                      dlong level, dlong Nrefine,
                                      memory<hlong>& new_v_id,
                                      hlong* new_vertex){

  
  dlong counter = 0;
  memory<hlong> elemList(Nrefine*3,-1);
  memory<dlong> conf_face(Nrefine,-1);
  memory<hlong> hang_v_id(Nrefine,-1);
  memory<dlong> Face_id(Nrefine,-1);
 
  printf("ConfFlag List\n");
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      if (ConfFlag[e]>-1 && EToRefLevel[e]<(level+3))
      {
      const dlong id = e*mesh.Nverts;
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2];
      const dlong sib_e2 = PToC[ef*(level+3)+EToRefLevel[ef]-1]; 
      printf("e=%d,ef=%d\n",e,ef );
      const dlong idf = ConfFlag[e]*mesh.Nverts;
      const dlong sib_id = sib_e*mesh.Nverts; 
      const dlong sib_id2 = sib_e2*mesh.Nverts;
      const float tol = 1e-10;
      dlong e_c =-1;        
     
      dfloat new_x;
      dfloat new_y;
// Face 0
if (e_c == -1 && mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
{
  new_x = 0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]);
  new_y = 0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]);

  if(abs(mesh.EX[id+2]-new_x)<tol&&abs(mesh.EY[id+2]-new_y)<tol) {   
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  } 
  else if (abs(mesh.EX[id+0]-new_x)<tol&&abs(mesh.EY[id+0]-new_y)<tol) {
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-new_x)<tol&&abs(mesh.EY[id+1]-new_y)<tol) {
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+0]==-1 && mesh.EToB[sib_id2+0]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }       
}

// Face 1        
if (e_c == -1 && mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
{    
  if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+0]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+1]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+1]==-1 && mesh.EToB[sib_id+1]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }     
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+1]==-1 && mesh.EToB[sib_id2+1]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

// Face 2  
if (e_c == -1 && mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
{  
  if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+0]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+1]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+2]==-1 && mesh.EToB[sib_id+2]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+2]==-1 && mesh.EToB[sib_id2+2]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id +0]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id +1]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}
if (ConfFlag[e] > -1 && e_c == -1) {
  printf("FAILED CONF REQUEST: e=%d ef=%d sib_e=%d sib_e2=%d level_e=%d level_ef=%d\n",
         e, ef, sib_e, sib_e2, EToRefLevel[e], EToRefLevel[ef]);

  for (int k = 0; k < 3; ++k) {
    printf("Bisected Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[id+k], mesh.EX[id+k], mesh.EY[id+k]);
  }

  for (int k = 0; k < 3; ++k) {
    printf("Bisected Element Sibling local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[PCS[e*3+2]*3+k], mesh.EX[PCS[e*3+2]*3+k], mesh.EY[PCS[e*3+2]*3+k]);
  }

    for (int k = 0; k < 3; ++k) {
    printf("Targeted Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[idf+k], mesh.EX[idf+k], mesh.EY[idf+k]);
  }

      for (int k = 0; k < 3; ++k) {
    printf("Targeted Sib Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[sib_id+k], mesh.EX[sib_id+k], mesh.EY[sib_id+k]);
  }

      if (sib_id2>-1)
       {
         for (int k = 0; k < 3; ++k) {
    printf("Targeted Child Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[sib_id2+k], mesh.EX[sib_id2+k], mesh.EY[sib_id2+k]);
  }
       } 

    for (int k = 0; k < 3; ++k) {
    printf("  Neighbor[%d]: %d\n",
           k, (long long)mesh.EToE[id+k]);
  }

  // print face midpoints of ef / sib_e / sib_e2
}
if (e_c > -1)
      {
        const dlong idc = e_c*mesh.Nverts;

        const dfloat x0 = mesh.EX[idc+0];
        const dfloat x1 = mesh.EX[idc+1];
        const dfloat x2 = mesh.EX[idc+2];

        const dfloat y0 = mesh.EY[idc+0];
        const dfloat y1 = mesh.EY[idc+1];
        const dfloat y2 = mesh.EY[idc+2];

        const dfloat mag0 = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        const dfloat mag1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        const dfloat mag2 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));

        hlong vid0, vid1;

        if (mag0 >= mag1 && mag0 >= mag2) {
          Face_id[counter] = e_c*mesh.Nfaces + 0;
          vid0 = mesh.EToV[idc+0];
          vid1 = mesh.EToV[idc+1];
        }
        else if (mag1 >= mag2) {
          Face_id[counter] = e_c*mesh.Nfaces + 1;
          vid0 = mesh.EToV[idc+1];
          vid1 = mesh.EToV[idc+2];
        }
        else {
          Face_id[counter] = e_c*mesh.Nfaces + 2;
          vid0 = mesh.EToV[idc+2];
          vid1 = mesh.EToV[idc+0];
        }

        if (vid0 > vid1)
        {
          hlong c = vid0;
          vid0 = vid1;
          vid1 = c;
        }

        elemList[counter*3+0] = e_c;
        elemList[counter*3+1] = vid0;
        elemList[counter*3+2] = vid1;
        FaceFlag[Face_id[counter]] = 1; 
        if (Face_id[counter] == conf_face[counter])
        {
          ConfFlag[e]=-1;
          printf("elements %lld and %lld are conformed! \n",e,e_c );
        }
        /*printf("raw e=%d, vid0=%lld, vid1=%lld, conf_face=%d, Face_id=%d, hang=%lld\n",
               e_c,
               (long long)vid0,
               (long long)vid1,
               conf_face[counter],
               Face_id[counter],
               (long long)hang_v_id[counter]);*/

        counter++;
      } else {ConfFlag[e]=-1;}
    }
  }
// Eliminate duplicates
for (int i = 0; i < counter; ++i)
{
  if (elemList[i*3+0] < 0) continue;

  for (int j = i+1; j < counter; ++j)
  {
    if (elemList[j*3+0] < 0) continue;

    if (elemList[j*3+0] == elemList[i*3+0])
    {
      bool iMatchesLE = (conf_face[i] == Face_id[i]);
      bool jMatchesLE = (conf_face[j] == Face_id[j]);

      if (!iMatchesLE && jMatchesLE)
      {
        // j is better, so discard i and replace it with j
        elemList[i*3+0] = elemList[j*3+0];
        elemList[i*3+1] = elemList[j*3+1];
        elemList[i*3+2] = elemList[j*3+2];
        conf_face[i]    = conf_face[j];
        hang_v_id[i]    = hang_v_id[j];
        Face_id[i]      = Face_id[j];

        // discard j
        elemList[j*3+0] = -1;
        elemList[j*3+1] = -1;
        elemList[j*3+2] = -1;
        conf_face[j]    = -1;
        hang_v_id[j]    = -1;
        Face_id[j]      = -1;
      }
      else
      {
        // keep i (either i already matches LE, or neither matches)
        elemList[j*3+0] = -1;
        elemList[j*3+1] = -1;
        elemList[j*3+2] = -1;
        conf_face[j]    = -1;
        hang_v_id[j]    = -1;
        Face_id[j]      = -1;
      }
    }
  }
}
// compact
dlong compactCounter = 0;

  memory<hlong> elemList2(counter*3,-1);
  memory<dlong> conf_face2(counter,-1);
  memory<hlong> hang_v_id2(counter,-1);
  memory<dlong> Face_id2(counter,-1);

  for (int i = 0; i < counter; ++i)
  {
    if (elemList[i*3+0] > -1)
    {
      elemList2[compactCounter*3+0] = elemList[i*3+0];
      elemList2[compactCounter*3+1] = elemList[i*3+1];
      elemList2[compactCounter*3+2] = elemList[i*3+2];
      conf_face2[compactCounter]    = conf_face[i];
      hang_v_id2[compactCounter]    = hang_v_id[i];
      Face_id2[compactCounter]   = Face_id[i];
      compactCounter++;
    }
  }

  counter = compactCounter;

  // sorting algorithm 
    bool swp;
  for (int i = 0; i < counter-1; ++i)
  {
    swp = false;
    for (int j = 0; j < counter-i-1; ++j)
    {
      if (elemList2[j*3+1] > elemList2[(j+1)*3+1] ||
         (elemList2[j*3+1] == elemList2[(j+1)*3+1] &&
          elemList2[j*3+2] > elemList2[(j+1)*3+2]))
      {
        std::swap(elemList2[j*3+0], elemList2[(j+1)*3+0]);
        std::swap(elemList2[j*3+1], elemList2[(j+1)*3+1]);
        std::swap(elemList2[j*3+2], elemList2[(j+1)*3+2]);

        std::swap(conf_face2[j],  conf_face2[j+1]);
        std::swap(hang_v_id2[j],  hang_v_id2[j+1]);
        std::swap(Face_id2[j], Face_id2[j+1]);

        swp = true;
      }
    }
    if (!swp) break;
  }

 dlong newCount = 0;

  for (int i = 0; i < counter; ++i)
  {
    if (i > 0)
    {
      bool sameEdge =
        (elemList2[i*3+1] == elemList2[(i-1)*3+1]) &&
        (elemList2[i*3+2] == elemList2[(i-1)*3+2]);

      if (sameEdge)
      {
        // inherit same node id as previous request
        new_v_id[i*2+0] = elemList2[i*3+0];
        new_v_id[i*2+1] = new_v_id[(i-1)*2+1];
        continue;
      }
    }

    // first row of a new edge group
    new_v_id[i*2+0] = elemList2[i*3+0];

    if (conf_face2[i] == Face_id2[i])
    {
      // hanging node already lies on split edge
      new_v_id[i*2+1] = hang_v_id2[i];
    }
    else
    {
      // create a new midpoint node for this edge group
      new_v_id[i*2+1] = mesh.Nnodes + newCount;
      newCount++;
    }
  }

  *new_vertex = newCount;

  printf("new_vertex=%lld,counter=%d\n", (long long)*new_vertex,counter);

  for (int i = 0; i < counter; ++i)
  {
    printf("e=%lld, vid0=%lld, vid1=%lld, conf_face=%d, split_face=%d, newid=%lld\n",
           (long long)elemList2[i*3+0],
           (long long)elemList2[i*3+1],
           (long long)elemList2[i*3+2],
           conf_face2[i],
           Face_id2[i],
           (long long)new_v_id[i*2+1]);
  }

}

void adaptivity_t::NewestVertexConform(memory<dlong>& FaceFlag, 
                                      memory<dlong>& RefFlag,
                                      memory<dlong>& ConfFlag, 
                                      dlong level, dlong Nrefine,
                                      memory<hlong>& new_v_id,
                                      hlong* new_vertex){

  
  dlong counter = 0;
  memory<hlong> elemList(Nrefine*3,-1);
  memory<dlong> conf_face(Nrefine,-1);
  memory<hlong> hang_v_id(Nrefine,-1);
  memory<dlong> Face_id(Nrefine,-1);
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      if (ConfFlag[e]>-1 && EToRefLevel[e]<(level+3))
      {
      const dlong id = e*mesh.Nverts;
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2];
      const dlong sib_e2 = PToC[ef*(level+3)+EToRefLevel[ef]-1]; 
      printf("e=%d,ef=%d\n",e,ef );
      const dlong idf = ConfFlag[e]*mesh.Nverts;
      const dlong sib_id = sib_e*mesh.Nverts; 
      const dlong sib_id2 = sib_e2*mesh.Nverts;
      const float tol = 1e-10;
      dlong e_c =-1;        
     
      dfloat new_x;
      dfloat new_y;
// Face 0
if (e_c == -1 && mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
{
  new_x = 0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]);
  new_y = 0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]);

  if(abs(mesh.EX[id+2]-new_x)<tol&&abs(mesh.EY[id+2]-new_y)<tol) {   
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  } 
  else if (abs(mesh.EX[id+0]-new_x)<tol&&abs(mesh.EY[id+0]-new_y)<tol) {
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-new_x)<tol&&abs(mesh.EY[id+1]-new_y)<tol) {
    e_c = ef;
    conf_face[counter] = idf+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+0]==-1 && mesh.EToB[sib_id2+0]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id2+0]+mesh.EX[sib_id2+1]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+0]+mesh.EY[sib_id2+1]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+0;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }       
}

// Face 1        
if (e_c == -1 && mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
{    
  if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+0]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+1]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<tol) { 
    e_c = ef;
    conf_face[counter] = idf+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+1]==-1 && mesh.EToB[sib_id+1]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<tol) { 
    e_c = sib_e;
    conf_face[counter] = sib_id+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }     
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+1]==-1 && mesh.EToB[sib_id2+1]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id2+1]+mesh.EX[sib_id2+2]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+1]+mesh.EY[sib_id2+2]))<tol) { 
    e_c = sib_e2;
    conf_face[counter] = sib_id2+1;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

// Face 2  
if (e_c == -1 && mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
{  
  if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+0]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if(abs(mesh.EX[id+1]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<tol) {   
    e_c = ef;
    conf_face[counter] = idf+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && mesh.EToF[sib_id+2]==-1 && mesh.EToB[sib_id+2]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+0]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id+1]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<tol) { 
    e_c = sib_e;  
    conf_face[counter] = sib_id+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}

if (e_c == -1 && sib_e2 > -1 && mesh.EToF[sib_id2+2]==-1 && mesh.EToB[sib_id2+2]==-1){
  if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+2];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id +0]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+0]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+0];
    RefFlag[e_c] = 1;
  }
  else if (abs(mesh.EX[id +1]-0.5*(mesh.EX[sib_id2+2]+mesh.EX[sib_id2+0]))<tol&&abs(mesh.EY[id+1]-0.5*(mesh.EY[sib_id2+2]+mesh.EY[sib_id2+0]))<tol) { 
    e_c = sib_e2;  
    conf_face[counter] = sib_id2+2;
    hang_v_id[counter] = mesh.EToV[id+1];
    RefFlag[e_c] = 1;
  }
}
if (ConfFlag[e] > -1 && e_c == -1) {
  printf("FAILED CONF REQUEST: e=%d ef=%d sib_e=%d sib_e2=%d level_e=%d level_ef=%d\n",
         e, ef, sib_e, sib_e2, EToRefLevel[e], EToRefLevel[ef]);

  for (int k = 0; k < 3; ++k) {
    printf("Bisected Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[id+k], mesh.EX[id+k], mesh.EY[id+k]);
  }

  for (int k = 0; k < 3; ++k) {
    printf("Bisected Element Sibling local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[PCS[e*3+2]*3+k], mesh.EX[PCS[e*3+2]*3+k], mesh.EY[PCS[e*3+2]*3+k]);
  }

    for (int k = 0; k < 3; ++k) {
    printf("Targeted Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[idf+k], mesh.EX[idf+k], mesh.EY[idf+k]);
  }

      for (int k = 0; k < 3; ++k) {
    printf("Targeted Sib Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[sib_id+k], mesh.EX[sib_id+k], mesh.EY[sib_id+k]);
  }

      if (sib_id2>-1)
       {
         for (int k = 0; k < 3; ++k) {
    printf("Targeted Child Element  local v%d: id=%lld x=% .16e y=% .16e\n",
           k, (long long)mesh.EToV[sib_id2+k], mesh.EX[sib_id2+k], mesh.EY[sib_id2+k]);
  }
       } 

    for (int k = 0; k < 3; ++k) {
    printf("  Neighbor[%d]: %d\n",
           k, (long long)mesh.EToE[id+k]);
  }
}

if (e_c > -1)
      {


        const dlong idc = e_c*mesh.Nverts;
        hlong vid0, vid1;

       if (EToRefLevel[e_c]==0)
      {

        const dfloat x0 = mesh.EX[idc+0];
        const dfloat x1 = mesh.EX[idc+1];
        const dfloat x2 = mesh.EX[idc+2];

        const dfloat y0 = mesh.EY[idc+0];
        const dfloat y1 = mesh.EY[idc+1];
        const dfloat y2 = mesh.EY[idc+2];

        const dfloat mag0 = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        const dfloat mag1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        const dfloat mag2 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));

        

        if (mag0 >= mag1 && mag0 >= mag2) {
          Face_id[counter] = e_c*mesh.Nfaces + 0;
          vid0 = mesh.EToV[idc+0];
          vid1 = mesh.EToV[idc+1];
        }
        else if (mag1 >= mag2) {
          Face_id[counter] = e_c*mesh.Nfaces + 1;
          vid0 = mesh.EToV[idc+1];
          vid1 = mesh.EToV[idc+2];
        }
        else {
          Face_id[counter] = e_c*mesh.Nfaces + 2;
          vid0 = mesh.EToV[idc+2];
          vid1 = mesh.EToV[idc+0];
        }
      } else 
      {      
       
       Face_id[counter] = e_c*mesh.Nfaces + 0; 
       vid0 = mesh.EToV[idc+0];
       vid1 = mesh.EToV[idc+1];     
      }       

        if (vid0 > vid1)
        {
          hlong c = vid0;
          vid0 = vid1;
          vid1 = c;
        }

        elemList[counter*3+0] = e_c;
        elemList[counter*3+1] = vid0;
        elemList[counter*3+2] = vid1;
        FaceFlag[Face_id[counter]] = 1; 
        printf("raw e=%d, vid0=%lld, vid1=%lld, conf_face=%d, Face_id=%d, hang=%lld\n",
               e_c,
               (long long)vid0,
               (long long)vid1,
               conf_face[counter],
               Face_id[counter],
               (long long)hang_v_id[counter]);

                if (Face_id[counter] == conf_face[counter])
        {
          ConfFlag[e]=-1;
          printf("elements %lld and %lld are conformed! \n",e,e_c );
        }

               counter++;

        
      }else {ConfFlag[e]=-1;}
    }
  }
// Eliminate duplicates
for (int i = 0; i < counter; ++i)
{
  if (elemList[i*3+0] < 0) continue;

  for (int j = i+1; j < counter; ++j)
  {
    if (elemList[j*3+0] < 0) continue;

    if (elemList[j*3+0] == elemList[i*3+0])
    {
      bool iMatchesLE = (conf_face[i] == Face_id[i]);
      bool jMatchesLE = (conf_face[j] == Face_id[j]);

      if (!iMatchesLE && jMatchesLE)
      {
        // j is better, so discard i and replace it with j
        elemList[i*3+0] = elemList[j*3+0];
        elemList[i*3+1] = elemList[j*3+1];
        elemList[i*3+2] = elemList[j*3+2];
        conf_face[i]    = conf_face[j];
        hang_v_id[i]    = hang_v_id[j];
        Face_id[i]      = Face_id[j];

        // discard j
        elemList[j*3+0] = -1;
        elemList[j*3+1] = -1;
        elemList[j*3+2] = -1;
        conf_face[j]    = -1;
        hang_v_id[j]    = -1;
        Face_id[j]      = -1;
      }
      else
      {
        // keep i (either i already matches LE, or neither matches)
        elemList[j*3+0] = -1;
        elemList[j*3+1] = -1;
        elemList[j*3+2] = -1;
        conf_face[j]    = -1;
        hang_v_id[j]    = -1;
        Face_id[j]      = -1;
      }
    }
  }
}
// compact
dlong compactCounter = 0;

  memory<hlong> elemList2(counter*3,-1);
  memory<dlong> conf_face2(counter,-1);
  memory<hlong> hang_v_id2(counter,-1);
  memory<dlong> Face_id2(counter,-1);

  for (int i = 0; i < counter; ++i)
  {
    if (elemList[i*3+0] > -1)
    {
      elemList2[compactCounter*3+0] = elemList[i*3+0];
      elemList2[compactCounter*3+1] = elemList[i*3+1];
      elemList2[compactCounter*3+2] = elemList[i*3+2];
      conf_face2[compactCounter]    = conf_face[i];
      hang_v_id2[compactCounter]    = hang_v_id[i];
      Face_id2[compactCounter]   = Face_id[i];
      compactCounter++;
    }
  }

  counter = compactCounter;

  // sorting algorithm 
    bool swp;
  for (int i = 0; i < counter-1; ++i)
  {
    swp = false;
    for (int j = 0; j < counter-i-1; ++j)
    {
      if (elemList2[j*3+1] > elemList2[(j+1)*3+1] ||
         (elemList2[j*3+1] == elemList2[(j+1)*3+1] &&
          elemList2[j*3+2] > elemList2[(j+1)*3+2]))
      {
        std::swap(elemList2[j*3+0], elemList2[(j+1)*3+0]);
        std::swap(elemList2[j*3+1], elemList2[(j+1)*3+1]);
        std::swap(elemList2[j*3+2], elemList2[(j+1)*3+2]);

        std::swap(conf_face2[j],  conf_face2[j+1]);
        std::swap(hang_v_id2[j],  hang_v_id2[j+1]);
        std::swap(Face_id2[j], Face_id2[j+1]);

        swp = true;
      }
    }
    if (!swp) break;
  }

 dlong newCount = 0;

  for (int i = 0; i < counter; ++i)
  {
    if (i > 0)
    {
      bool sameEdge =
        (elemList2[i*3+1] == elemList2[(i-1)*3+1]) &&
        (elemList2[i*3+2] == elemList2[(i-1)*3+2]);

      if (sameEdge)
      {
        // inherit same node id as previous request
        new_v_id[i*2+0] = elemList2[i*3+0];
        new_v_id[i*2+1] = new_v_id[(i-1)*2+1];
        continue;
      }
    }

    // first row of a new edge group
    new_v_id[i*2+0] = elemList2[i*3+0];

    if (conf_face2[i] == Face_id2[i])
    {
      // hanging node already lies on split edge
      new_v_id[i*2+1] = hang_v_id2[i];
    }
    else
    {
      // create a new midpoint node for this edge group
      new_v_id[i*2+1] = mesh.Nnodes + newCount;
      newCount++;
    }
  }

  *new_vertex = newCount;

  printf("new_vertex=%lld,counter=%d\n", (long long)*new_vertex,counter);

  for (int i = 0; i < counter; ++i)
  {
    printf("e=%lld, vid0=%lld, vid1=%lld, conf_face=%d, split_face=%d, newid=%lld\n",
           (long long)elemList2[i*3+0],
           (long long)elemList2[i*3+1],
           (long long)elemList2[i*3+2],
           conf_face2[i],
           Face_id2[i],
           (long long)new_v_id[i*2+1]);
  }

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

void adaptivity_t::NewestVertex2(memory<dlong>& FaceFlag, 
                               memory<dlong>& RefFlag, 
                               dlong level, dlong Nrefine,
                               memory<hlong>& new_v_id,
                               hlong* new_vertex){

  dlong counter = 0;
  memory<hlong> elemList(Nrefine*3,-1);
  for (int e = 0; e < mesh.Nelements; ++e)
  {

      if (RefFlag[e]==1&& EToRefLevel[e]<level+3)
      {     
      dlong Face_id;
      hlong vid0; hlong vid1;
      const dlong id = e*mesh.Nverts;
      const dlong idf = e*mesh.Nfaces; 
       if (EToRefLevel[e]<1)
      {
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
      

      if (mag0 >= mag1 && mag0 >= mag2) {Face_id = idf+0; 
      vid0= mesh.EToV[idf+0]; vid1= mesh.EToV[idf+1];} 
      else if (mag1 >= mag2) {Face_id = idf+1; vid0= mesh.EToV[idf+1]; vid1= mesh.EToV[idf+2];} 
      else {Face_id = idf+2; vid0= mesh.EToV[idf+2]; vid1= mesh.EToV[idf+0];}
    }
     else
      {      
      Face_id = idf+0; 
      vid0= mesh.EToV[idf+0]; vid1= mesh.EToV[idf+1];      
      }    
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

}

} // namespace libp
