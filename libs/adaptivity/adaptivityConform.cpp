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
void adaptivity_t::Conform(memory<dlong>& RefFlag,                        
                         memory<dlong>& FaceFlag,
                         dlong Nrefine){

  
  dlong const MAX_REFINEMENT_LEVEL = 1;
  
  // Determine elements to be refined by using Refine Flag
  /*
  memory<dlong> Ref(Nrefine*2,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming
  dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==1)
    {
      Ref[ii] = e;
      ii = ii + 1;
    }
  }
  printf("ii=%d\n",ii );*/
  //int ii = Nrefine;
  // Conforming loop
  /*#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {
                
        // Find vertex locations of elements to be refined
        const dfloat x0 = mesh.EX[id+0]; const dfloat x1 = mesh.EX[id+1]; const dfloat x2 = mesh.EX[id+2];    
        const dfloat y0 = mesh.EY[id+0]; const dfloat y1 = mesh.EY[id+1]; const dfloat y2 = mesh.EY[id+2];

        // Find Longest Edge
        const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
        const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));

        dlong Face_id;
        if (mag0 >= mag1 && mag0 >= mag2) {
          Face_id = id+0;
        } else if (mag1 >= mag2) {
          Face_id = id+1;
        } else {
          Face_id = id+2;
        }

        if (Face_id==id+0)
        {
          hlong const ne_id = mesh.EToE[id+0];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+0]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             Nrefine++; 
             //Ref[ii+i] = neighbor_id;
                          //new_vertex--;
          }
        }

        if (Face_id==id+1)
        {         
          hlong const ne_id = mesh.EToE[id+1];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+1]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             Nrefine++; 
             //Ref[ii+i] = neighbor_id;
             //new_vertex--;
          }
        }

        if (Face_id==id+2)
        {          
          hlong const ne_id = mesh.EToE[id+2];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+2]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1; 
             Nrefine++;
             //Ref[ii+i] = neighbor_id;
             //new_vertex--;

          }        
        } 
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }             
  }
 */
  
    #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {

        if (FaceFlag[id+0]==1)
        {
          hlong const ne_id = mesh.EToE[id+0];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+0]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[id+0];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
                          //new_vertex--;
          }
        }

        if (FaceFlag[id+1]==1)
        {         
          hlong const ne_id = mesh.EToE[id+1];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+1]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[id+1];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
             //new_vertex--;
          }
        }

        if (FaceFlag[id+2]==1)
        {          
          hlong const ne_id = mesh.EToE[id+2];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+2]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[id+2];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1; 
             Nrefine++;
             
             //new_vertex--;

          }        
        } 
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }             
  }

  /*#pragma omp parallel for
  for (int i = 0; i < ii; ++i)
  {
      int e = Ref[i];
      const hlong id = e*mesh.Nverts; 
      if (RefFlag[e]==1 && EToRefLevel[e]==MAX_REFINEMENT_LEVEL)
      {
        
        hlong sib_e ;
        hlong sib0_id;
        hlong sib_id;
        hlong v0 ;
        hlong v1 ;
        hlong v2 ;
        hlong v0_sib;
        hlong v1_sib;
        hlong v2_sib;     
  
        // Extract Vertex Number of Elements to be Combined
        sib0_id = PToC[2*e+0]*mesh.Nverts; 
        sib_e = PToC[2*e+1];
        sib_id = sib_e*mesh.Nverts;
        v0 = mesh.EToV[sib0_id+0]; 
        v1 = mesh.EToV[sib0_id+1]; 
        v2 = mesh.EToV[sib0_id+2];
        v0_sib = mesh.EToV[sib_id+0]; 
        v1_sib = mesh.EToV[sib_id+1]; 
        v2_sib = mesh.EToV[sib_id+2];
        

        
        hlong id_del; // Id of the node that will be unnecessary
        if (v0 >= v1 && v0 >= v2) {
        id_del = sib0_id+0;
        } else if (v1 >= v2) {
        id_del = sib0_id+1;
        } else {
        id_del = sib0_id+2;
        }

        hlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
        if (v0_sib >= v1_sib && v0_sib >= v2_sib) {
        id_del2 = sib_id+0;
        } else if (v1_sib >= v2_sib) {
        id_del2 = sib_id+1;
        } else {
        id_del2 = sib_id+2;
        }
        
        if (id_del==sib0_id+1&&id_del2==sib_id+0)
        {
          hlong const neighbor_id = mesh.EToE[sib0_id+0];
           
          if(RefFlag[neighbor_id]!=1 && mesh.EToE[sib0_id+0]!=-1){
             RefFlag[neighbor_id]=1;
             Nrefine++; 
             printf("in the conforming loop,neighbor_id=%d\n",neighbor_id);
             
          }
        }

        if (id_del==sib0_id+2&&id_del2==sib_id+1)
        {         
          hlong const neighbor_id = mesh.EToE[sib0_id+1];
          //printf("mesh.EToE[id+1]=%d\n",mesh.EToE[id+1]);
        
          if(RefFlag[neighbor_id]!=1 && mesh.EToE[id+1]!=-1){

             RefFlag[neighbor_id]=1;
             Nrefine++; 
             
          }
        }

        if (id_del==sib0_id+2&&id_del2==sib_id+0)
        {          
          hlong const neighbor_id = mesh.EToE[sib0_id+2];
     
          if(RefFlag[neighbor_id]!=1 && mesh.EToE[id+2]!=-1){
             RefFlag[neighbor_id]=1; 
             Nrefine++;
             
          }        
        } 

      }         
  }

  Refine(q,RefFlag,Nrefine);*/

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}
}