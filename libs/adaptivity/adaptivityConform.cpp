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
void adaptivity_t::Conform(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong& Nrefine){ 
dlong const MAX_REFINEMENT_LEVEL = 2;
  
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
          //printf("f1=%d,f2=%lld,f3=%d,\n",RefFlag[ne_id],mesh.EToE[id+1],EToRefLevel[ne_id]);
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

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}
  
void adaptivity_t::ConformLE(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong& Nrefine){ 
dlong const MAX_REFINEMENT_LEVEL = 2;
  
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      if (RefFlag[e]==1)
      {

        if (FaceFlag[id+0]==1)
        {

          hlong const ne_id = mesh.EToE[id+0];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+0]!=-1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id+2])
            {
             RefFlag[ne_id]=1;          
             Nrefine++;
            }
 
             
                          //new_vertex--;
          }
        }

        if (FaceFlag[id+1]==1)
        {         
          hlong const ne_id = mesh.EToE[id+1];
          //printf("f1=%d,f2=%lld,f3=%d,\n",RefFlag[ne_id],mesh.EToE[id+1],EToRefLevel[ne_id]);
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+1]!=-1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id+2])
            {
             RefFlag[ne_id]=1;          
             Nrefine++;
            }
          }
        }

        if (FaceFlag[id+2]==1)
        {          
          hlong const ne_id = mesh.EToE[id+2];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+2]!=-1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id+2])
            {
             RefFlag[ne_id]=1;          
             Nrefine++;
            }

          }        
        } 
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }             
  }

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}
}