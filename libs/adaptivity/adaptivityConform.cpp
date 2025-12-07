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
dlong const MAX_REFINEMENT_LEVEL = 1;
  
    #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      const dlong idf = e*mesh.Nfaces;
      
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {

        if (FaceFlag[id+0]==1)
        {

          hlong const ne_id = mesh.EToE[idf+0];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+0]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+0];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
                          //new_vertex--;
          }
        }

        if (FaceFlag[id+1]==1)
        {         
          hlong const ne_id = mesh.EToE[idf+1];
          //printf("f1=%d,f2=%lld,f3=%d,\n",RefFlag[ne_id],mesh.EToE[id+1],EToRefLevel[ne_id]);
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+1]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+1];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
             //new_vertex--;
          }
        }

        if (FaceFlag[id+2]==1)
        {          
          hlong const ne_id = mesh.EToE[idf+2];
          
          if(RefFlag[ne_id]!=1 && mesh.EToE[id+2]!=-1 && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+2];
             
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
  for (dlong e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces;

      if (RefFlag[e]==1)
      {

        if (FaceFlag[idf+0]==1)
        {

          hlong const ne_id = mesh.EToE[idf+0];
          
          if(ne_id!=-1 && RefFlag[ne_id]!=1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id*mesh.Nfaces+2])
            {
             RefFlag[ne_id]=1;
             Nrefine++;
            }
                          //new_vertex--;
          }
        }

        if (FaceFlag[idf+1]==1)
        {         
          hlong const ne_id = mesh.EToE[idf+1];
          //printf("f1=%d,f2=%lld,f3=%d,\n",RefFlag[ne_id],mesh.EToE[id+1],EToRefLevel[ne_id]);
          if(ne_id!=-1 && RefFlag[ne_id]!=1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id*mesh.Nfaces+2])
            {
             RefFlag[ne_id]=1;          
             Nrefine++;
            }
          }
        }

        if (FaceFlag[idf+2]==1)
        {          
          hlong const ne_id = mesh.EToE[idf+2];
          
          if(ne_id!=-1 && RefFlag[ne_id]!=1){
            if (mesh.EToV[id+2]!=mesh.EToV[ne_id*mesh.Nfaces+2])
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

void adaptivity_t::ConformByVertex(memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag, dlong& Nrefine){ 
dlong const MAX_REFINEMENT_LEVEL = 2;
printf("Conforming Start!, Nrefine=%d\n",Nrefine); 
    #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      const dlong idf = e*mesh.Nfaces;
      
      if (ConfFlag[e]==1)
      {
        
       hlong newest_vertex = mesh.EToV[id+2];

        if (FaceFlag[id+0]==1)
        {
          
          hlong const ne_id = mesh.EToE[idf+0];
         hlong newest_vertex_ne = mesh.EToV[idf+2];
         printf("e= %d, ne_id= %d, ConfFlag=%d,  RefFlag=%d,FaceFlag[id+0]=%d\n",e,ne_id,ConfFlag[e],RefFlag[e],FaceFlag[id+0]); 
         printf("ne_id=%d\n",EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL);
          if(ne_id!=-1 && RefFlag[ne_id]!=1 && abs(mesh.EX[newest_vertex_ne]-mesh.EX[newest_vertex])<1e-5 
                                            && abs(mesh.EY[newest_vertex_ne]-mesh.EY[newest_vertex])<1e-5
                                            && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+0];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
                          //new_vertex--;
          }
        }

        else if (FaceFlag[id+1]==1)
        { 
                 
          hlong const ne_id = mesh.EToE[idf+1];
          hlong newest_vertex_ne = mesh.EToV[idf+2];
          printf("e= %d, ne_id= %d, ConfFlag=%d,  RefFlag=%d,FaceFlag[id+1]=%d\n",e,ne_id,ConfFlag[e],RefFlag[e],FaceFlag[id+1] ); 
          printf("ne_id=%d\n",EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL);

          if(ne_id!=-1 && RefFlag[ne_id]!=1 && abs(mesh.EX[newest_vertex_ne]-mesh.EX[newest_vertex])<1e-5 
                                            && abs(mesh.EY[newest_vertex_ne]-mesh.EY[newest_vertex])<1e-5
                                            && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){

             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+1];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1;
             Nrefine++; 
             
             //new_vertex--;
          }
        }

        else if (FaceFlag[id+2]==1)
        { 
                  
          hlong const ne_id = mesh.EToE[idf+2];
          hlong newest_vertex_ne = mesh.EToV[idf+2];
          printf("e= %d, ne_id= %d, ConfFlag=%d,  RefFlag=%d,FaceFlag[id+2]=%d\n",e,ne_id,ConfFlag[e],RefFlag[e],FaceFlag[id+2]); 
          if(ne_id!=-1 && RefFlag[ne_id]!=1 && abs(mesh.EX[newest_vertex_ne]-mesh.EX[newest_vertex])<1e-5 
                                            && abs(mesh.EY[newest_vertex_ne]-mesh.EY[newest_vertex])<1e-5
                                            && EToRefLevel[ne_id]<MAX_REFINEMENT_LEVEL){
             RefFlag[ne_id]=1;
             hlong const fN = mesh.EToF[idf+2];
             
             FaceFlag[ne_id*mesh.Nfaces+fN] = 1; 
             Nrefine++;
             
             //new_vertex--;

          }
        } 
        else  ConfFlag[e]= 0; 
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }             
  }

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

void adaptivity_t::ConformByID(memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,memory<dlong>& new_v_id, dlong& Nrefine){ 
dlong const MAX_REFINEMENT_LEVEL = 1;
printf("Conforming Start!, Nrefine=%d\n",Nrefine); 
    #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]!=0)
      {
          
      const dlong ef = ConfFlag[e];   
      const dlong idf = ConfFlag[e]*mesh.Nfaces;  
       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             printf("conform_e=%d\n",e );    
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+0];    
             FaceFlag[idf+0] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
                          //new_vertex--;
        }

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        {    
             printf("conform_e=%d\n",e );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
                          //new_vertex--;
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {    
             printf("conform_e=%d\n",e );      
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+2];    
             FaceFlag[idf+2] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
                          //new_vertex--;
        }      
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }             
  }

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

}