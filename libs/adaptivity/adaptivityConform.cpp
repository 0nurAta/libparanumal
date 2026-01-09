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

void adaptivity_t::ConformByBisect(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 1;
printf("Conforming Start!, Nrefine=%d\n",Nrefine);
hlong nn = 0;
hlong new_vertex = 0;
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
             printf("conform_e0=%d\n",e );    
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+0];    
             FaceFlag[idf+0] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
             if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL)
             {
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
              *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
             }

        }

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        {    
             printf("conform_e1=%d\n",e );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
                          if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL)
             {
             BisectLocal1(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
             }
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {    
             printf("conform_e2=%d\n",e );      
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+2];    
             FaceFlag[idf+2] = 1;
             new_v_id[ef] = mesh.EToV[id+2];
             Nrefine++;
                          if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL)
             {
             BisectLocal2(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
        }
        }  

        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      }




  }

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}


void adaptivity_t::ConformByBisectMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 5;
printf("Conforming Start!, Nrefine=%d, Elements=%d\n",Nrefine,mesh.Nelements);
hlong nn = 0;
hlong new_vertex = 0;
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]>2)
      {
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2]; 
      const dlong idf = ConfFlag[e]*mesh.Nfaces;
      const dlong sib_id = sib_e*mesh.Nfaces; 
      printf("conform=%d sib_e=%d,neigh=%d\n",ef,sib_e,e );     

       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]))<1e-5) {   
             RefFlag[ef]=1;  
             FaceFlag[idf+0] = 1;
             new_v_id[idf+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e0=%d,new_id=%d\n",e,ef,new_v_id[idf+0] ); 
             
             Nrefine++;
             }
        } else if (mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
          if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+0] = 1;
             new_v_id[sib_id+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[sib_id+0] ); 
             Nrefine++;
             }}

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        
        {    
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[idf+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             } else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+1] = 1;
             new_v_id[sib_id+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }
                          
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {  
            if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<1e-5) {   
             //printf("e=%d,conform_e2=%d\n",e,ef );      
             RefFlag[ef]=1;  
             FaceFlag[idf+2] = 1;
             new_v_id[idf+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+2] );
             Nrefine++;
             }
              else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+2] = 1;
             new_v_id[sib_id+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }           
        
        }  
}
}
for (int ef = 0; ef < mesh.Nelements; ++ef)
  {
    const dlong idf = ef*mesh.Nfaces;

            // RULE 0
            if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+0] == 1)
             {
              printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
              FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;

               if(FaceFlag[idf+1] == 1){
               dlong ef_second = PCS[ef*3+2];
               printf("Rule 01: ef=%d,nn=%d\n",ef_second,nn );
               new_v_id[ef_second*mesh.Nfaces+0] = new_v_id[idf+1];
               FaceFlag[ef_second*mesh.Nfaces+0] = FaceFlag[idf+1];
               FaceFlag[idf+1]=0;
               BisectLocal0(ef_second,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);              
               FaceFlag[ef_second*mesh.Nfaces+1] = 0;
               *NN =nn;
               *New_vertex=new_vertex; 
               }

               if(FaceFlag[idf+2] == 1){
               RefFlag[ef]=1; 
               FaceFlag[idf+0]=1;
               new_v_id[idf+0] = new_v_id[idf+2];
               FaceFlag[idf+2] = 0;
               printf("Rule 02: ef=%d,nn=%d\n",ef,nn );
               BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
               FaceFlag[idf+0] = 0;
               
               *NN =nn;
               *New_vertex=new_vertex; 
               }
             }

             // RULE 1
             else if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+1] == 1)
             {
              printf("Rule 1: ef=%d,nn=%d\n",ef,nn );
             BisectLocal1(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
             FaceFlag[idf+1]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;

               if(FaceFlag[idf+2] == 1){
               
               dlong ef_second = PCS[ef*3+2];
               new_v_id[ef_second*mesh.Nfaces+0] = new_v_id[idf+2];
               FaceFlag[ef_second*mesh.Nfaces+0] = FaceFlag[idf+2];
               new_v_id[idf+2]=0;
               FaceFlag[idf+2]=0;   
               printf("Rule 12: ef=%d,nn=%d\n",ef_second,nn );             
               BisectLocal0(ef_second,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
               FaceFlag[ef_second*mesh.Nfaces+0] = 0;

               *NN =nn;
               *New_vertex=new_vertex; 
               }
             }

             // RULE 2
              else if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+2] == 1)
             {
              printf("Rule 2: ef=%d,nn=%d\n",ef,nn );
             BisectLocal2(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
             FaceFlag[idf+2]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
               }
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      } 

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

void adaptivity_t::ConformByCoarse(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 5;
printf("Conforming Start!, Nrefine=%d, Elements=%d\n",Nrefine,mesh.Nelements);
hlong nn = 0;
hlong new_vertex = 0;
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]>2)
      {
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2]; 
      const dlong idf = ConfFlag[e]*mesh.Nfaces;
      const dlong sib_id = sib_e*mesh.Nfaces; 
      //printf("conform=%d sib_e=%d\n",ef,sib_e );     

       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]))<1e-5) {   
             RefFlag[ef]=1;  
             FaceFlag[idf+0] = 1;
             new_v_id[idf+0] = mesh.EToV[id+2];
             //printf("e=%d,conform_e0=%d,new_id=%d\n",e,ef,new_v_id[idf+0] ); 
             
             Nrefine++;
             }
        } else if (mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
          if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_id]=1;  
             FaceFlag[sib_id+0] = 1;
             new_v_id[sib_id+0] = mesh.EToV[id+2];
             //printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[sib_id+0] ); 
             Nrefine++;
             }}

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        
        {    
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[idf+1] = mesh.EToV[id+2];
             //printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             } else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_id]=1;  
             FaceFlag[sib_id+1] = 1;
             new_v_id[sib_id+1] = mesh.EToV[id+2];
             //printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }
                          
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {  
            if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<1e-5) {   
             //printf("e=%d,conform_e2=%d\n",e,ef );      
             RefFlag[ef]=1;  
             FaceFlag[idf+2] = 1;
             new_v_id[idf+2] = mesh.EToV[id+2];
             //printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+2] );
             Nrefine++;
             }
              else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_id]=1;  
             FaceFlag[sib_id+2] = 1;
             new_v_id[sib_id+2] = mesh.EToV[id+2];
             //printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }           
        
        }  
}
}
for (int ef = 0; ef < mesh.Nelements; ++ef)
  {
    const dlong idf = ef*mesh.Nfaces;

            // RULE 0
            if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+0] == 1)
             {
              //printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
              FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;

               if(FaceFlag[idf+1] == 1){
               dlong ef_second = PCS[ef*3+2];
               //printf("Rule 01: ef=%d,nn=%d\n",ef_second,nn );
               new_v_id[ef_second*mesh.Nfaces+0] = new_v_id[idf+1];
               FaceFlag[ef_second*mesh.Nfaces+0] = FaceFlag[idf+1];
               FaceFlag[idf+1]=0;
               BisectLocal0(ef_second,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);              
               FaceFlag[ef_second*mesh.Nfaces+1] = 0;
               *NN =nn;
               *New_vertex=new_vertex; 
               }

               if(FaceFlag[idf+2] == 1){
               RefFlag[ef]=1; 
               FaceFlag[idf+0]=1;
               new_v_id[idf+0] = new_v_id[idf+2];
               FaceFlag[idf+2] = 0;
               //printf("Rule 02: ef=%d,nn=%d\n",ef,nn );
               BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
               FaceFlag[idf+0] = 0;
               
               *NN =nn;
               *New_vertex=new_vertex; 
               }
             }

             // RULE 1
             if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+1] == 1)
             {
              //printf("Rule 1: ef=%d,nn=%d\n",ef,nn );
             BisectLocal1(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
             FaceFlag[idf+1]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;

               if(FaceFlag[idf+2] == 1){
               
               dlong ef_second = PCS[ef*3+2];
               new_v_id[ef_second*mesh.Nfaces+0] = new_v_id[idf+2];
               FaceFlag[ef_second*mesh.Nfaces+0] = FaceFlag[idf+2];
               new_v_id[idf+2]=0;
               FaceFlag[idf+2]=0;   
               //printf("Rule 12: ef=%d,nn=%d\n",ef_second,nn );             
               BisectLocal0(ef_second,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
               FaceFlag[ef_second*mesh.Nfaces+0] = 0;

               *NN =nn;
               *New_vertex=new_vertex; 
               }
             }

             // RULE 2
              if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+2] == 1)
             {
              //printf("Rule 2: ef=%d,nn=%d\n",ef,nn );
             BisectLocal2(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
             FaceFlag[idf+2]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
               }
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      } 

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}


void adaptivity_t::ConformByLEMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 5;
printf("Conforming Start!, Nrefine=%d, Elements=%d\n",Nrefine,mesh.Nelements);
hlong nn = 0;
hlong new_vertex = 0;
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]>2)
      {
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2]; 
      const dlong idf = ConfFlag[e]*mesh.Nfaces;
      const dlong sib_id = sib_e*mesh.Nfaces; 
      //printf("conform=%d sib_e=%d\n",ef,sib_e );     

       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]))<1e-5) {   
             RefFlag[ef]=1;  
             FaceFlag[idf+0] = 1;
             new_v_id[idf+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e0=%d,new_id=%d\n",e,ef,new_v_id[idf+0] ); 
             
             Nrefine++;
             }
        } else if (mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
          if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+0] = 1;
             new_v_id[sib_id+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[sib_id+0] ); 
             Nrefine++;
             }}

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        
        {    
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[idf+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             } else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+1] = 1;
             new_v_id[sib_id+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }
                          
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {  
            if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<1e-5) {   
             //printf("e=%d,conform_e2=%d\n",e,ef );      
             RefFlag[ef]=1;  
             FaceFlag[idf+2] = 1;
             new_v_id[idf+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+2] );
             Nrefine++;
             }
              else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+2] = 1;
             new_v_id[sib_id+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }           
        
        }  
}
}
for (int ef = 0; ef < mesh.Nelements; ++ef)
  {
    const dlong idf = ef*mesh.Nfaces;
    const dlong id = ef*mesh.Nverts;
            // RULE 0
            if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  (FaceFlag[idf+0] == 1 || FaceFlag[idf+1] == 1 || FaceFlag[idf+2] == 1))
             {
              
              FaceFlag[idf+0] = 0;
              FaceFlag[idf+1] = 0;
              FaceFlag[idf+2] = 0;

              // Find vertex locations of element to be refined
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

              if (FaceFlag[idf+0]==1)
              {
              hlong vA = mesh.EToV[id + 0]; 
              hlong vB = mesh.EToV[id + 1];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+0] = mesh.Nnodes + (vmin * 31 + vmax); 
              printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }
              
              if (FaceFlag[idf+1]==1)
              {
              hlong vA = mesh.EToV[id + 1]; 
              hlong vB = mesh.EToV[id + 2];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+1] = mesh.Nnodes + (vmin * 31 + vmax); 
              printf("Rule 1: ef=%d,nn=%d\n",ef,nn );
              BisectLocal1(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+1]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }

              if (FaceFlag[idf+2]==1)
              {
              hlong vA = mesh.EToV[id + 2]; 
              hlong vB = mesh.EToV[id + 0];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+2] = mesh.Nnodes + (vmin * 31 + vmax); 
              printf("Rule 2: ef=%d,nn=%d\n",ef,nn );
              BisectLocal2(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+2]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }
             }
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      } 

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

void adaptivity_t::ConformByNVBMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 5;
printf("Conforming Start!, Nrefine=%d, Elements=%d\n",Nrefine,mesh.Nelements);
hlong nn = 0;
hlong new_vertex = 0;
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]>2)
      {
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2]; 
      const dlong idf = ConfFlag[e]*mesh.Nfaces;
      const dlong sib_id = sib_e*mesh.Nfaces; 
      //printf("conform=%d sib_e=%d\n",ef,sib_e );     

       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]))<1e-5) {   
             RefFlag[ef]=1;  
             FaceFlag[idf+0] = 1;
             new_v_id[idf+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e0=%d,new_id=%d\n",e,ef,new_v_id[idf+0] ); 
             
             Nrefine++;
             }
        } else if (mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
          if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+0] = 1;
             new_v_id[sib_id+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[sib_id+0] ); 
             Nrefine++;
             }}

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        
        {    
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[idf+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             } else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+1] = 1;
             new_v_id[sib_id+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }
                          
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {  
            if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<1e-5) {   
             //printf("e=%d,conform_e2=%d\n",e,ef );      
             RefFlag[ef]=1;  
             FaceFlag[idf+2] = 1;
             new_v_id[idf+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+2] );
             Nrefine++;
             }
              else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+2] = 1;
             new_v_id[sib_id+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }           
        
        }  
}
}
for (int ef = 0; ef < mesh.Nelements; ++ef)
  {
    const dlong idf = ef*mesh.Nfaces;
    const dlong id = ef*mesh.Nverts;
                // RULE 0
            if (EToRefLevel[ef]<1 &&  (FaceFlag[idf+0] == 1 || FaceFlag[idf+1] == 1 || FaceFlag[idf+2] == 1))
             {
              
              FaceFlag[idf+0] = 0;
              FaceFlag[idf+1] = 0;
              FaceFlag[idf+2] = 0;

              // Find vertex locations of element to be refined
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

              if (FaceFlag[idf+0]==1)
              {
              hlong vA = mesh.EToV[id + 0]; 
              hlong vB = mesh.EToV[id + 1];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+0] = mesh.Nnodes + (vmin * 10000 + vmax); 
              printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }
              
              if (FaceFlag[idf+1]==1)
              {
              hlong vA = mesh.EToV[id + 1]; 
              hlong vB = mesh.EToV[id + 2];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+1] = mesh.Nnodes + (vmin * 10000 + vmax); 
              printf("Rule 1: ef=%d,nn=%d\n",ef,nn );
              BisectLocal1(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+1]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }

              if (FaceFlag[idf+2]==1)
              {
              hlong vA = mesh.EToV[id + 2]; 
              hlong vB = mesh.EToV[id + 0];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+2] = mesh.Nnodes + (vmin * 10000 + vmax); 
              printf("Rule 2: ef=%d,nn=%d\n",ef,nn );
              BisectLocal2(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+2]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              }
             }
            // RULE 0
            else if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  (FaceFlag[idf+0] == 1 || FaceFlag[idf+1] == 1 || FaceFlag[idf+2] == 1))
             {
              
              FaceFlag[idf+0] = 0;
              FaceFlag[idf+1] = 0;
              FaceFlag[idf+2] = 0;           
              
              hlong vA = mesh.EToV[id + 0]; 
              hlong vB = mesh.EToV[id + 1];

              hlong vmin = (vA < vB) ? vA : vB;
              hlong vmax = (vA > vB) ? vA : vB;
           
              new_v_id[id+0] = mesh.Nnodes + (vmin * 10000 + vmax); 
              printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              BisectLocal0(ef,Q,Qold,RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,RefLevel);
                FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
              

             }
        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      } 

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

void adaptivity_t::ConformByRedMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                                   memory<dlong>& new_v_id, 
                                                                                   dlong& Nrefine, 
                                                                                   memory<dfloat>& EX_new, 
                                                                                   memory<dfloat>& EY_new,
                                                                                   memory<hlong>& EToV_new,
                                                                                   memory<int>& EToB_new,
                                                                                   memory<dlong>& SplitFlag,                                                                           
                                                                                   hlong* NN,
                                                                                   hlong* New_vertex,
                                                                                   dlong RefLevel){ 
dlong const MAX_REFINEMENT_LEVEL = 1;
printf("Conforming Start!, Nrefine=%d, Elements=%d\n",Nrefine,mesh.Nelements);
hlong nn = 0;
hlong new_vertex = 0;
    //#pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      //int e = Ref[i];
      const dlong id = e*mesh.Nfaces; 
      
      //printf("e=%d,conform_e=%d\n",e, ConfFlag[e]);  
      if (ConfFlag[e]>2)
      {
      const dlong ef = ConfFlag[e];  
      const dlong sib_e = PCS[ef*3+2]; 
      const dlong idf = ConfFlag[e]*mesh.Nfaces;
      const dlong sib_id = sib_e*mesh.Nfaces; 
      //printf("conform=%d sib_e=%d\n",ef,sib_e );     

       // hlong const ne_id_0 = mesh.EToE[idf+0];
       // hlong const ne_id_1 = mesh.EToE[idf+1];
       // hlong const ne_id_2 = mesh.EToE[idf+2];

        if (mesh.EToF[idf+0]==-1 && mesh.EToB[idf+0]==-1)
        {
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+0]+mesh.EX[idf+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+0]+mesh.EY[idf+1]))<1e-5) {   
             RefFlag[ef]=1;  
             FaceFlag[idf+0] = 1;
             new_v_id[idf+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e0=%d,new_id=%d\n",e,ef,new_v_id[idf+0] ); 
             
             Nrefine++;
             }
        } else if (mesh.EToF[sib_id+0]==-1 && mesh.EToB[sib_id+0]==-1){
          if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+0]+mesh.EX[sib_id+1]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+0]+mesh.EY[sib_id+1]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+0] = 1;
             new_v_id[sib_id+0] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[sib_id+0] ); 
             Nrefine++;
             }}

                if (mesh.EToF[idf+1]==-1 && mesh.EToB[idf+1]==-1)
        
        {    
             if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+1]+mesh.EX[idf+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+1]+mesh.EY[idf+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[ef]=1;
             hlong const fN = mesh.EToF[idf+1];    
             FaceFlag[idf+1] = 1;
             new_v_id[idf+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             } else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+1]+mesh.EX[sib_id+2]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+1]+mesh.EY[sib_id+2]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+1] = 1;
             new_v_id[sib_id+1] = mesh.EToV[id+2];
             printf("e=%d,conform_e1=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }
                          
        }

                if (mesh.EToF[idf+2]==-1 && mesh.EToB[idf+2]==-1)
        {  
            if(abs(mesh.EX[id+2]-0.5*(mesh.EX[idf+2]+mesh.EX[idf+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[idf+2]+mesh.EY[idf+0]))<1e-5) {   
             //printf("e=%d,conform_e2=%d\n",e,ef );      
             RefFlag[ef]=1;  
             FaceFlag[idf+2] = 1;
             new_v_id[idf+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+2] );
             Nrefine++;
             }
              else if (abs(mesh.EX[id+2]-0.5*(mesh.EX[sib_id+2]+mesh.EX[sib_id+0]))<1e-5&&abs(mesh.EY[id+2]-0.5*(mesh.EY[sib_id+2]+mesh.EY[sib_id+0]))<1e-5) { 
             //printf("e=%d,conform_e1=%d\n",e,ef );     
             RefFlag[sib_e]=1;  
             FaceFlag[sib_id+2] = 1;
             new_v_id[sib_id+2] = mesh.EToV[id+2];
             printf("e=%d,conform_e2=%d,new_id=%d\n",e,ef,new_v_id[idf+1] ); 
             Nrefine++;
             }           
        
        } 

}
}
for (int ef = 0; ef < mesh.Nelements; ++ef)
  {
    const dlong idf = ef*mesh.Nfaces;

            // RULE 0
            if (EToRefLevel[ef]<MAX_REFINEMENT_LEVEL&&  FaceFlag[idf+0] == 1&&  FaceFlag[idf+1] == 1&&  FaceFlag[idf+2] == 1)
             {
              printf("Rule 0: ef=%d,nn=%d\n",ef,nn );
              RedLocal(ef,Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,&nn,&new_vertex,1);
              FaceFlag[idf+0]=0;
               *NN =nn;
               *New_vertex=new_vertex;
                          //new_vertex--;
             }

        //printf("new_vertex_count_conf=%lld\n",new_vertex);
      } 

printf("Conforming Done!, Nrefine=%d\n",Nrefine);        
}

}