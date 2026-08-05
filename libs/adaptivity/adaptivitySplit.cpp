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


// Bisection
namespace libp {

void adaptivity_t::BisectLocal0(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
  
      
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 

      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces

          // Uniquely number new vertex 
          //hlong Local_id = 0+idf+mesh.Nnodes;

          // To have unique global vertex number
          //hlong Neigh_id = mesh.EToF[idf+0]+mesh.Nfaces*mesh.EToE[idf+0]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+0];
          printf("newNode=%lld\n",newNode );
          // 1st child vertex ids
          EToV_new[id+0] = v2;
          EToV_new[id+1] = v0;
          EToV_new[id+2] = newNode;
          
          // 2nd child vertex ids
          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          // 2nd child vertex coordinates
          EX_new[id_new+0] = x1;
          EX_new[id_new+1] = x2;
          EX_new[id_new+2] = 0.5f*(x0+x1);

          // 1st child vertex coordinates
          EX_new[id+0] = x2;
          EX_new[id+1] = x0;
          EX_new[id+2] = 0.5f*(x0+x1);

          // 2nd child vertex coordinates
          EY_new[id_new+0] = y1;
          EY_new[id_new+1] = y2;
          EY_new[id_new+2] = 0.5f*(y0+y1);

          // 1st child vertex coordinates
          EY_new[id+0] = y2;
          EY_new[id+1] = y0;
          EY_new[id+2] = 0.5f*(y0+y1);
          

          // Boundary information

          // 1st child
          EToB_new[id+0] = b2; 
          EToB_new[id+1] = b0; // New face will be inner in any situation
          EToB_new[id+2] = -1;
          // 2nd child
          EToB_new[id_new+0] = b1;
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = b0; // New face will be inner in any situation

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];
          
          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)+0]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 2;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
          const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              //const dlong id1_int = IntFlag[e*(level+3)+EToRefLevel[e]-1]-1;
              const dlong id1_int = 0;
              //const dlong id2_int = IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1]-1;
              const dlong id2_int = 1;
              printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
            //    printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

            for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
            //printf("Q_inbisect=%f\n",Q[1621] );
          //
          //ConfFlag[e]=mesh.EToE[idf+0];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+0];   
          RefFlag[e]=0;
          nn++;
          nv++;
       
      *new_vertex = nv;
      *NN = nn;      
}

void adaptivity_t::BisectLocal1(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
  
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 


      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
          //hlong Local_id = 1+idf+mesh.Nnodes;

          //hlong Neigh_id = mesh.EToF[idf+1]+mesh.Nfaces*mesh.EToE[idf+1]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+1];
          printf("newNode=%lld\n",newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v2;
          EToV_new[id_new+1] = v0;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(x1+x2); 

          
          EX_new[id_new+0] = x2;
          EX_new[id_new+1] = x0; 
          EX_new[id_new+2] = 0.5f*(x1+x2);  

          EY_new[id+2] = 0.5f*(y1+y2); 

           
          EY_new[id_new+0] = y2;
          EY_new[id_new+1] = y0; 
          EY_new[id_new+2] = 0.5f*(y1+y2); 

          EToB_new[id+0] = b0;
          EToB_new[id+1] = b1;
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = b2;
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = b1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 5;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 4;
              const dlong id2_int = 5;
                 //printf("id1_int=%d,id2_int=%d\n", id1_int,id2_int);
              //printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
              //  printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

                        for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          RefFlag[e]=0;
          //ConfFlag[e]=mesh.EToE[idf+1];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+1];   
          nn++;
          nv++;            
     *new_vertex = nv;
     *NN = nn;      
}

void adaptivity_t::BisectLocal2(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 
      
      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
          //hlong Local_id = 2+idf+mesh.Nnodes;
          
          //hlong Neigh_id = mesh.EToF[idf+2]+mesh.Nfaces*mesh.EToE[idf+2]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+2];
          printf("newNode=%lld\n",newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(x2+x0);

          EX_new[id_new+0] = x1;
          EX_new[id_new+1] = x2;
          EX_new[id_new+2] = 0.5f*(x2+x0);

          EY_new[id+2] = 0.5f*(y2+y0);

          EY_new[id_new+0] = y1;
          EY_new[id_new+1] = y2;
          EY_new[id_new+2] = 0.5f*(y2+y0);

          EToB_new[id+0] = b0;
          EToB_new[id+1] = -1;
          EToB_new[id+2] = b2;

          EToB_new[id_new+0] = b1;
          EToB_new[id_new+1] = b2;
          EToB_new[id_new+2] = -1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 2;
              const dlong id2_int = 3;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
                //printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          RefFlag[e]=0;
          //ConfFlag[e]=mesh.EToE[idf+2];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+2];   
          nn++;
          nv++;
      *new_vertex = nv;
      *NN = nn;      
}

void adaptivity_t::Bisect(  memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = 0 ; // new vertex
  hlong nn = 0 ; // Counts each refinement
  dlong const level = RefLevel;
  for (dlong e = 0; e < mesh.Nelements; ++e)
  {
      
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 

      if (RefFlag[e]==1 && EToRefLevel[e]<level)
      {
        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
        if (FaceFlag[idf+0]==1)
        { 

          // To have unique global vertex number
      
          //hlong  Local_id = 0*(level+3)+(mesh.Nfaces+3*(level+3))*e+mesh.Nnodes+EToRefLevel[e];      
          //hlong Neigh_id = mesh.EToF[idf+0]*(level+3)+(mesh.Nfaces+3*(level+3))*mesh.EToE[idf+0]+mesh.Nnodes+EToRefLevel[mesh.EToE[idf+0]];
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          printf("e0=%d,id_new0=%lld,mesh.EToF=%d\n",e,mesh.Nelements+nn,mesh.EToF[idf+0] );
          // Second version
          hlong vA = mesh.EToV[id + 0]; 
          hlong vB = mesh.EToV[id + 1];

          hlong vmin = (vA < vB) ? vA : vB;
          hlong vmax = (vA > vB) ? vA : vB;

          hlong newNode = mesh.Nnodes + (vmin * 31 + vmax);
          //hlong newNode = mesh.Nnodes + ( (vA + vB) * (vA + vB + 1) / 2 + vB );
          //hlong newNode = v0+v1+mesh.Nnodes;
          //printf("Local_id=%lld,Neigh_id=%lld,newNode=%lld\n",Local_id,Neigh_id,newNode );
          printf("newNode=%lld,vA=%lld,vB=%lld\n",newNode,vA,vB );
          // 1st child vertex ids
          EToV_new[id+0] = mesh.EToV[id+2];
          EToV_new[id+1] = mesh.EToV[id+0];
          EToV_new[id+2] = newNode;
          
          // 2nd child vertex ids
          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          // 2nd child vertex coordinates
          EX_new[id_new+0] = mesh.EX[id+1];
          EX_new[id_new+1] = mesh.EX[id+2];
          EX_new[id_new+2] = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);

          // 1st child vertex coordinates
          EX_new[id+0] = mesh.EX[id+2];
          EX_new[id+1] = mesh.EX[id+0];
          EX_new[id+2] = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);

          // 2nd child vertex coordinates
          EY_new[id_new+0] = mesh.EY[id+1];
          EY_new[id_new+1] = mesh.EY[id+2];
          EY_new[id_new+2] = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);

          // 1st child vertex coordinates
          EY_new[id+0] = mesh.EY[id+2];
          EY_new[id+1] = mesh.EY[id+0];
          EY_new[id+2] = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);
          

          // Boundary information

          // 1st child
          EToB_new[id+0] = mesh.EToB[id+2]; 
          EToB_new[id+1] = mesh.EToB[id+0]; // New face will be inner in any situation
          EToB_new[id+2] = -1;
          // 2nd child
          EToB_new[id_new+0] = mesh.EToB[id+1];
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+0]; // New face will be inner in any situation

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];
          

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)+0]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 2;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 0;
              const dlong id2_int = 1;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
                //printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          ConfFlag[e]=mesh.EToE[idf+0]; 
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+0]; 
          //printf("confFlag[e]=%d\n",mesh.EToE[idf+0] );
          RefFlag[e]=0;
          nn++;
          nv++;
        }

        else if (FaceFlag[idf+1]==1)
        {
          
          // To have unique global vertex number
          //hlong  Local_id = 1*(level+3)+(mesh.Nfaces+3*(level+3))*e+mesh.Nnodes+EToRefLevel[e];              
          //hlong Neigh_id = mesh.EToF[idf+1]*(level+3)+(mesh.Nfaces+3*(level+3))*mesh.EToE[idf+1]+mesh.Nnodes+EToRefLevel[mesh.EToE[idf+1]];
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          printf("e1=%d,id_new1=%lld,mesh.EToF=%d\n",e,mesh.Nelements+nn,mesh.EToF[idf+1] );
                    // Second version
          hlong vA = mesh.EToV[id + 1]; 
          hlong vB = mesh.EToV[id + 2];

          hlong vmin = (vA < vB) ? vA : vB;
          hlong vmax = (vA > vB) ? vA : vB;

          hlong newNode = mesh.Nnodes + (vmin * 31 + vmax);
          // hlong newNode = mesh.Nnodes + ( (vA + vB) * (vA + vB + 1) / 2 + vB );
           printf("newNode=%d,vA=%d,vB=%d\n",newNode,vA,vB );
          //hlong newNode = v1+v2+mesh.Nnodes;
          // printf("Neigh_e=%lld,MeshLevel=%d,NeighLevel=%d\n",mesh.EToE[idf+1],EToRefLevel[e],EToRefLevel[mesh.EToE[idf+1]]);
          //printf("Local_id=%lld,Neigh_id=%lld,newNode=%lld\n",Local_id,Neigh_id,newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v2;
          EToV_new[id_new+1] = v0;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]); 

          
          EX_new[id_new+0] = mesh.EX[id+2];
          EX_new[id_new+1] = mesh.EX[id+0]; 
          EX_new[id_new+2] = EX_new[id+2];  

          EY_new[id+2] = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]); 

           
          EY_new[id_new+0] = mesh.EY[id+2];
          EY_new[id_new+1] = mesh.EY[id+0]; 
          EY_new[id_new+2] = EY_new[id+2]; 

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[id+1];
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = mesh.EToB[id+2];
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+1];

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 5;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 4;
              const dlong id2_int = 5;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
                //printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          ConfFlag[e]=mesh.EToE[idf+1]; 
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+1]; 
          //printf("confFlag[e]=%d\n",mesh.EToE[idf+1] );
          RefFlag[e]=0;
          nn++;
          nv++;
        }

        else if (FaceFlag[idf+2]==1)
        {
          // To have unique global vertex number
          //hlong  Local_id = 2*(level+3)+(mesh.Nfaces+3*(level+3))*e+mesh.Nnodes+EToRefLevel[e];          
          //hlong Neigh_id = mesh.EToF[idf+2]*(level+3)+(mesh.Nfaces+3*(level+3))*mesh.EToE[idf+2]+mesh.Nnodes+EToRefLevel[mesh.EToE[idf+2]];
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          // Second version
          hlong vA = mesh.EToV[id + 0]; 
          hlong vB = mesh.EToV[id + 2];

          hlong vmin = (vA < vB) ? vA : vB;
          hlong vmax = (vA > vB) ? vA : vB;

          hlong newNode = mesh.Nnodes + (vmin * 31 + vmax);
          //hlong newNode = mesh.Nnodes + ( (vA + vB) * (vA + vB + 1) / 2 + vB );
           printf("newNode=%d,vA=%d,vB=%d\n",newNode,vA,vB );
          printf("e2=%d,id_new2=%lld,mesh.EToF=%d\n",e,mesh.Nelements+nn,mesh.EToF[idf+1] );
          //hlong newNode = v0+v2+mesh.Nnodes;
           //printf("Neigh_e=%lld,MeshLevel=%d,NeighLevel=%d\n",mesh.EToE[idf+2],EToRefLevel[e],EToRefLevel[mesh.EToE[idf+2]]);
          //printf("Local_id=%lld,Neigh_id=%lld,newNode=%lld\n",Local_id,Neigh_id,newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5*(mesh.EX[id+2]+mesh.EX[id+0]);

          EX_new[id_new+0] = mesh.EX[id+1];
          EX_new[id_new+1] = mesh.EX[id+2];
          EX_new[id_new+2] = EX_new[id+2];

          EY_new[id+2] = 0.5*(mesh.EY[id+2]+mesh.EY[id+0]);

          EY_new[id_new+0] = mesh.EY[id+1];
          EY_new[id_new+1] = mesh.EY[id+2];
          EY_new[id_new+2] = EY_new[id+2];

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = -1;
          EToB_new[id+2] = mesh.EToB[id+2];

          EToB_new[id_new+0] = mesh.EToB[id+1];
          EToB_new[id_new+1] = mesh.EToB[id+2];
          EToB_new[id_new+2] = -1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          
          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 2;
              const dlong id2_int = 3;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
                //printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

            /* for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }*/
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          ConfFlag[e]=mesh.EToE[idf+2]; 
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+2]; 
          //printf("confFlag[e]=%d,LocalFace=%d\n",mesh.EToE[idf+2],mesh.EToF[idf+2] );
          RefFlag[e]=0;
          nn++;
          nv++;
          } 
      }
        
  }   
      *new_vertex += nv;
      *NN += nn;      
}

void adaptivity_t::BisectNew(  memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          memory<hlong>& new_v_id,
                          dlong Nrefine,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel,
                          dlong ConfLevel){
  printf("bisect New Starts!\n");
  hlong nv = 0 ; // new vertex
  hlong nn = 0 ; // Counts each refinement
  dlong const level = RefLevel;
  for (dlong i = 0; i < Nrefine; ++i)
  {
      
      //int e = Ref[i];
      const int e = new_v_id[i*2+0];


      if ( e!=-1 && EToRefLevel[e]<ConfLevel)
      {
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 
        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
        if (FaceFlag[idf+0]==1)
        { 

          // To have unique global vertex number
  
          hlong newNode = new_v_id[i*2+1];
          printf("e=%lld,new_e=%lld,newNode=%lld,Level = %d\n",e,mesh.Nelements+nn,newNode,EToRefLevel[e]+1);
          // 1st child vertex ids
          EToV_new[id+0] = mesh.EToV[id+2];
          EToV_new[id+1] = mesh.EToV[id+0];
          EToV_new[id+2] = newNode;
          
          // 2nd child vertex ids
          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          // 2nd child vertex coordinates
          EX_new[id_new+0] = mesh.EX[id+1];
          EX_new[id_new+1] = mesh.EX[id+2];
          EX_new[id_new+2] = 0.5f*(mesh.EX[id+0]+mesh.EX[id+1]);

          // 1st child vertex coordinates
          EX_new[id+0] = mesh.EX[id+2];
          EX_new[id+1] = mesh.EX[id+0];
          EX_new[id+2] = 0.5f*(mesh.EX[id+0]+mesh.EX[id+1]);

          // 2nd child vertex coordinates
          EY_new[id_new+0] = mesh.EY[id+1];
          EY_new[id_new+1] = mesh.EY[id+2];
          EY_new[id_new+2] = 0.5f*(mesh.EY[id+0]+mesh.EY[id+1]);

          // 1st child vertex coordinates
          EY_new[id+0] = mesh.EY[id+2];
          EY_new[id+1] = mesh.EY[id+0];
          EY_new[id+2] = 0.5f*(mesh.EY[id+0]+mesh.EY[id+1]);
          

          // Boundary information

          // 1st child
          EToB_new[id+0] = mesh.EToB[id+2]; 
          EToB_new[id+1] = mesh.EToB[id+0]; // New face will be inner in any situation
          EToB_new[id+2] = -1;
          // 2nd child
          EToB_new[id_new+0] = mesh.EToB[id+1];
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+0]; // New face will be inner in any situation

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];
          

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)+0]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;

          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 2;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
              // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 0;
              const dlong id2_int = 1;
           
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          if (ConfFlag[e] == -1)
          {ConfFlag[e]=mesh.EToE[idf+0];}
        
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+0]; 
           printf("confFlag[%d]=%d\n",mesh.Nelements+nn,mesh.EToE[idf+0] );
           printf("confFlag[%d]=%d\n",e,mesh.EToE[idf+0] );
          RefFlag[e]=0;
          nn++;
          nv++;
        }

        else if (FaceFlag[idf+1]==1)
        {
          
          // To have unique global vertex number
          hlong newNode = new_v_id[i*2+1];
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v2;
          EToV_new[id_new+1] = v0;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(mesh.EX[id+1]+mesh.EX[id+2]); 

          
          EX_new[id_new+0] = mesh.EX[id+2];
          EX_new[id_new+1] = mesh.EX[id+0]; 
          EX_new[id_new+2] = EX_new[id+2];  

          EY_new[id+2] = 0.5f*(mesh.EY[id+1]+mesh.EY[id+2]); 

           
          EY_new[id_new+0] = mesh.EY[id+2];
          EY_new[id_new+1] = mesh.EY[id+0]; 
          EY_new[id_new+2] = EY_new[id+2]; 

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[id+1];
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = mesh.EToB[id+2];
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+1];

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
      
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 5;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
              // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 4;
              const dlong id2_int = 5;
          
            //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];

            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          if (ConfFlag[e] == -1)
          {ConfFlag[e]=mesh.EToE[idf+1];}
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+1]; 
          printf("confFlag[%d]=%d\n",mesh.Nelements+nn,mesh.EToE[idf+1] );
          printf("confFlag[%d]=%d\n",e,mesh.EToE[idf+1] );
          RefFlag[e]=0;
          nn++;
          nv++;
        }

        else if (FaceFlag[idf+2]==1)
        {
          // To have unique global vertex number
          hlong newNode = new_v_id[i*2+1];
          
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(mesh.EX[id+2]+mesh.EX[id+0]);

          EX_new[id_new+0] = mesh.EX[id+1];
          EX_new[id_new+1] = mesh.EX[id+2];
          EX_new[id_new+2] = EX_new[id+2];

          EY_new[id+2] = 0.5f*(mesh.EY[id+2]+mesh.EY[id+0]);

          EY_new[id_new+0] = mesh.EY[id+1];
          EY_new[id_new+1] = mesh.EY[id+2];
          EY_new[id_new+2] = EY_new[id+2];

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = -1;
          EToB_new[id+2] = mesh.EToB[id+2];

          EToB_new[id_new+0] = mesh.EToB[id+1];
          EToB_new[id_new+1] = mesh.EToB[id+2];
          EToB_new[id_new+2] = -1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          
          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(level+3)]   = e;
          PToC[e*(level+3)+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*(level+3)]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 2;
              const dlong id2_int = 3;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          if (ConfFlag[e] == -1)
          {ConfFlag[e]=mesh.EToE[idf+2];}
          ConfFlag[mesh.Nelements+nn]=mesh.EToE[idf+2]; 
           printf("confFlag[%d]=%d\n",mesh.Nelements+nn,mesh.EToE[idf+2] );
           printf("confFlag[%d]=%d\n",e,mesh.EToE[idf+2] );
          RefFlag[e]=0;
          nn++;
          nv++;
          } 
      }
        
  }   
      //*new_vertex += nv;
      *NN += nn;      
      printf("BisectNew is done!\n");
}

// Red refinement of triangle i.e. regular refinement
void adaptivity_t::Red( deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          memory<dlong>& RedFlag,
                          memory<hlong>& new_v_id,
                          dlong Nrefine,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nn = 0 ; // Counts each refinement
  dlong const level = RefLevel;
  dlong const stride = 16*level;
  //dlong const MAX_REFINEMENT_LEVEL = 1;

      //int e = Ref[i];

  for (int i = 0; i < Nrefine; ++i)
  {
      const dlong e = new_v_id[i*4+0];
      const dlong id = e*mesh.Nverts; 
      if (e!=-1 && EToRefLevel[e]<RefLevel)
      {
        printf("Red Refined Element = %d\n",e );

        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces

          // Uniquely number each vertex 

          hlong vmin = (v0 < v1) ? v0 : v1;
          hlong vmax = (v0 > v1) ? v0 : v1;

          //hlong v3 = mesh.Nnodes + (vmin * 10000 + vmax);
          
          hlong v3 = new_v_id[i*4+1];

          vmin = (v1 < v2) ? v1 : v2;
          vmax = (v1 > v2) ? v1 : v2;

          //hlong v4 = mesh.Nnodes + (vmin * 10000 + vmax);

          hlong v4 = new_v_id[i*4+2];

          vmin = (v0 < v2) ? v0 : v2;
          vmax = (v0 > v2) ? v0 : v2;

          //hlong v5 = mesh.Nnodes + (vmin * 10000 + vmax);

          hlong v5 = new_v_id[i*4+3];
          printf("Red refined element e=%d with new ids v3=%d,v4=%d,v5=%d",e,v3,v4,v5);
          //hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_5 = 2+mesh.Nfaces*e+mesh.Nnodes;

          //hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          //hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;
          //hlong Neigh_id_5 = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;

          //hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          //hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;
          //hlong v5 = (Local_id_5>=Neigh_id_5)? Local_id_5:Neigh_id_5 ;

          // Replace original with triangle connecting edge centers
          EToV_new[id+0] = v3;
          EToV_new[id+1] = v4;
          EToV_new[id+2] = v5;

          // Add extra triangles to EToV
          EToV_new[id_new+0] = v0; EToV_new[id_new+1] = v3; EToV_new[id_new+2] = v5;
          EToV_new[id_new+3] = v1; EToV_new[id_new+4] = v4; EToV_new[id_new+5] = v3;
          EToV_new[id_new+6] = v2; EToV_new[id_new+7] = v5; EToV_new[id_new+8] = v4;
          
          // Find vertex locations of elements to be refined
          dfloat  Xv3 = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);  
          dfloat  Xv4 = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);
          dfloat  Xv5 = 0.5*(mesh.EX[id+2]+mesh.EX[id+0]);
          
          dfloat  Yv3 = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);  
          dfloat  Yv4 = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          dfloat  Yv5 = 0.5*(mesh.EY[id+2]+mesh.EY[id+0]);
          
          // Add coordinates for refined edge centers
          EX_new[id+0] = Xv3; EY_new[id+0] = Yv3;
          EX_new[id+1] = Xv4; EY_new[id+1] = Yv4;
          EX_new[id+2] = Xv5; EY_new[id+2] = Yv5;
          
          EX_new[id_new+0] = mesh.EX[id+0]; EX_new[id_new+1] = Xv3; EX_new[id_new+2] = Xv5;
          EX_new[id_new+3] = mesh.EX[id+1]; EX_new[id_new+4] = Xv4; EX_new[id_new+5] = Xv3;
          EX_new[id_new+6] = mesh.EX[id+2]; EX_new[id_new+7] = Xv5; EX_new[id_new+8] = Xv4;

          EY_new[id_new+0] = mesh.EY[id+0]; EY_new[id_new+1] = Yv3; EY_new[id_new+2] = Yv5;
          EY_new[id_new+3] = mesh.EY[id+1]; EY_new[id_new+4] = Yv4; EY_new[id_new+5] = Yv3;
          EY_new[id_new+6] = mesh.EY[id+2]; EY_new[id_new+7] = Yv5; EY_new[id_new+8] = Yv4;

          // Create boundary condition type for refined elements

          EToB_new[id+0] = -1; 
          EToB_new[id+1] = -1; 
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = mesh.EToB[id+0]; EToB_new[id_new+3] = mesh.EToB[id+1]; EToB_new[id_new+6] = mesh.EToB[id+2];  
          EToB_new[id_new+1] = -1;              EToB_new[id_new+4] = -1;              EToB_new[id_new+7] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+2]; EToB_new[id_new+5] = mesh.EToB[id+0]; EToB_new[id_new+8] = mesh.EToB[id+1];



          // Update lists related to AMR
          // Update cell levels
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn+0)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+1)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+2)] = EToRefLevel[e];

          // Parent to Child Connection
          PToC[e*stride+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[e*stride+(EToRefLevel[e]-1)*4+1] = mesh.Nelements+nn+0;
          PToC[e*stride+(EToRefLevel[e]-1)*4+2] = mesh.Nelements+nn+1;
          PToC[e*stride+(EToRefLevel[e]-1)*4+3] = mesh.Nelements+nn+2; 

          PToC[(mesh.Nelements+nn+0)*stride+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn+0)*stride+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn+0)*stride+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn+0)*stride+(EToRefLevel[e]-1)*4+3] = -1; 

          PToC[(mesh.Nelements+nn+1)*stride+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn+1)*stride+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn+1)*stride+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn+1)*stride+(EToRefLevel[e]-1)*4+3] = -1; 

          PToC[(mesh.Nelements+nn+2)*stride+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn+2)*stride+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn+2)*stride+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn+2)*stride+(EToRefLevel[e]-1)*4+3] = -1; 

          
          
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn+0)*(level+3)+EToRefLevel[e]-1] = 2;
          IntFlag[(mesh.Nelements+nn+1)*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn+2)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=2;
          SplitFlag[mesh.Nelements+nn+0]=2;
          SplitFlag[mesh.Nelements+nn+1]=2;
          SplitFlag[mesh.Nelements+nn+2]=2;


           const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(4*Np*Np,0);
          o_IMRed.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn+0)*Np;
              const dlong id3 = (mesh.Nelements+nn+1)*Np;
              const dlong id4 = (mesh.Nelements+nn+2)*Np;
              
              const dlong id1_int = 0;
              const dlong id2_int = 1;
              const dlong id3_int = 2;
              const dlong id4_int = 3;
                 //printf("id1_int=%d,id2_int=%d\n", id1_int,id2_int);
              //printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; dfloat qn3=0.; dfloat qn4=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];
                const dfloat Ii3 = IM[n+i*Np+id3_int*Np*Np];
                     qn3 += Ii3*Qold[id1+i];
                const dfloat Ii4 = IM[n+i*Np+id4_int*Np*Np];
                     qn4 += Ii4*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
              Q[id3+n] = qn3;
              Q[id4+n] = qn4;
              //  printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

                        for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
              Qold[id3+n] = Q[id3+n];
              Qold[id4+n] = Q[id4+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

          ConfFlag[(mesh.Nelements+nn+0)] = mesh.EToE[id+2];
          ConfFlag[(mesh.Nelements+nn+1)] = mesh.EToE[id+0];
          ConfFlag[(mesh.Nelements+nn+2)] = mesh.EToE[id+1];   

          RedFlag[e*(level+3)+(EToRefLevel[e])-1] = 1;
          RedFlag[(mesh.Nelements+nn+0)*(level+3)+(EToRefLevel[e])-1] = 1;
          RedFlag[(mesh.Nelements+nn+1)*(level+3)+(EToRefLevel[e])-1] = 1;
          RedFlag[(mesh.Nelements+nn+2)*(level+3)+(EToRefLevel[e])-1] = 1;     

          

          nn=nn+3;
          //new_vertex++;
      }

      *NN = nn;             
  }
}   

// Red refinement of triangle i.e. regular refinement
void adaptivity_t::RedOld( deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          memory<dlong>& RedFlag,
                          memory<hlong>& new_v_id,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nn = 0 ; // Counts each refinement
  dlong const level = RefLevel;
  //dlong const MAX_REFINEMENT_LEVEL = 1;

  for (int e = 0; e < mesh.Nelements; ++e)
  {
      const dlong id = e*mesh.Nverts; 
      if (RefFlag[e]==1 && (EToRefLevel[e]<RefLevel))
      {
        printf("Red Refined Element = %d\n",e );

        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces

          // Uniquely number each vertex 

          hlong vmin = (v0 < v1) ? v0 : v1;
          hlong vmax = (v0 > v1) ? v0 : v1;

          hlong v3 = mesh.Nnodes + (vmin * 10000 + vmax);
          
          //hlong v3 = new_v_id[e*4+1];

          vmin = (v1 < v2) ? v1 : v2;
          vmax = (v1 > v2) ? v1 : v2;

          hlong v4 = mesh.Nnodes + (vmin * 10000 + vmax);

          //hlong v4 = new_v_id[e*4+2];

          vmin = (v0 < v2) ? v0 : v2;
          vmax = (v0 > v2) ? v0 : v2;

          hlong v5 = mesh.Nnodes + (vmin * 10000 + vmax);

          //hlong v5 = new_v_id[e*4+3];

          //hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_5 = 2+mesh.Nfaces*e+mesh.Nnodes;

          //hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          //hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;
          //hlong Neigh_id_5 = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;

          //hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          //hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;
          //hlong v5 = (Local_id_5>=Neigh_id_5)? Local_id_5:Neigh_id_5 ;

          // Replace original with triangle connecting edge centers
          EToV_new[id+0] = v3;
          EToV_new[id+1] = v4;
          EToV_new[id+2] = v5;

          // Add extra triangles to EToV
          EToV_new[id_new+0] = v0; EToV_new[id_new+1] = v3; EToV_new[id_new+2] = v5;
          EToV_new[id_new+3] = v1; EToV_new[id_new+4] = v4; EToV_new[id_new+5] = v3;
          EToV_new[id_new+6] = v2; EToV_new[id_new+7] = v5; EToV_new[id_new+8] = v4;
          
          // Find vertex locations of elements to be refined
          dfloat  Xv3 = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);  
          dfloat  Xv4 = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);
          dfloat  Xv5 = 0.5*(mesh.EX[id+2]+mesh.EX[id+0]);
          
          dfloat  Yv3 = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);  
          dfloat  Yv4 = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          dfloat  Yv5 = 0.5*(mesh.EY[id+2]+mesh.EY[id+0]);
          
          // Add coordinates for refined edge centers
          EX_new[id+0] = Xv3; EY_new[id+0] = Yv3;
          EX_new[id+1] = Xv4; EY_new[id+1] = Yv4;
          EX_new[id+2] = Xv5; EY_new[id+2] = Yv5;
          
          EX_new[id_new+0] = mesh.EX[id+0]; EX_new[id_new+1] = Xv3; EX_new[id_new+2] = Xv5;
          EX_new[id_new+3] = mesh.EX[id+1]; EX_new[id_new+4] = Xv4; EX_new[id_new+5] = Xv3;
          EX_new[id_new+6] = mesh.EX[id+2]; EX_new[id_new+7] = Xv5; EX_new[id_new+8] = Xv4;

          EY_new[id_new+0] = mesh.EY[id+0]; EY_new[id_new+1] = Yv3; EY_new[id_new+2] = Yv5;
          EY_new[id_new+3] = mesh.EY[id+1]; EY_new[id_new+4] = Yv4; EY_new[id_new+5] = Yv3;
          EY_new[id_new+6] = mesh.EY[id+2]; EY_new[id_new+7] = Yv5; EY_new[id_new+8] = Yv4;

          // Create boundary condition type for refined elements

          EToB_new[id+0] = -1; 
          EToB_new[id+1] = -1; 
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = mesh.EToB[id+0]; EToB_new[id_new+3] = mesh.EToB[id+1]; EToB_new[id_new+6] = mesh.EToB[id+2];  
          EToB_new[id_new+1] = -1;              EToB_new[id_new+4] = -1;              EToB_new[id_new+7] = -1;
          EToB_new[id_new+2] = mesh.EToB[id+2]; EToB_new[id_new+5] = mesh.EToB[id+0]; EToB_new[id_new+8] = mesh.EToB[id+1];



          // Update lists related to AMR
          // Update cell levels
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn+0)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+1)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+2)] = EToRefLevel[e];

          // Parent to Child Connection
          PToC[e*(level+3)+EToRefLevel[e]*4+0]   = e; 
          PToC[e*(level+3)+EToRefLevel[e]*4+1] = mesh.Nelements+nn+0;
          PToC[e*(level+3)+EToRefLevel[e]*4+2] = mesh.Nelements+nn+1;
          PToC[e*(level+3)+EToRefLevel[e]*4+3] = mesh.Nelements+nn+2; 
          
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn+0)*(level+3)+EToRefLevel[e]-1] = 2;
          IntFlag[(mesh.Nelements+nn+1)*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn+2)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=2;
          SplitFlag[mesh.Nelements+nn+0]=2;
          SplitFlag[mesh.Nelements+nn+1]=2;
          SplitFlag[mesh.Nelements+nn+2]=2;


           const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(4*Np*Np,0);
          o_IMRed.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn+0)*Np;
              const dlong id3 = (mesh.Nelements+nn+1)*Np;
              const dlong id4 = (mesh.Nelements+nn+2)*Np;
              
              const dlong id1_int = 0;
              const dlong id2_int = 1;
              const dlong id3_int = 2;
              const dlong id4_int = 3;
                 //printf("id1_int=%d,id2_int=%d\n", id1_int,id2_int);
              //printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; dfloat qn3=0.; dfloat qn4=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];
                const dfloat Ii3 = IM[n+i*Np+id3_int*Np*Np];
                     qn3 += Ii3*Qold[id1+i];
                const dfloat Ii4 = IM[n+i*Np+id4_int*Np*Np];
                     qn4 += Ii4*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
              Q[id3+n] = qn3;
              Q[id4+n] = qn4;
              //  printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

                        for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
              Qold[id3+n] = Q[id3+n];
              Qold[id4+n] = Q[id4+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

          ConfFlag[(mesh.Nelements+nn+0)] = mesh.EToE[id+2];
          ConfFlag[(mesh.Nelements+nn+1)] = mesh.EToE[id+0];
          ConfFlag[(mesh.Nelements+nn+2)] = mesh.EToE[id+1];   

          RedFlag[e] = 1;
          RedFlag[(mesh.Nelements+nn+0)] = 1;
          RedFlag[(mesh.Nelements+nn+1)] = 1;
          RedFlag[(mesh.Nelements+nn+2)] = 1;     

          

          nn=nn+3;
          //new_vertex++;
      }

      *NN = nn;             
  }
}   

// Red refinement of triangle i.e. regular refinement
void adaptivity_t::RedLocal( dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          memory<dlong>& RedFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nn = 0 ; // Counts each refinement
  dlong const MAX_REFINEMENT_LEVEL = 1;
      const dlong id = e*mesh.Nverts; 
      if ((EToRefLevel[e]<MAX_REFINEMENT_LEVEL))
      {

        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces

          // Uniquely number each vertex 

          hlong vmin = (v0 < v1) ? v0 : v1;
          hlong vmax = (v0 > v1) ? v0 : v1;

          hlong v3 = mesh.Nnodes + (vmin * 10000 + vmax);
          
          vmin = (v1 < v2) ? v1 : v2;
          vmax = (v1 > v2) ? v1 : v2;

          hlong v4 = mesh.Nnodes + (vmin * 10000 + vmax);

          vmin = (v0 < v2) ? v0 : v2;
          vmax = (v0 > v2) ? v0 : v2;

          hlong v5 = mesh.Nnodes + (vmin * 10000 + vmax);

          //hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;
          //hlong Local_id_5 = 2+mesh.Nfaces*e+mesh.Nnodes;

          //hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          //hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;
          //hlong Neigh_id_5 = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;

          //hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          //hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;
          //hlong v5 = (Local_id_5>=Neigh_id_5)? Local_id_5:Neigh_id_5 ;

          // Replace original with triangle connecting edge centers
          EToV_new[id+0] = v3;
          EToV_new[id+1] = v4;
          EToV_new[id+2] = v5;

          // Add extra triangles to EToV
          EToV_new[id_new+0] = v0; EToV_new[id_new+1] = v3; EToV_new[id_new+2] = v5;
          EToV_new[id_new+3] = v1; EToV_new[id_new+4] = v4; EToV_new[id_new+5] = v3;
          EToV_new[id_new+6] = v2; EToV_new[id_new+7] = v5; EToV_new[id_new+8] = v4;
          
          // Find vertex locations of elements to be refined
          dfloat  Xv3 = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);  
          dfloat  Xv4 = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);
          dfloat  Xv5 = 0.5*(mesh.EX[id+2]+mesh.EX[id+0]);
          
          dfloat  Yv3 = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);  
          dfloat  Yv4 = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          dfloat  Yv5 = 0.5*(mesh.EY[id+2]+mesh.EY[id+0]);
          
          // Add coordinates for refined edge centers
          EX_new[id+0] = Xv3; EY_new[id+0] = Yv3;
          EX_new[id+1] = Xv4; EY_new[id+1] = Yv4;
          EX_new[id+2] = Xv5; EY_new[id+2] = Yv5;
          
          EX_new[id_new+0] = mesh.EX[id+0]; EX_new[id_new+1] = Xv3; EX_new[id_new+2] = Xv5;
          EX_new[id_new+3] = mesh.EX[id+1]; EX_new[id_new+4] = Xv4; EX_new[id_new+5] = Xv3;
          EX_new[id_new+6] = mesh.EX[id+2]; EX_new[id_new+7] = Xv5; EX_new[id_new+8] = Xv4;

          EY_new[id_new+0] = mesh.EY[id+0]; EY_new[id_new+1] = Yv3; EY_new[id_new+2] = Yv5;
          EY_new[id_new+3] = mesh.EY[id+1]; EY_new[id_new+4] = Yv4; EY_new[id_new+5] = Yv3;
          EY_new[id_new+6] = mesh.EY[id+2]; EY_new[id_new+7] = Yv5; EY_new[id_new+8] = Yv4;

          // Create boundary condition type for refined elements

          EToB_new[id+0] = -1; 
          EToB_new[id+1] = -1; 
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = mesh.EToB[id+0]; EToB_new[id_new+3] = mesh.EToB[id+0]; EToB_new[id_new+6] = mesh.EToB[id+1];  
          EToB_new[id_new+1] = -1;              EToB_new[id_new+4] = mesh.EToB[id+1]; EToB_new[id_new+7] = mesh.EToB[id+2];
          EToB_new[id_new+2] = mesh.EToB[id+2]; EToB_new[id_new+5] = -1;              EToB_new[id_new+8] = -1;

          // Update lists related to AMR
          // Update cell levels
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn+0)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+1)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+2)] = EToRefLevel[e];

          // Parent to Child Connection
          if (EToRefLevel[e]==1)
          { PToC[e*4+0]   = e; }
          PToC[e*4+1] = mesh.Nelements+nn+0;
          PToC[e*4+2] = mesh.Nelements+nn+1;
          PToC[e*4+3] = mesh.Nelements+nn+2; 
          
          IntFlag[e*(1+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn+0)*(1+3)+EToRefLevel[e]-1] = 2;
          IntFlag[(mesh.Nelements+nn+1)*(1+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn+2)*(1+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=2;
          SplitFlag[mesh.Nelements+nn+0]=2;
          SplitFlag[mesh.Nelements+nn+1]=2;
          SplitFlag[mesh.Nelements+nn+2]=2;


           const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(4*Np*Np,0);
          o_IMRed.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn+0)*Np;
              const dlong id3 = (mesh.Nelements+nn+1)*Np;
              const dlong id4 = (mesh.Nelements+nn+2)*Np;
              
              const dlong id1_int = 0;
              const dlong id2_int = 1;
              const dlong id3_int = 2;
              const dlong id4_int = 3;
                 //printf("id1_int=%d,id2_int=%d\n", id1_int,id2_int);
              //printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; dfloat qn3=0.; dfloat qn4=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];
                const dfloat Ii3 = IM[n+i*Np+id3_int*Np*Np];
                     qn3 += Ii3*Qold[id1+i];
                const dfloat Ii4 = IM[n+i*Np+id4_int*Np*Np];
                     qn4 += Ii4*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
              Q[id3+n] = qn3;
              Q[id4+n] = qn4;
              //  printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

                        for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
              Qold[id3+n] = Q[id3+n];
              Qold[id4+n] = Q[id4+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

          ConfFlag[(mesh.Nelements+nn+0)] = mesh.EToE[id+2];
          ConfFlag[(mesh.Nelements+nn+1)] = mesh.EToE[id+0];
          ConfFlag[(mesh.Nelements+nn+2)] = mesh.EToE[id+1];   

          RedFlag[e] = 1;
          RedFlag[(mesh.Nelements+nn+0)] = 1;
          RedFlag[(mesh.Nelements+nn+1)] = 1;
          RedFlag[(mesh.Nelements+nn+2)] = 1;     

          

          nn=nn+3;
          //new_vertex++;
      }

      *NN = nn;             
  
}   
// Blue refinement of triangle i.e. double bisect
void adaptivity_t::Blue(memory<dlong>& RefFlag,
                        memory<dlong>& FaceFlag,
                        memory<dfloat>& EX_new,
                        memory<dfloat>& EY_new,
                        memory<hlong>& EToV_new,
                        memory<int>& EToB_new,
                        memory<dlong>& SplitFlag,
                        hlong* NN){

  hlong nn = 0 ; // Counts each refinement
  dlong const MAX_REFINEMENT_LEVEL = 1;

  for (int e = 0; e < mesh.Nelements; ++e)
  {
      const dlong id = e*mesh.Nverts; 
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {
        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
        if (FaceFlag[id+0]==1 /*&& FaceFlag[id+1]==1 && FaceFlag[id+2]==0*/)
        { 

          const dfloat x0 = mesh.EX[id+0]; const dfloat x1 = mesh.EX[id+1]; const dfloat x2 = mesh.EX[id+2];
          const dfloat y0 = mesh.EY[id+0]; const dfloat y1 = mesh.EY[id+1]; const dfloat y2 = mesh.EY[id+2];

          // Find Longest Edge
          const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
          const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
          
          if (mag0 > mag1)
          {
            // Uniquely number each vertex 
          hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;

          hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;

          hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;

          // Replace original with triangle connecting first edge center
          EToV_new[id+1] = v3;

          // Add extra triangles to EToV
          EToV_new[id_new+0] = v3; EToV_new[id_new+1] = v1; EToV_new[id_new+2] = v4;
          EToV_new[id_new+3] = v3; EToV_new[id_new+4] = v4; EToV_new[id_new+5] = v2;
          
          // Find vertex locations of elements to be refined
          dfloat  Xv3 = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);  
          dfloat  Xv4 = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);
          
          dfloat  Yv3 = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);  
          dfloat  Yv4 = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          
          
          // Add coordinates for refined edge centers
          EX_new[id+1] = Xv3; EY_new[id+1] = Yv3;
          
          EX_new[id_new+0] = Xv3; EX_new[id_new+1] = Xv4; EX_new[id_new+2] = mesh.EX[id+2];
          EX_new[id_new+3] = Xv3; EX_new[id_new+4] = mesh.EX[id+1]; EX_new[id_new+5] = Xv4;
        
          EY_new[id_new+0] = Yv3; EY_new[id_new+1] = Yv4; EY_new[id_new+2] = mesh.EY[id+2];
          EY_new[id_new+3] = Yv3; EY_new[id_new+4] = mesh.EY[id+1]; EY_new[id_new+5] = Yv4;

          // Create boundary condition type for refined elements
 
          EToB_new[id+1] = -1; 
          
          EToB_new[id_new+0] = -1;                EToB_new[id_new+3] = mesh.EToB[id+0];   
          EToB_new[id_new+1] = mesh.EToB[id+1];   EToB_new[id_new+4] = mesh.EToB[id+1]; 
          EToB_new[id_new+2] = -1;                EToB_new[id_new+5] = -1;              
          }else
          {
                        // Uniquely number each vertex 
          hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;

          hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;

          hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;

          // Replace original with triangle connecting first edge center
          EToV_new[id+1] = v3;

          // Add extra triangles to EToV
          EToV_new[id_new+0] = v3; EToV_new[id_new+1] = v1; EToV_new[id_new+2] = v4;
          EToV_new[id_new+3] = v3; EToV_new[id_new+4] = v4; EToV_new[id_new+5] = v2;
          
          // Find vertex locations of elements to be refined
          dfloat  Xv4 = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);  
          dfloat  Xv3 = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);
          
          dfloat  Yv4 = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);  
          dfloat  Yv3 = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          
          
          // Add coordinates for refined edge centers
          EX_new[id+1] = Xv3; EY_new[id+1] = Yv3;
          
          EX_new[id_new+0] = Xv3; EX_new[id_new+1] = Xv4; EX_new[id_new+2] = mesh.EX[id+2];
          EX_new[id_new+3] = Xv3; EX_new[id_new+4] = mesh.EX[id+1]; EX_new[id_new+5] = Xv4;
        
          EY_new[id_new+0] = Yv3; EY_new[id_new+1] = Yv4; EY_new[id_new+2] = mesh.EY[id+2];
          EY_new[id_new+3] = Yv3; EY_new[id_new+4] = mesh.EY[id+1]; EY_new[id_new+5] = Yv4;

          // Create boundary condition type for refined elements
 
          EToB_new[id+1] = -1; 
          
          EToB_new[id_new+0] = -1;                EToB_new[id_new+3] = mesh.EToB[id+0];   
          EToB_new[id_new+1] = mesh.EToB[id+1];   EToB_new[id_new+4] = mesh.EToB[id+1]; 
          EToB_new[id_new+2] = -1;                EToB_new[id_new+5] = -1;        
          }

          // Update lists related to AMR
          // Update cell levels
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn+0)] = EToRefLevel[e];
          EToRefLevel[(mesh.Nelements+nn+1)] = EToRefLevel[e];

          // Parent to Child Connection
          if (EToRefLevel[e]==1)
          { PToC[e*4+0]   = e; }
          PToC[e*4+1] = mesh.Nelements+nn+0;
          PToC[e*4+2] = mesh.Nelements+nn+1;
          
          IntFlag[e] = 1;
          IntFlag[(mesh.Nelements+nn+0)] = 2;
          IntFlag[(mesh.Nelements+nn+1)] = 3;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn+0]=1;
          SplitFlag[mesh.Nelements+nn+1]=1;

          nn=nn+2;
          //new_vertex++;
        }
      }

      *NN = nn;             
  }    
        
}
     
    

void adaptivity_t::Green0(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
  dlong const stride = 16*level;
  
      
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 

      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces

          // Uniquely number new vertex 
          //hlong Local_id = 0+idf+mesh.Nnodes;

          // To have unique global vertex number
          //hlong Neigh_id = mesh.EToF[idf+0]+mesh.Nfaces*mesh.EToE[idf+0]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+0];
          printf("newNode=%lld\n",newNode );
          // 1st child vertex ids
          EToV_new[id+0] = v2;
          EToV_new[id+1] = v0;
          EToV_new[id+2] = newNode;
          
          // 2nd child vertex ids
          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          // 2nd child vertex coordinates
          EX_new[id_new+0] = x1;
          EX_new[id_new+1] = x2;
          EX_new[id_new+2] = 0.5f*(x0+x1);

          // 1st child vertex coordinates
          EX_new[id+0] = x2;
          EX_new[id+1] = x0;
          EX_new[id+2] = 0.5f*(x0+x1);

          // 2nd child vertex coordinates
          EY_new[id_new+0] = y1;
          EY_new[id_new+1] = y2;
          EY_new[id_new+2] = 0.5f*(y0+y1);

          // 1st child vertex coordinates
          EY_new[id+0] = y2;
          EY_new[id+1] = y0;
          EY_new[id+2] = 0.5f*(y0+y1);
          

          // Boundary information

          // 1st child
          EToB_new[id+0] = b2; 
          EToB_new[id+1] = b0; // New face will be inner in any situation
          EToB_new[id+2] = -1;
          // 2nd child
          EToB_new[id_new+0] = b1;
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = b0; // New face will be inner in any situation

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];
          
          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+1] = mesh.Nelements+nn;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 

          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 

          
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 1;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 2;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
          const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              //const dlong id1_int = IntFlag[e*(level+3)+EToRefLevel[e]-1]-1;
              const dlong id1_int = 0;
              //const dlong id2_int = IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1]-1;
              const dlong id2_int = 1;
              printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
            //    printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

            for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
            //printf("Q_inbisect=%f\n",Q[1621] );
          //
          //ConfFlag[e]=mesh.EToE[idf+0];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+0];   
          RefFlag[e]=0;
          nn++;
          nv++;
       
      *new_vertex = nv;
      *NN = nn;      
}

void adaptivity_t::Green1(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
  dlong const stride = 16*level;
  
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 


      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
          //hlong Local_id = 1+idf+mesh.Nnodes;

          //hlong Neigh_id = mesh.EToF[idf+1]+mesh.Nfaces*mesh.EToE[idf+1]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+1];
          printf("newNode=%lld\n",newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v2;
          EToV_new[id_new+1] = v0;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(x1+x2); 

          
          EX_new[id_new+0] = x2;
          EX_new[id_new+1] = x0; 
          EX_new[id_new+2] = 0.5f*(x1+x2);  

          EY_new[id+2] = 0.5f*(y1+y2); 

           
          EY_new[id_new+0] = y2;
          EY_new[id_new+1] = y0; 
          EY_new[id_new+2] = 0.5f*(y1+y2); 

          EToB_new[id+0] = b0;
          EToB_new[id+1] = b1;
          EToB_new[id+2] = -1;

          EToB_new[id_new+0] = b2;
          EToB_new[id_new+1] = -1;
          EToB_new[id_new+2] = b1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+1] = mesh.Nelements+nn;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 

          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 
    
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 5;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 4;
              const dlong id2_int = 5;
                 //printf("id1_int=%d,id2_int=%d\n", id1_int,id2_int);
              //printf("id1=%d,id2=%d\n", e*Np,(mesh.Nelements+nn)*Np);
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
              //  printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

                        for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          RefFlag[e]=0;
          //ConfFlag[e]=mesh.EToE[idf+1];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+1];   
          nn++;
          nv++;            
     *new_vertex = nv;
     *NN = nn;      
}

void adaptivity_t::Green2(dlong e,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& new_v_id,
                          memory<dfloat>& EX_new,
                          memory<dfloat>& EY_new,
                          memory<hlong>& EToV_new,
                          memory<int>& EToB_new,
                          memory<dlong>& SplitFlag,
                          hlong* NN,
                          hlong* new_vertex,
                          dlong RefLevel){

  hlong nv = *new_vertex ; // new vertex
  hlong nn = *NN ; // Counts each refinement
  dlong const level = RefLevel;
  dlong const stride = 16*level;
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
      const dlong idf = e*mesh.Nfaces; 
      
      // Parent vertex ids
      const hlong v0 = EToV_new[id + 0];
      const hlong v1 = EToV_new[id + 1];
      const hlong v2 = EToV_new[id + 2];
      
      // Parent coords
      const dfloat x0 = EX_new[id + 0];
      const dfloat x1 = EX_new[id + 1];
      const dfloat x2 = EX_new[id + 2];
      
      const dfloat y0 = EY_new[id + 0];
      const dfloat y1 = EY_new[id + 1];
      const dfloat y2 = EY_new[id + 2];
      
      // Parent BC
      const int b0 = EToB_new[id + 0];
      const int b1 = EToB_new[id + 1];
      const int b2 = EToB_new[id + 2];
        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        printf("Local bisect e=%d\n",e);
        
          //hlong Local_id = 2+idf+mesh.Nnodes;
          
          //hlong Neigh_id = mesh.EToF[idf+2]+mesh.Nfaces*mesh.EToE[idf+2]+mesh.Nnodes;
          //hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          //hlong newNode = new_v_id[e];
          hlong newNode = new_v_id[id+2];
          printf("newNode=%lld\n",newNode );
          EToV_new[id+2] = newNode;

          EToV_new[id_new+0] = v1;
          EToV_new[id_new+1] = v2;
          EToV_new[id_new+2] = newNode;

          EX_new[id+2] = 0.5f*(x2+x0);

          EX_new[id_new+0] = x1;
          EX_new[id_new+1] = x2;
          EX_new[id_new+2] = 0.5f*(x2+x0);

          EY_new[id+2] = 0.5f*(y2+y0);

          EY_new[id_new+0] = y1;
          EY_new[id_new+1] = y2;
          EY_new[id_new+2] = 0.5f*(y2+y0);

          EToB_new[id+0] = b0;
          EToB_new[id+1] = -1;
          EToB_new[id+2] = b2;

          EToB_new[id_new+0] = b1;
          EToB_new[id_new+1] = b2;
          EToB_new[id_new+2] = -1;

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e];

          PCS[e*3+0] = e;
          PCS[e*3+1] = e;
          PCS[e*3+2] = mesh.Nelements+nn;

          PCS[(mesh.Nelements+nn)*3+0] = e;
          PCS[(mesh.Nelements+nn)*3+1] = (mesh.Nelements+nn);
          PCS[(mesh.Nelements+nn)*3+2] = e;

          PToC[e*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+1] = mesh.Nelements+nn;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[e*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 

          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+0]   = e; 
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+1] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+2] = -1;
          PToC[(mesh.Nelements+nn)*(stride)+(EToRefLevel[e]-1)*4+3] = -1; 
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e*(level+3)+EToRefLevel[e]-1] = 3;
          IntFlag[(mesh.Nelements+nn)*(level+3)+EToRefLevel[e]-1] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;
            const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
          memory<dfloat> IM(6*Np*Np,0);
          o_IM.copyTo(IM);
          
                    // Update solution Local
              const dlong id1 = e*Np;
              const dlong id2 = (mesh.Nelements+nn)*Np;
              
              const dlong id1_int = 2;
              const dlong id2_int = 3;
          
              //const dlong n=0;        
            for(int n=0;n<Np;++n){
              dfloat qn1=0.; dfloat qn2=0.; 
              for(int i=0;i<Np;++i){
                const dfloat Ii1 = IM[n+i*Np+id1_int*Np*Np];
                     qn1 += Ii1*Qold[id1+i];
                const dfloat Ii2 = IM[n+i*Np+id2_int*Np*Np];
                     qn2 += Ii2*Qold[id1+i];                         
                
              }                  
              Q[id1+n] = qn1;
              Q[id2+n] = qn2;
                //printf("Q[%d]=%f,qn1[%d]=%f\n",id2+n,Q[id2+n],id2+n,qn2);
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

             for(int n=0;n<Np;++n){
                        
              Qold[id1+n] = Q[id1+n];
              Qold[id2+n] = Q[id2+n];
                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
            SplitFlag[e] = 0;
            SplitFlag[mesh.Nelements+nn]=0;
          //
          RefFlag[e]=0;
          //ConfFlag[e]=mesh.EToE[idf+2];   
          //ConfFlag[(mesh.Nelements+nn)]=mesh.EToE[idf+2];   
          nn++;
          nv++;
      *new_vertex = nv;
      *NN = nn;      
}




}