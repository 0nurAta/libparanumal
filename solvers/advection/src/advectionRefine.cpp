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

#include "advection.hpp"

// 
void advection_t::Refine(memory<dfloat>& Q,
                         memory<dlong>& RefFlag,
                         dlong Nrefine){

  
  dlong MAX_REFINEMENT_LEVEL = 1;
  // extract q halo on DEVICE
  //traceHalo.ExchangeStart(o_Q, 1);

  // Store old info & Allocate new arrays
  
  // Solution info
  //memory <dfloat> qold = Q;
  
  // Element to vertex & Element to boundary connectivity 
  //memory<hlong>EToV_old(mesh.Nelements*mesh.Nverts);
  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  //memory<dfloat>EX_old(mesh.Nelements*mesh.Nverts);
  //memory<dfloat>EY_old(mesh.Nelements*mesh.Nverts);

  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);
   //printf("EToVsize=%ld\n",EToV_new.size());
  memory<dlong> SplitFlag(2*mesh.Nelements);
  // Copy old Element to Vertex Connectivity to the New One
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    for (int n = 0; n < 3; ++n)
    {
    const dlong id = e*mesh.Nverts+n;
    EToB_new[id] = mesh.EToB[id];
    EToV_new[id] = mesh.EToV[id];
    EX_new[id] = mesh.EX[id]; 
    EY_new[id] = mesh.EY[id]; 
    //printf("EToV=%lld\n",mesh.EToV[id]);
    //printf("EX=%g\n",mesh.EX[id]);
    //printf("EY=%f\n",mesh.EY[id]);
    }
  }
  
  hlong new_vertex = 0; // Counts each new_vertex that will be created
  // Determine Triangles To be Refined
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      
      const dlong id = e*mesh.Nverts; 
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {
        

        // Extract Vertex Number of Element to Refine
        //const hlong v0 = mesh.EToV[id+0]; 
        //const hlong v1 = mesh.EToV[id+1]; 
        //const hlong v2 = mesh.EToV[id+2];
        
        // Find vertex locations of elements to be refined
        const dfloat x0 = mesh.EX[id+0]; const dfloat x1 = mesh.EX[id+1]; const dfloat x2 = mesh.EX[id+2];    
        const dfloat y0 = mesh.EY[id+0]; const dfloat y1 = mesh.EY[id+1]; const dfloat y2 = mesh.EY[id+2];

        // Find Longest Edge
        const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
        const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));

        dlong Le;
        if (mag0 >= mag1 && mag0 >= mag2) {
          Le = id+0;
        } else if (mag1 >= mag2) {
          Le = id+1;
        } else {
          Le = id+2;
        }

        if (Le==id+0)
        {
          hlong const neighbor_id = mesh.EToE[id+0];
          //printf("hello1\n");
          //printf("neighbor_id%lld\n",RefFlag[neighbor_id]==0);
          //printf("mesh.EToE%lld\n",mesh.EToE[id+0]);
          if(RefFlag[neighbor_id]!=1 && mesh.EToE[id+0]!=-1){
             RefFlag[neighbor_id]=1;
             Nrefine++; 
             new_vertex--;
          }
        }

        if (Le==id+1)
        {         
          hlong const neighbor_id = mesh.EToE[id+1];

          if(RefFlag[neighbor_id]!=1 && mesh.EToE[id+1]!=-1){
             RefFlag[neighbor_id]=1;
             Nrefine++; 
             new_vertex--;
          }
        }

        if (Le==id+2)
        {          
          hlong const neighbor_id = mesh.EToE[id+2];

          if(RefFlag[neighbor_id]!=1 && mesh.EToE[id+2]!=-1){
             RefFlag[neighbor_id]=1; 
             Nrefine++;
             new_vertex--;
          }        
        } 

      }             
  }   

  // Refinement Loop
  // Determine ids of new vertices and EToV
  
  hlong nn = 0 ; // Counts each refinement
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);

  for (int e = 0; e < mesh.Nelements; ++e)
  {


      const dlong id = e*mesh.Nverts; 
      if (RefFlag[e]==1 && EToRefLevel[e]<MAX_REFINEMENT_LEVEL)
      {
        // Extract Vertex Number of Element to Refine
        const hlong v0 = mesh.EToV[id+0]; 
        const hlong v1 = mesh.EToV[id+1]; 
        const hlong v2 = mesh.EToV[id+2];
        
        // Find vertex locations of elements to be refined
        const dfloat x0 = mesh.EX[id+0]; const dfloat x1 = mesh.EX[id+1]; const dfloat x2 = mesh.EX[id+2];    
        const dfloat y0 = mesh.EY[id+0]; const dfloat y1 = mesh.EY[id+1]; const dfloat y2 = mesh.EY[id+2];

        // Find Longest Edge
        const dfloat mag0 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)); 
        const dfloat mag1 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        const dfloat mag2 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));

        dlong Le;
        if (mag0 >= mag1 && mag0 >= mag2) {
          Le = id+0;
        } else if (mag1 >= mag2) {
          Le = id+1;
        } else {
          Le = id+2;
        }

        
        // Number the new Vertex at the Longest Edge
        const hlong id_new = (mesh.Nelements+nn)*mesh.Nverts;
        
        //hlong newNode = mesh.Nnodes;
        //mesh.Nnodes++;
        
        //printf("e=%d\n",e);
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
        if (Le==id+0)
        { 
          // Uniquely number new vertex 
          hlong Local_id = 0+mesh.Nfaces*e+mesh.Nnodes;

          hlong Neigh_id = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          
          EToV_new[id+1] = newNode;
          EToV_new[id_new+0] = newNode;
          EToV_new[id_new+1] = v1;
          EToV_new[id_new+2] = v2;

          EX_new[id+1] = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);
          EX_new[id_new+0] = 0.5*(mesh.EX[id+0]+mesh.EX[id+1]);
          EX_new[id_new+1] = mesh.EX[id+1];
          EX_new[id_new+2] = mesh.EX[id+2];

          EY_new[id+1] = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);
          EY_new[id_new+0] = 0.5*(mesh.EY[id+0]+mesh.EY[id+1]);
          EY_new[id_new+1] = mesh.EY[id+1];
          EY_new[id_new+2] = mesh.EY[id+2];

          EToB_new[id+0] = mesh.EToB[id+0]; 
          EToB_new[id+1] = -1; // New face will be inner in any situation
          EToB_new[id+2] = mesh.EToB[id+2]; 
          EToB_new[id_new+0] = mesh.EToB[id+0];
          EToB_new[id_new+1] = mesh.EToB[id+1];
          EToB_new[id_new+2] = -1; // New face will be inner in any situation

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[(mesh.Nelements+nn)]+1;
          PToC[e*2]   = e;
          PToC[e*2+1] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*2]   = e;
          PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e] = 1;
          IntFlag[(mesh.Nelements+nn)] = 1;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=2;

          nn++;
          new_vertex++;
        }

        if (Le==id+1)
        {
          hlong Local_id = 1+mesh.Nfaces*e+mesh.Nnodes;
          hlong Neigh_id = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;
          hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          
          EToV_new[id+2] = newNode;
          EToV_new[id_new+0] = v0;
          EToV_new[id_new+1] = newNode;
          EToV_new[id_new+2] = v2;

          EX_new[id_new+1] = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);  
          EX_new[id+2] = 0.5*(mesh.EX[id+1]+mesh.EX[id+2]);  
          EX_new[id_new+0] = mesh.EX[id+0];
          EX_new[id_new+2] = mesh.EX[id+2]; 
          
          EY_new[id_new+1] = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          EY_new[id+2] = 0.5*(mesh.EY[id+1]+mesh.EY[id+2]);
          EY_new[id_new+0] = mesh.EY[id+0];
          EY_new[id_new+2] = mesh.EY[id+2];

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[id+1];
          EToB_new[id+2] = -1;
          EToB_new[id_new+0] = -1;
          EToB_new[id_new+1] = mesh.EToB[id+1];
          EToB_new[id_new+2] = mesh.EToB[id+2];

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[(mesh.Nelements+nn)]+1;
          PToC[e*2]   = e;
          PToC[e*2+1] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*2]   = e;
          PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e] = 5;
          IntFlag[(mesh.Nelements+nn)] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          nn++;
          new_vertex++;
        }

        if (Le==id+2)
        {
          hlong Local_id = 2+mesh.Nfaces*e+mesh.Nnodes;
          hlong Neigh_id = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;
          hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          
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
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[(mesh.Nelements+nn)]+1;
          PToC[e*2]   = e;
          PToC[e*2+1] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*2]   = e;
          PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e] = 3;
          IntFlag[(mesh.Nelements+nn)] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          nn++;
          new_vertex++;
          } 

        
        //id_new_node++;
        //printf("nn=%d\n",nn);
      }             
  }   

  // Coarsening Loop
  

      if (Nrefine!=0&& nn!=0)
      {

        mesh.EToV.free();
        mesh.EToV = EToV_new;
        EToV_new.free();

        mesh.EToB.free();
        mesh.EToB = EToB_new;
        EToB_new.free();
        
        mesh.EX.free();
        mesh.EY.free();
        mesh.EX = EX_new;
        mesh.EY = EY_new;
        EX_new.free();
        EY_new.free();

        mesh.Nelements = mesh.Nelements + nn;
        mesh.Nnodes = mesh.Nnodes+new_vertex;
        
        o_PToC = platform.malloc<dlong>(PToC);
        
        o_IntFlag = platform.malloc<dlong>(IntFlag);
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        
        //mesh.SetupUpdate(Nrefine);
        mesh_t updatedMeshPtr = mesh.SetupUpdate(Nrefine);
        mesh = updatedMeshPtr;

        memory<dfloat> Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np);
        Qold = Q;
        deviceMemory<dfloat> o_Qold = platform.malloc<dfloat>(Qold);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Qold ,o_q, o_splitFlag,o_IntFlag,o_PToC,mesh.o_IM);
        SplitFlag.free();
        printf("Refinement Done!\n");
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }//else if(Nrefine==0 || nn==0)
       // {
       //   for (int i = 0; i < 2*mesh.Nelements; ++i)
       //   {
       //     IntFlag[i]=0;
       //   }
       //   o_IntFlag = platform.malloc<dlong>(IntFlag);
       // }
         
}
