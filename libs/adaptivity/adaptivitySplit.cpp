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
//      /|\
//     / | \
//    /  |  \
//   /   |   \
//  /    |    \
// /_____|_____\

void bns_t::Bisect(memory<dlong>& RefFlag,
                         memory<dlong>& FaceFlag,
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

  for (int e = 0; e < mesh.Nelements; ++e)
  {
      
      //int e = Ref[i];
      const dlong id = e*mesh.Nverts; 
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
        if (FaceFlag[id+0]==1)
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
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e]+1;
          

          PToC[e*4+0]   = e;
          PToC[e*4+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*4]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          
          IntFlag[e] = 1;
          IntFlag[(mesh.Nelements+nn)] = 2;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          nn++;
          nv++;
        }

        if (FaceFlag[id+1]==1)
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
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e]+1;

          
          PToC[e*4]   = e;
          PToC[e*4+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*4]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e] = 5;
          IntFlag[(mesh.Nelements+nn)] = 6;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          nn++;
          nv++;
        }

        if (FaceFlag[id+2]==1)
        {
          hlong Local_id = 2+mesh.Nfaces*e+mesh.Nnodes;
          hlong Neigh_id = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;
          hlong newNode = (Local_id>=Neigh_id)? Local_id:Neigh_id ;
          
          EToV_new[id+2] = newNode;
          EToV_new[id_new+0] = newNode;
          EToV_new[id_new+1] = v1;
          EToV_new[id_new+2] = v2;

          EX_new[id+2] = 0.5*(mesh.EX[id+2]+mesh.EX[id+0]);
          EX_new[id_new+0] = EX_new[id+2];
          EX_new[id_new+1] = mesh.EX[id+1];
          EX_new[id_new+2] = mesh.EX[id+2];

          EY_new[id+2] = 0.5*(mesh.EY[id+2]+mesh.EY[id+0]);
          EY_new[id_new+0] = EY_new[id+2];
          EY_new[id_new+1] = mesh.EY[id+1];
          EY_new[id_new+2] = mesh.EY[id+2];

          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = -1;
          EToB_new[id+2] = mesh.EToB[id+2];
          EToB_new[id_new+0] = -1;
          EToB_new[id_new+1] = mesh.EToB[id+1];
          EToB_new[id_new+2] = mesh.EToB[id+2];

          // Update lists related to AMR
          EToRefLevel[e] = EToRefLevel[e]+1;
          EToRefLevel[(mesh.Nelements+nn)] = EToRefLevel[e]+1;

          
          PToC[e*4]   = e;
          PToC[e*4+EToRefLevel[e]] = mesh.Nelements+nn;
          PToC[(mesh.Nelements+nn)*4]   = e;
          //PToC[(mesh.Nelements+nn)*2+1] = mesh.Nelements+nn;
          IntFlag[e] = 3;
          IntFlag[(mesh.Nelements+nn)] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn]=1;

          nn++;
          nv++;
          } 
      }
      *new_vertex = nv;
      *NN = nn;             
  }   
        
}


// Red refinement of triangle i.e. regular refinement
//      /\
//     /  \
//    /____\
//   /\    /\
//  /  \  /  \
// /____\/____\

void bns_t::Red(memory<dlong>& RefFlag,
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
        if (FaceFlag[id+0]==1)
        { 
          // Uniquely number each vertex 
          hlong Local_id_3 = 0+mesh.Nfaces*e+mesh.Nnodes;
          hlong Local_id_4 = 1+mesh.Nfaces*e+mesh.Nnodes;
          hlong Local_id_5 = 2+mesh.Nfaces*e+mesh.Nnodes;

          hlong Neigh_id_3 = mesh.EToF[id+0]+mesh.Nfaces*mesh.EToE[id+0]+mesh.Nnodes;
          hlong Neigh_id_4 = mesh.EToF[id+1]+mesh.Nfaces*mesh.EToE[id+1]+mesh.Nnodes;
          hlong Neigh_id_5 = mesh.EToF[id+2]+mesh.Nfaces*mesh.EToE[id+2]+mesh.Nnodes;

          hlong v3 = (Local_id_3>=Neigh_id_3)? Local_id_3:Neigh_id_3 ;
          hlong v4 = (Local_id_4>=Neigh_id_4)? Local_id_4:Neigh_id_4 ;
          hlong v5 = (Local_id_5>=Neigh_id_5)? Local_id_5:Neigh_id_5 ;

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
          
          IntFlag[e] = 1;
          IntFlag[(mesh.Nelements+nn+0)] = 2;
          IntFlag[(mesh.Nelements+nn+1)] = 3;
          IntFlag[(mesh.Nelements+nn+2)] = 4;

          SplitFlag[e]=1;
          SplitFlag[mesh.Nelements+nn+0]=1;
          SplitFlag[mesh.Nelements+nn+1]=1;
          SplitFlag[mesh.Nelements+nn+2]=1;

          nn=nn+3;
          //new_vertex++;
        }
      }

      *NN = nn;             
  }
}   

// Blue refinement of triangle i.e. double bisect
//      /|\
//     / | \
//    /  |  \
//   /   |  /\
//  /    | /  \
// /_____|/____\


void bns_t::Blue(memory<dlong>& RefFlag,
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