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


void bns_t::Refine(memory<dfloat>& Q,
                         memory<dlong>& RefFlag,
                         memory<dlong>& FaceFlag,
                         dlong Nrefine){

  
  //dlong const MAX_REFINEMENT_LEVEL = 1;

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(2*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    for (int n = 0; n < 3; ++n)
    {
    const dlong id = e*mesh.Nverts+n;
    EToB_new[id] = mesh.EToB[id];
    EToV_new[id] = mesh.EToV[id];
    EX_new[id] = mesh.EX[id]; 
    EY_new[id] = mesh.EY[id]; 
    }
  }
  
  hlong new_vertex = 0; // Counts each new_vertex that will be created
  
  // Determine elements to be refined by using Refine Flag
  
  LongestEdge(FaceFlag,RefFlag);
  Conform(RefFlag,FaceFlag,Nrefine);

  // Refinement Loop
  // Determine ids of new vertices and EToV

  hlong nn = 0 ; // Counts each refinement
  
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);

  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn);
  Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
    
      if (Nrefine!=0 && nn!=0)
      {
        // Update mesh connectivity and physical coordinates
        mesh.EToV = EToV_new;
        mesh.EToB = EToB_new;
        mesh.EX = EX_new;
        mesh.EY = EY_new;

        // Update total number of elements and nodes
        mesh.Nelements = mesh.Nelements + nn;
        //mesh.Nnodes = mesh.Nnodes + new_vertex;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);
        mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,mesh.o_IM);

        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }


        
}
