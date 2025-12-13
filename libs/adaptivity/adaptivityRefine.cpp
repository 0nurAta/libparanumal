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
void adaptivity_t::Refine(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine)
{

  
  //dlong const MAX_REFINEMENT_LEVEL = 1;

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("EToE[26*mesh.Nverts+0]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+0],mesh.EToF[26*mesh.Nverts+0]);
  printf("EToE[26*mesh.Nverts+1]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+1],mesh.EToF[26*mesh.Nverts+1]);
  printf("EToE[26*mesh.Nverts+2]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+2],mesh.EToF[26*mesh.Nverts+2]);

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(2*mesh.Nelements,0);

  // A flag for conforming
  memory<dlong> ConfFlag(2*mesh.Nelements,0);
  RefFlag.copyTo(ConfFlag);
  // Copy old Element to Vertex Connectivity to the New One
  #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    //printf("ConfFlag=%d,  RefFlag=%d\n",ConfFlag[e],RefFlag[e] );
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
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  LongestEdge(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  //ConformByVertex(RefFlag,FaceFlag,Nrefine);
  //LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV

  //hlong nn = 0 ; // Counts each refinement
  
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  BisectNEW(RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
   printf("Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);

     for (int e = 0; e < mesh.Nelements; ++e)
  {
    RefFlag[e] = 0 ;

  }
          // Update mesh connectivity and physical coordinates
        // Update total number of elements and nodes
         mesh.Nelements = mesh.Nelements + nn;
        //mesh.Nnodes = mesh.Nnodes + new_vertex;
        dlong update;
        update = nn;// new_vertex=0;
        nn = 0;
         // for (int e = 0; e < mesh.Nelements; ++e)
         // {
         //   for (int i = 0; i < mesh.Nverts; ++i)
         //   {
         //      mesh.EToV[e*mesh.Nverts+i] = EToV_new[e*mesh.Nverts+i];
         //      mesh.EToB[e*mesh.Nverts+i] = EToB_new[e*mesh.Nverts+i];
         //      mesh.EX  [e*mesh.Nverts+i] = EX_new[e*mesh.Nverts+i];
         //      mesh.EY  [e*mesh.Nverts+i] = EY_new[e*mesh.Nverts+i];
//
         //   }
         // }
          printf("mesh.EToV=%lld\n",mesh.EToV[(mesh.Nelements-1)*mesh.Nverts+2]  );
        //  EToV_new.copyTo(mesh.EToV);
        //  EToB_new.copyTo(mesh.EToB);
        //  EX_new.copyTo(mesh.EX);
        //  EY_new.copyTo(mesh.EY);
  ConformByVertex(RefFlag,ConfFlag,FaceFlag,Nrefine);
  BisectNEW(RefFlag,FaceFlag,ConfFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  printf("Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);

      if (Nrefine!=0 && update!=0)
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
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);
        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }

  // Longest Edge Bisection Loop
  /*dlong Nrefine_old = Nrefine;
  dlong Counter = 1;
  while (Counter) {
    nn=0;
  printf("Counter=%d\n", Counter);
  Nrefine_old = Nrefine;
  ConformLE(RefFlag,FaceFlag,Nrefine);
  LongestEdge(FaceFlag,RefFlag);
  BisectNEW(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  printf("nn=%lld\n", nn);

  Counter = Counter-1;
        if (Nrefine!=0 && nn!=0)
      {
        //Reset nn

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
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);

        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  }*/
        
}


void adaptivity_t::RefineLE(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dlong>& RefFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine)
{

  // Longest Edge Refinement algorithm by Rivara 1984
  // Algorithm 2

  //dlong const MAX_REFINEMENT_LEVEL = 1;
  printf("EToE[26*mesh.Nverts+0]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+0],mesh.EToF[26*mesh.Nverts+0]);
  printf("EToE[26*mesh.Nverts+1]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+1],mesh.EToF[26*mesh.Nverts+1]);
  printf("EToE[26*mesh.Nverts+2]=%lld,EToF[26*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+2],mesh.EToF[26*mesh.Nverts+2]);

  printf("EToB[26*mesh.Nverts+0]=%d\n",mesh.EToB[26*mesh.Nverts+0]);
  printf("EToB[26*mesh.Nverts+1]=%d\n",mesh.EToB[26*mesh.Nverts+1]);
  printf("EToB[26*mesh.Nverts+2]=%d\n",mesh.EToB[26*mesh.Nverts+2]);
  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(2*mesh.Nelements,0);
  // A flag for conforming
  memory<dlong> ConfFlag(2*mesh.Nelements,0);
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
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  LongestEdge(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  ConformLE(RefFlag,FaceFlag,Nrefine);
  LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV

  //hlong nn = 0 ; // Counts each refinement
  
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  BisectNEW(RefFlag,FaceFlag,ConfFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("Bisect Done!\n");
      if (Nrefine!=0 && nn!=0)
      {
        // Update mesh connectivity and physical coordinates
        mesh.EToV = EToV_new;
        mesh.EToB = EToB_new;
        mesh.EX = EX_new;
        mesh.EY = EY_new;

        // Update total number of elements and nodes
        mesh.Nelements = mesh.Nelements + nn;
        mesh.Nnodes = mesh.Nnodes + new_vertex;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);
        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }

  // Longest Edge Bisection Loop
  /*dlong Nrefine_old = Nrefine;
  dlong Counter = 1;
  while (Counter) {
    nn=0;
  printf("Counter=%d\n", Counter);
  Nrefine_old = Nrefine;
  ConformLE(RefFlag,FaceFlag,Nrefine);
  LongestEdge(FaceFlag,RefFlag);
  BisectNEW(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  printf("nn=%lld\n", nn);

  Counter = Counter-1;
        if (Nrefine!=0 && nn!=0)
      {
        //Reset nn

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
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);

        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  }*/
        
}

void adaptivity_t::RefinebyID(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine)
{

  
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

  // A flag for conforming
  //memory<dlong> ConfFlag(2*mesh.Nelements,0);
  //RefFlag.copyTo(ConfFlag);
  // Copy old Element to Vertex Connectivity to the New One
  #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    //printf("ConfFlag=%d,  RefFlag=%d\n",ConfFlag[e],RefFlag[e] );
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
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  LongestEdge(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  //ConformByVertex(RefFlag,FaceFlag,Nrefine);
  //LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV

  //hlong nn = 0 ; // Counts each refinement
  
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  BisectbyID(RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,3);
   printf("Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);

     for (int e = 0; e < mesh.Nelements; ++e)
  {
    RefFlag[e] = 0 ;

  }

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
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);
        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  
    
}


void adaptivity_t::RefinebyID2(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old)
{

  
  //dlong const MAX_REFINEMENT_LEVEL = 1;

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("CONFORMING_STARTS, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  printf("EToE[31*mesh.Nverts+0]=%lld,EToF[31*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+0],mesh.EToF[26*mesh.Nverts+0]);
  printf("EToE[31*mesh.Nverts+1]=%lld,EToF[31*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+1],mesh.EToF[26*mesh.Nverts+1]);
  printf("EToE[31*mesh.Nverts+2]=%lld,EToF[31*mesh.Nverts+0]=%d\n",mesh.EToE[26*mesh.Nverts+2],mesh.EToF[26*mesh.Nverts+2]);
  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(2*mesh.Nelements,0);
  memory<dlong> new_v_id(2*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  Q.copyTo(Qold);
  // A flag for conforming
  //memory<dlong> ConfFlag(2*mesh.Nelements,0);
  //RefFlag.copyTo(ConfFlag);
  // Copy old Element to Vertex Connectivity to the New One
  #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    //printf("ConfFlag=%d,  RefFlag=%d\n",ConfFlag[e],RefFlag[e] );
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
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  //LongestEdge(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  //ConformByVertex(RefFlag,FaceFlag,Nrefine);
  //LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  ConformByBisectMultiLvl(Q,Qold,RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("2nd Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;
  printf("Q_inrefine=%f\n",Q[1690] );
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
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
      
        o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_PToC,o_IM);

        o_q.copyTo(Q);
         printf("Q_inrefine=%f\n",Q[1690] );
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }

  printf("EToE[163*mesh.Nverts+0]=%lld,EToF[163*mesh.Nverts+0]=%d\n",mesh.EToE[163*mesh.Nverts+0],mesh.EToF[163*mesh.Nverts+0]);
  printf("EToE[163*mesh.Nverts+1]=%lld,EToF[163*mesh.Nverts+0]=%d\n",mesh.EToE[163*mesh.Nverts+1],mesh.EToF[163*mesh.Nverts+1]);
  printf("EToE[163*mesh.Nverts+2]=%lld,EToF[163*mesh.Nverts+0]=%d\n",mesh.EToE[163*mesh.Nverts+2],mesh.EToF[163*mesh.Nverts+2]);
  printf("EToV[171*mesh.Nverts+0]=%d,\n",mesh.EToV[171*mesh.Nverts+0],mesh.EToF[171*mesh.Nverts+0]);
  printf("EToV[171*mesh.Nverts+1]=%d,\n",mesh.EToV[171*mesh.Nverts+1],mesh.EToF[171*mesh.Nverts+1]);
  printf("EToV[171*mesh.Nverts+2]=%d,\n",mesh.EToV[171*mesh.Nverts+2],mesh.EToF[171*mesh.Nverts+2]);
  
  printf("EToE[171*mesh.Nverts+0]=%lld,EToF[171*mesh.Nverts+0]=%d\n",mesh.EToE[171*mesh.Nverts+0],mesh.EToF[171*mesh.Nverts+0]);
  printf("EToE[171*mesh.Nverts+1]=%lld,EToF[171*mesh.Nverts+0]=%d\n",mesh.EToE[171*mesh.Nverts+1],mesh.EToF[171*mesh.Nverts+1]);
  printf("EToE[171*mesh.Nverts+2]=%lld,EToF[171*mesh.Nverts+0]=%d\n",mesh.EToE[171*mesh.Nverts+2],mesh.EToF[171*mesh.Nverts+2]);
  printf("EX[164*mesh.Nverts+0]=%f,EY[164*mesh.Nverts+0]=%f\n",mesh.EX[171*mesh.Nverts+0],mesh.EY[171*mesh.Nverts+0]);
  printf("EX[164*mesh.Nverts+1]=%f,EY[164*mesh.Nverts+0]=%f\n",mesh.EX[171*mesh.Nverts+1],mesh.EY[171*mesh.Nverts+1]);
  printf("EX[164*mesh.Nverts+2]=%f,EY[164*mesh.Nverts+0]=%f\n",mesh.EX[171*mesh.Nverts+2],mesh.EY[171*mesh.Nverts+2]);
 }
}