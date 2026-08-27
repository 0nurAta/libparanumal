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

void adaptivity_t::RefinebyBisect(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong level)
{

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(16*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(16*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(16*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(16*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(16*mesh.Nelements,0);

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

  hlong nn = 0 ; // Counts each refinement
  
  //LongestEdge(FaceFlag,RefFlag,level);
  memory<hlong> new_v_id(16*Nrefine,-1);
  LongestEdge2(FaceFlag,RefFlag,level,Nrefine,new_v_id,&new_vertex);
 // NewestVertex2(FaceFlag,RefFlag,level,Nrefine,new_v_id,&new_vertex);
  //dlong* _Nrefine = &Nrefine;
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  const dlong conflevel = level;
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,level);
   
  BisectNew(Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,new_v_id,Nrefine,&nn,&new_vertex,level,conflevel);
   printf("Bisect Done! Nrefine=%d, nn=%lld, new_vertex=%lld\n", Nrefine,nn,new_vertex);

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
       mesh.Nnodes = mesh.Nnodes + new_vertex;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);

        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        //deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        //deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        printf("Nodes=%d\n",mesh.Nnodes );
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);
        o_q.copyFrom(Q);

        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  
    
}


void adaptivity_t::RefinebyBisectGPU(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          deviceMemory<dlong>& o_RefFlag,
                          deviceMemory<dlong>& o_ConfFlag,
                          deviceMemory<dlong>& o_FaceFlag,
                          dlong Nrefine,
                          dlong level)
{

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  
  memory<long long int>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);

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

  hlong nn = 2 ; // Counts each refinement
  
  // 
  deviceMemory<dfloat> o_qold = platform.malloc<dfloat>(Q);
  deviceMemory<long long int> o_new_v_id = platform.reserve<long long int>(16*Nrefine);
  deviceMemory<long long int> o_elemList = platform.reserve<long long int>(16*Nrefine);
  deviceMemory<long long int> o_elemListCompact = platform.reserve<long long int>(16*Nrefine);
  deviceMemory<dlong> o_NrefineOut = platform.reserve<dlong>(1);
  deviceMemory<dlong> o_new_vertex = platform.reserve<dlong>(1);
  deviceMemory<dlong> o_activeFlag = platform.reserve<dlong>(mesh.Nelements);
  deviceMemory<dlong> o_PCS = platform.reserve<dlong>(128*mesh.Nelements*3);
  deviceMemory<long long int> o_EToV_new = platform.reserve<long long int>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dlong> o_EToB_new = platform.reserve<int>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dlong> o_SplitFlag = platform.reserve<dlong>(16*mesh.Nelements);
  o_EToE = platform.malloc<long long int>(mesh.EToE);
  o_EToV_new = platform.malloc<long long int>(EToV_new);
  o_EToV = platform.malloc<long long int>(EToV_new);

  o_EToB_new = platform.malloc<int>(EToB_new);
  o_EToRefLevel = platform.malloc<dlong>(EToRefLevel);
  o_EX = platform.malloc<dfloat>(mesh.EX);
  o_EY = platform.malloc<dfloat>(mesh.EY);
  deviceMemory<dfloat> o_EX_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dfloat> o_EY_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  o_EX_new = platform.malloc<dfloat>(EX_new);
  o_EY_new = platform.malloc<dfloat>(EY_new);
  dlong Nnodes = mesh.Nnodes;
  const dlong conflevel = level;
  candidateKernel(mesh.Nelements,level,o_FaceFlag,o_RefFlag,
                  o_EToRefLevel,o_EX,o_EY,o_EToV,o_elemList,
                  o_activeFlag);
                  printf("candidateKernel is completed\n");
  assignKernel(mesh.Nelements,Nnodes,o_elemList,o_activeFlag,
               o_elemListCompact,o_new_v_id,o_NrefineOut,o_new_vertex);
                  printf("assignKernel is completed\n");
  memory<dlong>NrefineActual(1,0);
memory<dlong>newVertexActual(1,0);
  
  o_NrefineOut.copyTo(NrefineActual);
  o_new_vertex.copyTo(newVertexActual);              

  const dlong actualNrefine = NrefineActual[0];
  const dlong newVertexCount = newVertexActual[0];  


  

  bisectKernel(actualNrefine,mesh.Nelements,level,conflevel,
               o_q,o_qold,o_RefFlag,o_FaceFlag,o_ConfFlag,
               o_EX,o_EY,o_EToV,mesh.o_EToB, o_EToE,
               o_EX_new,o_EY_new,o_EToV_new,o_EToB_new,
               o_SplitFlag,o_new_v_id, 
               o_EToRefLevel,o_PCS,o_PToC,o_IntFlag,o_IM);
                  printf("bisectKernel is completed\n");

                  memory<dlong> checkLevel(EToRefLevel.length(), 0);
memory<dlong> checkInt(IntFlag.length(), 0);

o_EToRefLevel.copyTo(checkLevel);
o_IntFlag.copyTo(checkInt);

for(dlong e=0; e<mesh.Nelements+actualNrefine; ++e) {
  const dlong L = checkLevel[e];

  if(L > 0) {
    printf("e=%d level=%d rule=%d\n",
           e,
           L,
           checkInt[e*(level+3)+L-1]);
  }
}
  printf("Number_of_Elements_refined= %d\n",NrefineActual[0]);
  printf("Number_of_Vertices_created= %d\n",newVertexActual[0]);




  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  
  printf("Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
 
  /*BisectNew(Q,Qold,RefFlag,FaceFlag, ConfFlag, 
              EX_new,EY_new,EToV_new,EToB_new,
              SplitFlag,new_v_id,Nrefine,&nn,
              &new_vertex,level,conflevel);*/
  
  printf("Bisect Done! Nrefine=%d, nn=%d, new_vertex=%d\n", actualNrefine,newVertexCount,new_vertex);

 /*    for (int e = 0; e < mesh.Nelements; ++e)
  {
    RefFlag[e] = 0 ;
 
  }*/

      if (actualNrefine!=0)
      {
        // Update mesh connectivity and physical coordinates
        o_EToV_new.copyTo(EToV_new); 
        o_EToB_new.copyTo(EToB_new); 
        o_EX_new.copyTo(EX_new);
        o_EY_new.copyTo(EY_new);
        o_EToRefLevel.copyTo(EToRefLevel);
        mesh.EToV = EToV_new;
        
        mesh.EToB = EToB_new;
        mesh.EX = EX_new;
        mesh.EY = EY_new;
        
        // Update total number of elements and nodes
        mesh.Nelements = mesh.Nelements + actualNrefine;
        mesh.Nnodes = mesh.Nnodes + newVertexCount;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);

        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        //o_PToC = platform.malloc<dlong>(PToC);  
      //  o_IntFlag = platform.malloc<dlong>(IntFlag);   
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        //deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        //deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        printf("Nodes=%d\n",mesh.Nnodes );
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);
        //o_q.copyFrom(Q);

        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  
    
}


void adaptivity_t::RefinebyID2(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old,
                          dlong level)
{

  

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  memory<hlong>EToV_new(4*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(4*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(4*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(4*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(4*mesh.Nelements,0);
  memory<dlong> new_v_id(4*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);

   printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
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
  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  //hlong nn = 0 ; // Counts each refinement
  ConformByBisectMultiLvl(Q,Qold,RefFlag,ConfFlag,FaceFlag,
    new_v_id,Nrefine,
    EX_new,EY_new,EToV_new,EToB_new,
    SplitFlag,
    &nn,&new_vertex,level);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("2nd Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;

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
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        o_q.copyFrom(Q);
        //deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        //deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        //o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }


 void adaptivity_t::RefinebyLE(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old,
                          dlong level)
{

  
  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyLE starts, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  memory<hlong>EToV_new(32*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(32*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(32*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(32*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  printf("Nrefine in Conform = %d\n", Nrefine );
  memory<dlong> SplitFlag(32*mesh.Nelements,0);
  memory<hlong> new_v_id(32*Nrefine,-1);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);

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
  dlong const reflevel= level;
  dlong const conflevel= level+3;
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  printf("ConfFlag Last elem =%d \n",ConfFlag[mesh.Nelements-1]);
   LongestEdgeConform(FaceFlag,RefFlag,ConfFlag,level,Nrefine,new_v_id,&new_vertex);
   BisectNew(Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,new_v_id,Nrefine,&nn,&new_vertex,reflevel,conflevel);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  
  //dlong const count = mesh.Nelements+nn-Nelements_old;

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
         printf("2nd Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }

 void adaptivity_t::RefinebyNV(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old,
                          dlong level)
{

  
  

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);
  //memory<dlong> new_v_id(2*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);
  memory<hlong> new_v_id(8*Nrefine,-1);
   printf("RefinebyID3_Starts, number of total elements=%d\n",mesh.Nelements);
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
    dlong const reflevel= level;
  dlong const conflevel= level+3;
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  NewestVertexConform(FaceFlag,RefFlag,ConfFlag,level,Nrefine,new_v_id,&new_vertex);
  BisectNew(Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,new_v_id,Nrefine,&nn,&new_vertex,reflevel,conflevel);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("2nd Bisect Done! Nrefine=%d, nn=%lld\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;

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
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }

 void adaptivity_t::RefinebyNVGPU(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          deviceMemory<dlong>& o_RefFlag,
                          deviceMemory<dlong>& o_ConfFlag,
                          deviceMemory<dlong>& o_FaceFlag,
                          dlong Nrefine,
                          dlong level)
{

  
  

  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
//for (dlong e = 0; e < mesh.Nelements; ++e){
//  dlong id =  e*mesh.Nfaces;
//  FaceFlag[id+0] = 0;
//  FaceFlag[id+1]=0;
//  FaceFlag[id+2]=0;}
  memory<long long int>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);
  //memory<dlong> new_v_id(2*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);
  memory<hlong> new_v_id(8*Nrefine,-1);
   printf("RefinebyID3_Starts, number of total elements=%d\n",mesh.Nelements);
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
    dlong const reflevel= level;
  dlong const conflevel= level+3;
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  //NewestVertexConform(FaceFlag,RefFlag,ConfFlag,level,Nrefine,new_v_id,&new_vertex);
  //BisectNew(Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,new_v_id,Nrefine,&nn,&new_vertex,reflevel,conflevel);
  
  // Declare and define kernel inputs
   const dlong Nelements = mesh.Nelements;
   const dlong Nfaces = mesh.Nfaces;
   dlong Nnodes = mesh.Nnodes;
  deviceMemory<dfloat> o_qold = platform.malloc<dfloat>(Q);
  deviceMemory<long long int> o_new_v_id = platform.reserve<long long int>(2*Nelements);
  deviceMemory<long long int> o_elemList = platform.reserve<long long int>(3*Nelements);
  deviceMemory<long long int> o_elemListCompact = platform.reserve<long long int>(3*Nelements);
  deviceMemory<dlong> o_confFace = platform.reserve<dlong>(Nelements);
  deviceMemory<dlong> o_splitFace = platform.reserve<dlong>(Nelements);
  deviceMemory<long long int> o_hangVertex = platform.reserve<long long int>(Nelements);
  deviceMemory<dlong> o_activeFlag = platform.reserve<dlong>(Nelements);
  deviceMemory<dlong> o_SplitFlag = platform.reserve<dlong>(16*mesh.Nelements);

  // NewestVertexConform expects fresh face flags.
  memory<dlong> FaceFlagZero(Nfaces*Nelements, 0);
  o_FaceFlag.copyFrom(FaceFlagZero, Nfaces*Nelements);

  deviceMemory<dlong> o_PCS = platform.malloc<dlong>(PCS);
  deviceMemory<int> o_EToB_new = platform.reserve<int>(16*mesh.Nelements*mesh.Nverts);
 // o_PToC = platform.malloc<dlong>(PToC);

  o_EToRefLevel = platform.malloc<dlong>(EToRefLevel);


  // Use the current mesh data.
  o_EX = platform.malloc<dfloat>(mesh.EX);
  o_EY = platform.malloc<dfloat>(mesh.EY);
  deviceMemory<dfloat> o_EX_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dfloat> o_EY_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<long long int> o_EToV_new = platform.reserve<long long int>(16*mesh.Nelements*mesh.Nverts);
  o_EX_new = platform.malloc<dfloat>(EX_new);
  o_EY_new = platform.malloc<dfloat>(EY_new);
  o_EToE = platform.malloc<long long int>(mesh.EToE);
  o_EToV_new = platform.malloc<long long int>(EToV_new);
  o_EToV = platform.malloc<long long int>(EToV_new);
  o_EToB_new = platform.malloc<int>(EToB_new);
  o_EToRefLevel = platform.malloc<dlong>(EToRefLevel);
  deviceMemory<int> o_EToF = platform.malloc<int>(mesh.EToF);

  // Compact workspaces
  deviceMemory<dlong> o_confFaceCompact = platform.reserve<dlong>(mesh.Nelements);
  deviceMemory<dlong> o_splitFaceCompact = platform.reserve<dlong>(mesh.Nelements);
  deviceMemory<long long> o_hangVertexCompact = platform.reserve<long long>(mesh.Nelements);
  deviceMemory<dlong> o_NrefineOut = platform.reserve<dlong>(1);
  deviceMemory<dlong> o_new_vertex = platform.reserve<dlong>(1);
  

  conformKernel(Nelements,level,o_RefFlag,o_ConfFlag,
                o_FaceFlag,o_EToRefLevel,o_PCS,o_PToC,
                o_EToF,mesh.o_EToB,o_EX,o_EY,o_EToV,
                o_elemList,o_confFace,o_splitFace,o_hangVertex,
                o_activeFlag);



  assignconformKernel(mesh.Nelements,Nnodes,o_elemList,o_activeFlag,
                      o_confFace,o_splitFace,o_hangVertex,o_elemListCompact,
                      o_confFaceCompact,o_splitFaceCompact,o_hangVertexCompact,o_new_v_id,
                      o_NrefineOut,o_new_vertex);

  memory<dlong>NrefineActual(1,0);
  memory<dlong>newVertexActual(1,0);
  
  o_NrefineOut.copyTo(NrefineActual);
  o_new_vertex.copyTo(newVertexActual);              

  const dlong actualNrefine = NrefineActual[0];
  const dlong newVertexCount = newVertexActual[0];  
                    
  
  bisectKernel(actualNrefine,mesh.Nelements,level,conflevel,
               o_q,o_qold,o_RefFlag,o_FaceFlag,o_ConfFlag,
               o_EX,o_EY,o_EToV,mesh.o_EToB, o_EToE,
               o_EX_new,o_EY_new,o_EToV_new,o_EToB_new,
               o_SplitFlag,o_new_v_id, 
               o_EToRefLevel,o_PCS,o_PToC,o_IntFlag,o_IM);

           
memory<dlong> checkLevel(EToRefLevel.length(), 0);
memory<dlong> checkInt(IntFlag.length(), 0);

o_EToRefLevel.copyTo(checkLevel);
o_IntFlag.copyTo(checkInt);

for(dlong e=0; e<mesh.Nelements+actualNrefine; ++e) {
  const dlong L = checkLevel[e];

  if(L > 0) {
    printf("e=%d level=%d rule=%d\n",
           e,
           L,
           checkInt[e*(level+3)+L-1]);
  }
}
  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("2nd Bisect Done! Nrefine=%d, nn=%lld\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;

      if (Nrefine!=0 && actualNrefine > 0)
      {
        // Update mesh connectivity and physical coordinates
        o_EToV_new.copyTo(EToV_new); 
        o_EToB_new.copyTo(EToB_new); 
        o_EX_new.copyTo(EX_new);
        o_EY_new.copyTo(EY_new);
        o_EToRefLevel.copyTo(EToRefLevel);
        mesh.EToV = EToV_new;
        
        mesh.EToB = EToB_new;
        mesh.EX = EX_new;
        mesh.EY = EY_new;
        printf("mesh.EToV=%d,EToV_new=%d\n",mesh.EToV[4*3+1],EToV_new[4*3+1]);
        // Update total number of elements and nodes
        mesh.Nelements = mesh.Nelements + actualNrefine;
        mesh.Nnodes = mesh.Nnodes + newVertexCount;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
      //  o_PToC = platform.malloc<dlong>(PToC);  
//        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        //o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }


void adaptivity_t::RefineRGB(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& confGreen,
                          memory<dlong>& FaceFlag,
                          memory<hlong>& EToNewV,
                          dlong Nrefine,
                          dlong level)
{

  


  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);

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
  memory<hlong> new_v_id(8*mesh.Nelements*3,-1);
  memory<dlong> ConfRed(8*mesh.Nelements*3,-1);
  // Determine elements to be refined by using Refine Flag
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  //LongestEdge(FaceFlag,RefFlag);
  //NewestVertex(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  //ConformByVertex(RefFlag,FaceFlag,Nrefine);
  //LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV

  //hlong nn = 0 ; // Counts each refinement
  
  printf("Red Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  RGB_flag(FaceFlag, RefFlag, level, Nrefine, new_v_id, ConfRed,confGreen,EToNewV,&new_vertex);
  Red(o_q,Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,new_v_id,Nrefine,&nn,&new_vertex,level);
  //RedOld(o_q,Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,new_v_id,&nn,&new_vertex,level);
   printf("Red refinement Done! Nrefine=%d, nn=%d\n", Nrefine,nn);

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
        mesh.Nnodes = mesh.Nnodes + new_vertex;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);
        o_q.copyFrom(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  
    
}

void adaptivity_t::RefineRGB2(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& RefFlag2,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old,
                          dlong level)
{

  


  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);
  memory<dlong> new_v_id(8*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);

   printf("RefineRGB2 starts, number of total elements=%d\n",mesh.Nelements);
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
  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  
  // Green and Blue bisection to conform Red Refined Elements
  ConformByBisectGB(Q,Qold,RefFlag,RefFlag2,ConfFlag,FaceFlag,new_v_id,Nrefine,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,&nn,&new_vertex,level);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("2nd Bisect Done! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;

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
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }

void adaptivity_t::RefineRGB3(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          dlong Nrefine,
                          dlong Nelements_old,
                          dlong level)
{

  


  // Store old info & Allocate new arrays
  // Element to vertex & Element to boundary connectivity 
  printf("RefinebyID2_Starts, number of total elements=%d\n",mesh.Nelements);
for (dlong e = 0; e < mesh.Nelements; ++e){
  dlong id =  e*mesh.Nfaces;
  FaceFlag[id+0] = 0;
  FaceFlag[id+1]=0;
  FaceFlag[id+2]=0;}
  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);
  
  // A flag to be used in split kernel, initialized with zero values.
  memory<dlong> SplitFlag(8*mesh.Nelements,0);
  memory<dlong> new_v_id(8*mesh.Nelements*mesh.Nverts,0);
  // For local interpolation
  //memory<dfloat>Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
  //Q.copyTo(Qold);

   printf("RefineRGB3 starts, number of total elements=%d\n",mesh.Nelements);
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
  
  // Refinement Loop
  // Determine ids of new vertices and EToV
  // printf("Q_inrefinebefore=%f\n",Q[1621] );
  //hlong nn = 0 ; // Counts each refinement
  
  // Green and Blue bisection to conform Red Refined Elements
  ConformByBisectRToR(Q,Qold,RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,&nn,&new_vertex,level);

  //ConformByID(RefFlag,ConfFlag,FaceFlag,new_v_id,Nrefine);
  //printf("Bisect Start! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //BisectbyID2(RefFlag,FaceFlag,ConfFlag,new_v_id,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,1);
  printf("RefineRGB3 Done! Nrefine=%d, nn=%d\n", Nrefine,nn);
  //dlong const count = mesh.Nelements+nn-Nelements_old;

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
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        o_q.copyFrom(Q);
        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
       
        // Interpolate Solution
        splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);

        o_q.copyTo(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
 }



void adaptivity_t::RefineGreentoRed(deviceMemory<dfloat>& o_q,
                          memory<dfloat>& Q,
                          memory<dfloat>& Qold,
                          memory<dlong>& RefFlag,
                          memory<dlong>& ConfFlag,
                          memory<dlong>& FaceFlag,
                          memory<hlong>& EToNewV,
                          dlong Nrefine,
                          dlong level)
{

  


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
  memory<hlong> new_v_id(4*mesh.Nelements*3,-1);
  memory<dlong> ConfRed(4*mesh.Nelements*3,-1);
  // Determine elements to be refined by using Refine Flag
  /*for (dlong i = 0; i < mesh.Nelements*mesh.Nfaces; ++i){
  FaceFlag[i] = 0;}*/

   hlong nn = 0 ; // Counts each refinement
  //LongestEdge(FaceFlag,RefFlag);
  //NewestVertex(FaceFlag,RefFlag);
  //dlong* _Nrefine = &Nrefine;
  //ConformByVertex(RefFlag,FaceFlag,Nrefine);
  //LongestEdge(FaceFlag,RefFlag);
  //Nrefine* = _Nrefine;
  printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
  
  // Refinement Loop
  // Determine ids of new vertices and EToV

  //hlong nn = 0 ; // Counts each refinement
  
  printf("Red Refinement Start!\n");
  printf("Old Element Number=%d\n",mesh.Nelements);
  //Bisect(RefFlag,FaceFlag,EX_new,EY_new,EToV_new,EToB_new,SplitFlag,&nn,&new_vertex,2);
  //RGB_flag(FaceFlag, RefFlag, level, Nrefine, new_v_id, ConfRed,&new_vertex);
  Red(o_q,Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,EToNewV,Nrefine,&nn,&new_vertex,level);
  RedOld(o_q,Q,Qold,RefFlag,FaceFlag, ConfFlag, EX_new,EY_new,EToV_new,EToB_new,SplitFlag,RedFlag,new_v_id,&nn,&new_vertex,level);
   printf("Red refinement Done! Nrefine=%d, nn=%d\n", Nrefine,nn);

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
        mesh.Nnodes = mesh.Nnodes + new_vertex;
        
        // Update mesh
        mesh = mesh.SetupUpdate(Nrefine);
        // mesh.PmlSetup();
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);  // NEW!!
        o_PToC = platform.malloc<dlong>(PToC);  
        o_IntFlag = platform.malloc<dlong>(IntFlag);   
        o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        deviceMemory<dlong> o_splitFlag = platform.malloc<dlong>(SplitFlag);
        
        // Interpolate Solution
        //splitKernel(mesh.Nelements,o_Q ,o_q, o_splitFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_IM,level);
        o_q.copyFrom(Q);
        printf("Refinement Done!, Nrefine=%d\n",Nrefine);
        printf("new_vertex_count=%lld\n",new_vertex);
        printf("New Element Number=%d\n",mesh.Nelements);
      }
  
    
}



}