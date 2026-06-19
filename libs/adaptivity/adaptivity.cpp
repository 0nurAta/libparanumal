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
// Conduct Adaptive Mesh Refinement
void adaptivity_t::adaptivity(deviceMemory<dfloat>& o_q,dlong* _N){

    dlong const level = 2;

    // Construct Refinement Flag
    deviceMemory<dlong> o_refFlag = platform.reserve<dlong>(2*mesh.Nelements);

    // Compute Vorticity for Indicator
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    //vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    // Indicator 
    indicatorKernel(mesh.Nelements, o_q, o_refFlag);
    
    // Array holds element ids for refining. Holds some extra mem. for 
    // conforming
    //memory<dlong> Ref(mesh.Nelements,0); 
                                    
    // A flag to address that which face will be used for bisection
    memory<dlong> FaceFlag(4*mesh.Nfaces*mesh.Nelements,0);
    memory<dlong> refFlag(4*mesh.Nelements,0);
    memory<dlong> refFlag2(4*mesh.Nelements,0);
    memory<dlong> confFlag(4*mesh.Nelements,-1);
    memory<dlong> confFlag2(4*mesh.Nelements,-1);
    memory<dlong> coarsefFlag(4*mesh.Nelements,0);
    o_refFlag.copyTo(refFlag); // copy data back to host
    refFlag.copyTo(refFlag2);
    //refFlag.copyTo(confFlag);
    dlong Nrefine = 0;   
    dlong Ncoarse = 0; 
    dlong Nelements_old = mesh.Nelements ;

    memory<dfloat> q(mesh.Np*4*mesh.Nelements,0);
    // copy data back to host
    o_q.copyTo(q);
    //refFlag[18] = 1;
    //refFlag[19] = 1;
    //refFlag[22] = 1;
    //refFlag[21] = 1;
    for (int e = 0; e < mesh.Nelements; ++e)
    {
           /* if (EToRefLevel[e]==5)
      {
        printf("Level 5 e = %d\n",e );
        printf("V0=%d,V1=%d,V2=%d\n",mesh.EToV[e*3+0],mesh.EToV[e*3+1],mesh.EToV[e*3+2]);
        printf("E0=%d,E1=%d,E2=%d\n",mesh.EToE[e*3+0],mesh.EToE[e*3+1],mesh.EToE[e*3+2]);
        printf("B0=%d,B1=%d,B2=%d\n",mesh.EToB[e*3+0],mesh.EToB[e*3+1],mesh.EToB[e*3+2]);
        printf("X0=%f,X1=%f,X2=%f\n",mesh.EX[e*3+0],mesh.EX[e*3+1],mesh.EX[e*3+2]);
        printf("Y0=%f,Y1=%f,Y2=%f\n",mesh.EY[e*3+0],mesh.EY[e*3+1],mesh.EY[e*3+2]);
        printf("field=%f\n",q[e*mesh.Np+0]);

      }*/
      if (refFlag[e]==1)
      {
        //Ref[Nrefine] = e;
        Nrefine = Nrefine + 1;
      }
      if (refFlag[e]==-1)
      {
        Ncoarse = Ncoarse + 1;
      }   
    }
    printf("_N_before=%d\n",mesh.Nelements*1*mesh.Np );
    printf("Nnodes=%lld\n",mesh.Nnodes );
    printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
    printf("Number_of_Elements_to_be_coarsened= %d\n",Ncoarse);

    // Selection of algorithm LE2
    

    // 1 LVL Refinement & Coarsening (Runs with old interpolation method)
    //Coarse(o_q,q,refFlag,Ncoarse);
    //Refine(o_q,q,refFlag,FaceFlag,Nrefine);
    //
    // Conform
    //Conform(refFlag,FaceFlag,Nrefine);
    #if 1
    // RGB Refinement
    memory<dfloat>qold_coarse(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_coarse1(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_coarse2(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_coarse3(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_coarse4(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine1(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine2(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);

    memory<hlong>EToNewV(4*mesh.Nelements*3,-1);
    //RGB
    // Coarse Red refined element
    q.copyTo(qold_coarse);
    dlong const L = EToRefLevel[86];
    dlong e= 86;
    dlong stride= level*4;
        printf(" %d in the main loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    
    //CoarseRed(o_q,q,qold_coarse,refFlag,RedFlag,confFlag,Ncoarse,level);
    // Coarse bisected element

        printf("%d  after the red coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    q.copyTo(qold_coarse1);
    //CoarseGreen(o_q,q,qold_coarse1,refFlag,RedFlag,confFlag,EToNewV,Ncoarse,level);
   
        printf("%d  after the green coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    // Refine a bisected element if it is bisected from its all edges
       q.copyTo(qold_coarse2);
    CoarseGreentoRed(o_q,q,qold_coarse2,refFlag,RedFlag,confFlag,EToNewV,Ncoarse,level);
    q.copyTo(qold_coarse3);
    CoarseGreentoRed(o_q,q,qold_coarse3,refFlag,RedFlag,confFlag,EToNewV,Ncoarse,level);
 
        printf("%d  after the green to red cyc 1 coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    q.copyTo(qold_coarse4);
    CoarseGreentoRed(o_q,q,qold_coarse4,refFlag,RedFlag,confFlag,EToNewV,Ncoarse,level);

        printf("%d  after the green to red cyc 2 coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    q.copyTo(qold_refine1);
    RefineRGB(o_q,q,qold_refine1,refFlag,confFlag,FaceFlag,EToNewV,Nrefine,level);
   
        printf("%d after the refRGB coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
    // Conform by bisecting (Green, Blue refinements)
    q.copyTo(qold_refine2);
    RefineRGB2(o_q,q,qold_refine2,refFlag,confFlag,FaceFlag,Nrefine,Nelements_old,level);
  
        printf("%d  after the refRGB2 coarsening loop with rep=%d and sib=%d",e,PToC[e*(stride)+(L-1)*4+0],PToC[e*(stride)+(L-1)*4+1]);
     #endif

    #if 0
    // Newest Vertex Bisection
    memory<dfloat>qold_coarse(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine1(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine2(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine3(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
      Nrefine =1000;
    q.copyTo(qold_coarse);
    // Rivara 2nd Algorithm (Global)
    CoarsebyID(o_q,q,qold_coarse,refFlag,Ncoarse,level);
    q.copyTo(qold_refine1);
    //printf("FIRST STAGE STARTS\n");
    RefinebyBisect(o_q,q,qold_refine1,refFlag,confFlag,FaceFlag,Nrefine,level);
    //printf("SECOND STAGE STARTS\n");
    q.copyTo(qold_refine2);
    int i = 0;
    while (i<5){
    RefinebyNV(o_q,q,qold_refine2,refFlag,confFlag,FaceFlag,Nrefine,Nelements_old,level);
    i++;
    }
    #endif
    #if 0
    // Algorithm 1 (Rivara,1994)
    memory<dfloat>qold_coarse(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine1(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine2(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
   Nrefine =1000;
    q.copyTo(qold_coarse);
    // Rivara 2nd Algorithm (Global)
    CoarsebyID(o_q,q,qold_coarse,refFlag,Ncoarse,level);
    q.copyTo(qold_refine1);
    printf("FIRST STAGE STARTS\n");
    RefinebyBisect(o_q,q,qold_refine1,refFlag,confFlag,FaceFlag,Nrefine,level);
    printf("SECOND STAGE STARTS\n");
    q.copyTo(qold_refine2);
    int i = 0;
    while (i<20){
    RefinebyLE(o_q,q,qold_refine2,refFlag,confFlag,FaceFlag,Nrefine,Nelements_old,level);
    i++;
    }
    #endif
      #if 0
    // Algorithm 2 (Rivara,1994)
   memory<dfloat>qold_coarse(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
   memory<dfloat>qold_refine1(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
   memory<dfloat>qold_refine2(4*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
   q.copyTo(qold_coarse);
   // Rivara 2nd Algorithm (Global)
   CoarsebyID(o_q,q,qold_coarse,refFlag,Ncoarse,level);
   q.copyTo(qold_refine1);
   //printf("FIRST STAGE STARTS\n");
   RefinebyBisect(o_q,q,qold_refine1,refFlag,confFlag,FaceFlag,Nrefine,level);
   //printf("SECOND STAGE STARTS\n");
   q.copyTo(qold_refine2);
   RefinebyID2(o_q,q,qold_refine2,refFlag,confFlag,FaceFlag,Nrefine,Nelements_old,level);
    #endif
    //printf("q[0]=%d,q[1]=%d,q[2]=%d\n",mesh.EToE[686*mesh.Nverts+0],mesh.EToE[686*mesh.Nverts+1],mesh.EToE[686*mesh.Nverts+2] );
    *_N = mesh.Nelements*1*mesh.Np;

    printf("_N_after=%d\n",*_N );
    printf("Number_of_Elements_after= %d\n",mesh.Nelements);
    
}
}
