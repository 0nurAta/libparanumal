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
    memory<dlong> FaceFlag(2*mesh.Nfaces*mesh.Nelements,0);
    memory<dlong> refFlag(2*mesh.Nelements,0);
    memory<dlong> confFlag(2*mesh.Nelements,0);
    memory<dlong> coarsefFlag(2*mesh.Nelements,0);
    o_refFlag.copyTo(refFlag); // copy data back to host
    //refFlag.copyTo(confFlag);
    dlong Nrefine = 0;   
    dlong Ncoarse = 0; 
    dlong Nelements_old = mesh.Nelements ;

    for (int e = 0; e < mesh.Nelements; ++e)
    {
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
    printf("Nnodes=%d\n",mesh.Nnodes );
    printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
    printf("Number_of_Elements_to_be_coarsened= %d\n",Ncoarse);
    memory<dfloat> q(mesh.Np*2*mesh.Nelements,0);
    // copy data back to host
    o_q.copyTo(q);

    // 1 LVL Refinement & Coarsening (Runs with old interpolation method)
    //Coarse(o_q,q,refFlag,Ncoarse);
    //Refine(o_q,q,refFlag,FaceFlag,Nrefine);
    //

    // Conform
    //Conform(refFlag,FaceFlag,Nrefine);
    // Refine
    memory<dfloat>qold_coarse(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine1(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine2(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    q.copyTo(qold_coarse);
    // Rivara 2nd Algorithm (Global)
    CoarsebyID(o_q,q,qold_coarse,refFlag,Ncoarse,level);
    q.copyTo(qold_refine1);
    //printf("FIRST STAGE STARTS\n");
    RefinebyID(o_q,q,qold_refine1,refFlag,confFlag,FaceFlag,Nrefine);
    //printf("SECOND STAGE STARTS\n");
    q.copyTo(qold_refine2);
    RefinebyID2(o_q,q,qold_refine2,refFlag,confFlag,FaceFlag,Nrefine,Nelements_old);
    //
    
    *_N = mesh.Nelements*1*mesh.Np;
            printf("x1=%f,x2=%f,x3=%f\n",mesh.EX[486*mesh.Nverts+0],mesh.EX[486*mesh.Nverts+1],mesh.EX[486*mesh.Nverts+2]);
        printf("y1=%f,y2=%f,y3=%f\n",mesh.EY[486*mesh.Nverts+0],mesh.EY[486*mesh.Nverts+1],mesh.EY[486*mesh.Nverts+2]);
    printf("_N_after=%d\n",*_N );
    printf("Number_of_Elements_after= %d\n",mesh.Nelements);
    
}
}