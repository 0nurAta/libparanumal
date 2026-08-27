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
    deviceMemory<dlong> o_refFlag = platform.reserve<dlong>(mesh.Nelements);

    // Compute Vorticity for Indicator
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    //vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    // Indicator 
    indicatorKernel(mesh.Nelements, o_q, o_refFlag);    
                                 
    // A flag to address that which face will be used for bisection
    deviceMemory<dlong> o_FaceFlag = platform.reserve<dlong>(16*mesh.Nfaces*mesh.Nelements); 
    deviceMemory<dlong> o_confFlag = platform.reserve<dlong>(16*mesh.Nfaces*mesh.Nelements); 
    deviceMemory<dlong> o_coarseFlag = platform.reserve<dlong>(16*mesh.Nfaces*mesh.Nelements); 

    memory<dlong> FaceFlag(16*mesh.Nfaces*mesh.Nelements,0);
    memory<dlong> refFlag(16*mesh.Nelements,0);
    memory<dlong> confFlag(16*mesh.Nelements,-1);
    memory<dlong> coarsefFlag(16*mesh.Nelements,0);
   // o_refFlag[4] = 1;
    o_refFlag.copyTo(refFlag); // copy data back to host

    dlong Nrefine = 0;   
    dlong Ncoarse = 0; 
    dlong Nelements_old = mesh.Nelements ;

    memory<dfloat> q(mesh.Np*128*mesh.Nelements,0);
    // copy data back to host
    o_q.copyTo(q);

    for (int e = 0; e < mesh.Nelements; ++e)
    {
      if (refFlag[e]==1)
      {
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


    std::ofstream nochange_file("RefinebyNVB_nochange.txt", std::ios::app);
    // Newest Vertex Bisection
    memory<dfloat>qold_coarse(128*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine1(128*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine2(128*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
    memory<dfloat>qold_refine3(128*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0);
      Nrefine =3000;
    q.copyTo(qold_coarse);
    //CoarseGPU(o_q,q,qold_coarse,o_refFlag,Ncoarse,level);
    q.copyTo(qold_refine1);
    //printf("FIRST STAGE STARTS\n");
    RefinebyBisectGPU(o_q,q,qold_refine1,o_refFlag,o_confFlag,o_FaceFlag,Nrefine,level);
  

    //printf("SECOND STAGE STARTS\n");
    q.copyTo(qold_refine2);
    int i = 0;
      dlong counter = 0;
    while (i<5){
      dlong elem_before = mesh.Nelements;
    RefinebyNVGPU(o_q,q,qold_refine2,o_refFlag,o_confFlag,o_FaceFlag,Nrefine,level);
    i++;
    dlong elem_after = mesh.Nelements;
    if (elem_after == elem_before && counter==0) {
    counter =1;  
    nochange_file << i << "\n";
    // optional debug print
    printf("%d\n", i);
    }
    }    
    *_N = mesh.Nelements*1*mesh.Np;

    printf("_N_after=%d\n",*_N );
    printf("Number_of_Elements_after= %d\n",mesh.Nelements);
    
}
}
