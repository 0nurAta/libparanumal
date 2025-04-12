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

// Conduct Adaptive Mesh Refinement
void advection_t::Amr(deviceMemory<dfloat>& o_Q,dlong* _N){



    // Construct Refinement Flag
    deviceMemory<dlong> o_refFlag = platform.reserve<dlong>(2*mesh.Nelements);

    // Indicator 
    indicatorKernel(mesh.Nelements, o_q, o_refFlag);
    
    memory<dlong> refFlag(2*mesh.Nelements);
    o_refFlag.copyTo(refFlag); // copy data back to host
    
    dlong Nrefine = 0;   
    Ncoarse = 0; 
    
    for (int i = 0; i < mesh.Nelements; ++i)
    {
      if (refFlag[i]==1)
      {
        Nrefine = Nrefine + 1;
      }   
    }
    printf("%d\n",Nrefine);

    //    for (int i = 0; i < mesh.Nelements; ++i)
    //{
    //  if (refFlag[i]==-1)
    //  {
    //    Ncoarse = Ncoarse + 1;
    //  }
    //   
    //}
    printf("Ncoarse_outside=%d\n",Ncoarse);

        // Store Interpolated Solution
    o_q.copyTo(q);
    //Qold = q;
    //deviceMemory<dfloat> o_Qold1 = platform.malloc<dfloat>(Qold);

    //Ncoarse =2;
    //refFlag[1]=-1;
    //refFlag[2]=-1;
    //refFlag[34]=-1;
    //refFlag[40]=-1;

    //refFlag[6]=-1;

    // Coarse
    Coarse(q,refFlag,Nrefine);
    // copy data back to host
    o_q.copyTo(q);

    //for test
    //Nrefine =3;
    //refFlag[1]=1;
    //refFlag[2]=1;
    //refFlag[34]=1;
    //refFlag[37]=1;
    //refFlag[40]=1;
    //refFlag[126]=1;
    //refFlag[23]=1;
    //refFlag[232]=1;
    //refFlag[246]=1;
    //refFlag[85]=1;
    //refFlag[77]=1;
    // Refine
    Refine(q,refFlag,Nrefine); 

    // Store old solution
//    memory<dfloat> Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np);
//    Qold = q;
//    deviceMemory<dfloat> o_Qold = platform.malloc<dfloat>(Qold);
//
//    // Interpolate Solution
//    splitKernel(mesh.Nelements,o_Qold ,o_q, o_refFlag, o_IntFlag,o_PToC,mesh.o_IM);



    //printf("Ncoarse=%d\n",Ncoarse);

    // Restrict Solution
    //combineKernel(mesh.Nelements,Ncoarse,o_Qold1 ,o_q, o_refFlag, o_IntFlag,o_PToC,mesh.o_RM);

    //o_q.copyTo(q);
    /*    for (int i = 0; i < mesh.Nelements*mesh.Np; ++i)
    {
      printf("q2=%f\n",q[i]);
    }*/
    //PlotFields(q, refFlag,std::string(fname));

    *_N = mesh.Nelements*1*mesh.Np;
    //printf("N=%d\n",*_N );
  
}
