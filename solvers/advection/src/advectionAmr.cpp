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
    dlong Ncoarse = 0; 
    
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

    printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
    printf("Number_of_Elements_to_be_coarsened= %d\n",Ncoarse);

    // copy data back to host
    o_q.copyTo(q);
    
    //for test
    //Nrefine =1;
    //refFlag[1]=1;
    //Ncoarse =1;
    //refFlag[1]=-1;

    // Coarse
    Coarse(q,refFlag,Ncoarse);
    // copy data back to host
    //o_q.copyTo(q);
    // Conform
    Conform(q,refFlag,Nrefine);
    // Refine
    Refine(q,refFlag,Nrefine);

    
    

    *_N = mesh.Nelements*1*mesh.Np;
}
