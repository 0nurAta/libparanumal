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

#include "bns.hpp"

// Conduct Adaptive Mesh Refinement
void bns_t::Amr(deviceMemory<dfloat>& o_Q,dlong* _N){



    // Construct Refinement Flag
    deviceMemory<dlong> o_refFlag = platform.reserve<dlong>(2*mesh.Nelements);

    // Compute Vorticity for Indicator
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    // Indicator 
    indicatorKernel(mesh.Nelements, o_Vort, o_refFlag);
    
    // Array holds element ids for refining. Holds some extra mem. for 
    // conforming
    //memory<dlong> Ref(mesh.Nelements,0); 
                                    
    // A flag to address that which face will be used for bisection
    memory<dlong> FaceFlag(mesh.Nfaces*mesh.Nelements,0);
    
    memory<dlong> refFlag(2*mesh.Nelements);
    o_refFlag.copyTo(refFlag); // copy data back to host
    
    dlong Nrefine = 0;   
    dlong Ncoarse = 0; 
    
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
    printf("_N_before=%d\n",mesh.Nelements*6*mesh.Np );
    printf("Number_of_Elements_to_be_refined= %d\n",Nrefine);
    printf("Number_of_Elements_to_be_coarsened= %d\n",Ncoarse);

    // copy data back to host
    o_q.copyTo(q);

    // Coarse
    Coarse(q,refFlag,Ncoarse);
    // Conform
    //Conform(refFlag,FaceFlag,Nrefine);
    // Refine
    Refine(q,refFlag,FaceFlag,Nrefine);
    
    *_N = mesh.Nelements*Nfields*mesh.Np;
    printf("_N_after=%d\n",*_N );

    
}
