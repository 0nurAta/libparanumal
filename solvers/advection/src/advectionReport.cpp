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

void advection_t::Report(dfloat time, int tstep){

  static int frame=0;

  //compute q.M*q
  dlong Nentries = mesh.Nelements*mesh.Np;

   
        
  deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
  printf("New Element Number in report=%d\n",mesh.Nelements);
  mesh.MassMatrixApply(o_q, o_Mq);

  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));



  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    // Refinement Flag
    //deviceMemory<dlong> o_refFlag = platform.reserve<dlong>(2*mesh.Nelements);
    //indicatorKernel(mesh.Nelements, o_q, o_refFlag);
    //memory<dlong> refFlag(2*mesh.Nelements);
    //o_refFlag.copyTo(refFlag); // copy data back to host
    //
    //dlong Nrefine = 0;   
//
    //for (int i = 0; i < mesh.Nelements; ++i)
    //{
    //  if (refFlag[i]==1)
    //  {
    //    Nrefine = Nrefine + 1;
    //  }
    //   
    //}
    //printf("%d\n",Nrefine);

    printf("mesh.Nelements=%d\n",mesh.Nelements );
    // copy data back to host
    o_q.copyTo(q);
    /*for (int i = 0; i < mesh.Nelements*mesh.Np; ++i)
    {
      printf("q1=%f\n",q[i]);
    }*/
    //for test
    //Nrefine =1;
    //refFlag[3]=1;

    //refFlag[1]=1;
    //refFlag[1]=1;
    //refFlag[2]=1;
    //refFlag[5]=1;
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    //Refine(q,refFlag,Nrefine); 
    //memory<dfloat> Qold(2*mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np);
    //Qold = q;
    //deviceMemory<dfloat> o_Qold = platform.malloc<dfloat>(Qold);
    //splitKernel(mesh.Nelements,o_Qold ,o_q, o_refFlag, o_IntFlag,o_PToC,mesh.o_IM);
    // copy data back to host
    //o_q.copyTo(q);
    //Qold = q;
    //deviceMemory<dfloat> o_Qold1 = platform.malloc<dfloat>(Qold);
    //Ncoarse =1;
    //refFlag[3]=-1;
    //Coarse(q,refFlag,Nrefine);
    //printf("Ncoarse=%d\n",Ncoarse);
    //combineKernel(mesh.Nelements,Ncoarse,o_Qold1 ,o_q, o_refFlag, o_IntFlag,o_PToC,mesh.o_RM);

    //  o_q.copyTo(q);
    /*    for (int i = 0; i < mesh.Nelements*mesh.Np; ++i)
    {
      printf("q2=%f\n",q[i]);
    }*/
    PlotFields(q, std::string(fname));

    
  }
}
