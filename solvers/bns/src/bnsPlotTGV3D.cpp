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

// interpolate data to plot nodes and save to file (one per process)
void bns_t::PlotTGV3D(memory<dfloat>& Q, memory<dfloat>& V, std::string fileName, dfloat time){

  FILE *fp;
  dfloat KE=0.0;
  dfloat Eps=0.0;

  if(mesh.rank==0)
    fp = fopen("KE.txt", "a+");

  //fprintf(fp, "\n Kinetic Energy>\n");

  memory<dfloat> u(mesh.Np);
  memory<dfloat> v(mesh.Np);
  memory<dfloat> w(mesh.Np);

    // write out velocity
    
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        const dfloat rm = Q[e*mesh.Np*Nfields+n];
        u[n] = c*Q[e*mesh.Np*Nfields+n+mesh.Np*1]/rm;
        v[n] = c*Q[e*mesh.Np*Nfields+n+mesh.Np*2]/rm;
        if(mesh.dim==3){
          w[n] = c*Q[e*mesh.Np*Nfields+n+mesh.Np*3]/rm;
                
        // write kinetic energy
        KE += mesh.wJ[mesh.Np*e+n]*(u[n]*u[n]+v[n]*v[n]+w[n]*w[n])*rm;

        // write enstrophy
        Eps += mesh.wJ[mesh.Np*e+n]*(V[e*mesh.Np*3+n+mesh.Np*0]*V[e*mesh.Np*3+n+mesh.Np*0]
              +V[e*mesh.Np*3+n+mesh.Np*1]*V[e*mesh.Np*3+n+mesh.Np*1]
              +V[e*mesh.Np*3+n+mesh.Np*2]*V[e*mesh.Np*3+n+mesh.Np*2])*rm; 

      } else {

                // write kinetic energy
        KE += mesh.wJ[mesh.Np*e+n]*(u[n]*u[n]);
                // write enstrophy
        Eps += mesh.wJ[mesh.Np*e+n]*(V[e*mesh.Np+n]*V[e*mesh.Np+n]); 

      }
      }
    }

<<<<<<< HEAD
comm.Allreduce(KE, Comm::Sum);
comm.Allreduce(Eps, Comm::Sum);


const dfloat scale = mesh.dim==2 ? 0.5/(6.28318530718*6.28318530718): 0.5/(6.28318530718*6.28318530718*6.28318530718); 


KE *= scale; Eps *= scale;  

if(mesh.rank==0){
  fprintf(fp,"%5.2f %lf %lf\n",time, KE, Eps);
  fclose(fp);
}
  //std::cout << KE << std::endl;
  //std::cout << "Total Kinetic Energy " << KE << std::endl;
}
