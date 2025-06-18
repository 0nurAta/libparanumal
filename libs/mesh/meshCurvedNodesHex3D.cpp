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

#include "mesh.hpp"

namespace libp {

void mesh_t::CurvedNodesHex3D(){

  //x.malloc(Nelements*Np);
  //y.malloc(Nelements*Np);

  /* Origin coordinates and radius of the sphere */
  dfloat x0 = 0.0;
  dfloat y0 = 0.0;
  dfloat z0 = 0.0;
  dfloat r_cyl =0.5;

  //LIBP_FORCE_ABORT("Stopped");


  //#pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToB[e*Nfaces+f]==1) {        

        // Store coordinates of vertices 
        dfloat x0 = x[e*Np];         dfloat y0 = y[e*Np];         dfloat z0 = z[e*Np];
        dfloat x1 = x[e*Np+N];       dfloat y1 = y[e*Np+N];       dfloat z1 = z[e*Np+N];
        dfloat x2 = x[e*Np+(N+1)*N]; dfloat y2 = y[e*Np+(N+1)*N]; dfloat z2 = z[e*Np+(N+1)*N];  
        dfloat x3 = x[e*Np+(N+2)*N]; dfloat y3 = x[e*Np+(N+2)*N]; dfloat z3 = x[e*Np+(N+2)*N];
        
        //
        memory<dfloat> A(9);
        A[0] = x0-x3; A[1] = y0-y3; A[2] = z0-z3;
        A[3] = x1-x3; A[4] = y1-y3; A[5] = z1-z3;
        A[6] = x2-x3; A[7] = y2-y3; A[8] = z2-z3;

      
        //
        dfloat b[3];
        dfloat mag0 = x0*x0+y0*y0+z0*z0;
        dfloat mag1 = x1*x1+y1*y1+z1*z1;
        dfloat mag2 = x2*x2+y2*y2+z2*z2;
        dfloat mag3 = x3*x3+y3*y3+z3*z3;

        b[0] = 0.5*(mag0-mag3);
        b[1] = 0.5*(mag1-mag3);
        b[2] = 0.5*(mag2-mag3);


        dfloat c[3] = {};
        //linAlg_t::matrixInverse((3), A);

               for (int i = 0; i < 9; ++i)
       {
          //printf("c = %f \n", A[i]); 
       }
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 9; ++j)
          {
            c[i] += A[j]*b[i];
            
          }
         // printf("c = %f \n", c[i]);
        }

        //LIBP_FORCE_ABORT("Stopped");
        
        // Deformation array
        dfloat fdx[(N+1)*(N+1)];
        dfloat fdy[(N+1)*(N+1)];
        dfloat fdz[(N+1)*(N+1)];

        // Correcting boundary face points by normalizing them with their radius and calculating deformations
        for (int i=0;i<(N+1)*(N+1);++i)
        {
        
          dfloat R = sqrt(x[e*Np + i]*x[e*Np + i]+y[e*Np + i]*y[e*Np + i]+z[e*Np + i]*z[e*Np + i])/r_cyl;
          //printf("R=%f\n",R); 
          fdx[i] = x[e*Np + i]/R - x[e*Np + i];
          fdy[i] = y[e*Np + i]/R - y[e*Np + i];
          fdz[i] = z[e*Np + i]/R - z[e*Np + i];

          //printf("fdx=%f\n",fdx[i]);
        }

        
        // Storing calculated deformation at face to be used for volume nodes
        dfloat     fDX[Np];
        dfloat     fDY[Np];
        dfloat     fDZ[Np];
        for (int i = 0; i < (N+1); ++i)
        {
          for (int m = 0; m <((N+1)*(N+1)) ; ++m)
          {
            fDX[i*(N+1)*(N+1)+m] = fdx[m];
            fDY[i*(N+1)*(N+1)+m] = fdy[m];
            fDZ[i*(N+1)*(N+1)+m] = fdz[m];
            
            
            
          }
        }
        // Blending deformation
        for (int n=0;n<Np;++n)
        {
          
          x[e*Np + n] = x[e*Np + n]+0.5*(1.0-t[n])*fDX[n]; 
          y[e*Np + n] = y[e*Np + n]+0.5*(1.0-t[n])*fDY[n];
          z[e*Np + n] = z[e*Np + n]+0.5*(1.0-t[n])*fDZ[n];
        //  if (e==37066)
        //{
          //   x[e*Np + n]=xold;
          //   y[e*Np + n]=yold;
          //   z[e*Np + n]=zold;
             //printf("x2= %f y2= %f z2= %f\n",x[e*Np + n],y[e*Np + n],z[e*Np + n] );
        //}
          //printf("x2= %f\n",0.5*(1.0-r[n])*Vdx[n]);
        }  


        
         //LIBP_FORCE_ABORT("Stopped");
      }
    }
  }

}


} //namespace libp
