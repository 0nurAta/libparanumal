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

void mesh_t::CurvedNodesQuad2D(){

  //x.malloc(Nelements*Np);
  //y.malloc(Nelements*Np);

  /* Origin coordinates and radius of the cylinder */
  dfloat x0 = 0.0;
  dfloat y0 = 0.0;
  dfloat r_cyl =1.0;

  //LIBP_FORCE_ABORT("Stopped");


  //#pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToB[e*Nfaces+f]==1) {
      
       

        dlong id = e*Nverts;

        dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
        dfloat xe4 = EX[id+3];
        dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
        dfloat ye4 = EY[id+3];

        dfloat theta1 = atan2(ye1-y0,xe1-x0); 
        dfloat theta2 = atan2(ye4-y0,xe4-x0);

        //check to make sure they are in the same quadrant
        if ((theta2 > 0) && (theta1 < 0)){
          printf("theta1 =%f,theta2=%f",theta1,theta2);
          theta1 = theta1;
          printf("theta1 =%f,theta2=%f",theta1,theta2);} 
        if ((theta1 > 0) && (theta2 < 0)){ 
          printf("theta1 =%f,theta2=%f",theta1,theta2);
          theta2 = theta2 + 2*M_PI;
          printf("theta1 =%f,theta2=%f",theta1,theta2);}

        dfloat theta[N+1];
        dfloat fdx[N+1];
        dfloat fdy[N+1];

        for (int n=0;n<(N+1);++n)
        {
          theta[n]=0.5*theta1*(1-r[n]) + 0.5*theta2*(1+r[n]);
          fdx[n] = x0 + r_cyl*cos(theta[n])-x[e*Np + (N+1)*n]; 
          fdy[n] = y0 + r_cyl*sin(theta[n])-y[e*Np + (N+1)*n];
          //printf("r =%f",r[n]);
        }

      

       memory<dfloat> V_vol;
       mesh_t::Vandermonde1D(N, s, V_vol);
       
       memory<dfloat> rf(N+1);
       //dfloat r_face[N];
       


        for (int i = 0; i < N+1; ++i)
        {
         rf[i] = r[i];
         //printf("%f \n",s[i]);
         //printf("%f \n",r[i]);
        }
        /*for (int i = 0; i < (16); ++i)
        {
         printf("s= %f\n",r[i]);
        }*/
        

        memory<dfloat> V_face;
        mesh_t::Vandermonde1D(N, rf, V_face);
        




       linAlg_t::matrixInverse((N+1), V_face);

        /*for (int i = 0; i < (64); ++i)
        {
         printf("x2= %f\n",V_vol[i]);
        }
        */
       
 
        //printf("HEREEEE\n");
        //LIBP_FORCE_ABORT("Stopped");


        //const int Npoints = rf.length();
        dfloat vdx[N+1]= { };
        dfloat vdy[N+1]= { };
        dfloat Vdx[Np] = { };
        dfloat Vdy[Np] = { };
      
        const int Npoints = V_face.length();
        //printf("Npoints=%d \n",Np);

       /*for (int i = 0; i < (Np); ++i)
        {
         printf("x2= %f\n",Vdx[i]);
        }*/


        for (int n=0;n<(N+1);++n)
        {
        for(int m=0;m<(N+1);++m){
      
        vdx[n] += V_face[n*(N+1)+m]*fdx[m];
        vdy[n] += V_face[n*(N+1)+m]*fdy[m];
        //printf("x2= %f\n",V_face[n*(N+1)+m]);
        //printf("x2= %d\n",n*(N+1)+m);
        }
        }


 
        for (int n=0;n<Np;++n)
        {
        for(int m=0;m<(N+1);++m){
      
        Vdx[n] += V_vol[n*(N+1)+m]*vdx[m];
        Vdy[n] += V_vol[n*(N+1)+m]*vdy[m];
        //printf("x2= %f\n",V_vol[n*(N+1)+m]*vdx[m] );
        //printf("Vdx= %f\n",Vdx[n]);
        }

        }

       //printf("Npoints=%d \n",Np);

     /*    for (int i = 0; i < (16); ++i)
        {
         printf("V_face= %f\n",V_face[i]);
        }

        for (int i = 0; i < (4); ++i)
        {
         printf("vdx= %f\n",vdx[i]);
        }
*/
        /*for (int i = 0; i < (16); ++i)
        {
         printf("Vvol= %f\n",Vdx[i]);
        }*/

        if (e==196)
        {
                  for (int i = 0; i < (16); ++i)
        {
         //printf("fx=%f \n",fdx[i]);
         //printf("x=%f \n",x[e*Np +i]);
         //printf("fx=%f \n",y[e*Np + 4*i]);
         
        }
          //LIBP_FORCE_ABORT("Stopped");
        }
        dfloat     xold;
        dfloat     yold;
        for (int n=0;n<Np;++n)
        {
          //printf("x1= %f\n",y[e*Np + n] );
          //printf("x1= %f\n",x[e*Np + n] );
         /* if (e==196)
        {
             xold = x[e*Np + n];
             yold = y[e*Np + n];
        }*/
          x[e*Np + n] = x[e*Np + n]+0.5*(1.0-r[n])*Vdx[n]; 
          y[e*Np + n] = y[e*Np + n]+0.5*(1.0-r[n])*Vdy[n];

           /*         if (e==196)
        {
             x[e*Np + n]=xold;
             y[e*Np + n]=yold;
        }*/
          //printf("x2= %f\n",0.5*(1.0-r[n])*Vdx[n]);
        }  


        
         //LIBP_FORCE_ABORT("Stopped");
      }
    }
  }

}


} //namespace libp
