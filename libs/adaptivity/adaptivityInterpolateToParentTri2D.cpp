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

/*void mesh_t::InterpolateToParentTri2D(){

  // Purpose: Constructing 6 different IM to interpolate solution from child to parent,
  
  Np = (N+1)*(N+2)/2;
  const memory<dfloat>r_child(Np);
  const memory<dfloat>s_child(Np);

  memory<dfloat>R1(Np*Np);
  memory<dfloat>R2(Np*Np);
  memory<dfloat>R3(Np*Np);
  memory<dfloat>R4(Np*Np);
  memory<dfloat>R5(Np*Np);
  memory<dfloat>R6(Np*Np); 

  memory<dfloat>R1T(Np*Np);
  memory<dfloat>R2T(Np*Np);
  memory<dfloat>R3T(Np*Np);
  memory<dfloat>R4T(Np*Np);
  memory<dfloat>R5T(Np*Np);
  memory<dfloat>R6T(Np*Np); 

  memory<dfloat>MI1(Np*Np);
  memory<dfloat>MI2(Np*Np);
  memory<dfloat>MI3(Np*Np);
  memory<dfloat>MI4(Np*Np);
  memory<dfloat>MI5(Np*Np);
  memory<dfloat>MI6(Np*Np); 

  for (int i = 0; i < Np*Np; ++i)
  {
    MI1[i]= 0.0;
    MI2[i]= 0.0;
    MI3[i]= 0.0;
    MI4[i]= 0.0;
    MI5[i]= 0.0;
    MI6[i]= 0.0;
  }

  memory<dfloat>MRMT(6*Np*Np);

  // First: r_child = (r-1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-1);
    s_child[i] = s[i];
    printf("r[%d]=%f,s[%d]=%f\n", r[i],s[i]);
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R1);

  //printf("r[%d]=%f,s[%d]=%f\n", r[i],s[i]);
  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-s[i]);
    s_child[i] = s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-1.0);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-r[i]);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R4);

  // Fifth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(s[i]+1)+r[i];
    s_child[i] = 0.5*(s[i]-1);
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R5);

  // Sixth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-1.0);
    s_child[i] = 0.5*(r[i]+1.0)+s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R6);

  linAlg_t::matrixTranspose(Np, Np, R1, Np, R1T, Np);
  linAlg_t::matrixTranspose(Np, Np, R2, Np, R2T, Np);
  linAlg_t::matrixTranspose(Np, Np, R3, Np, R3T, Np);
  linAlg_t::matrixTranspose(Np, Np, R4, Np, R4T, Np);
  linAlg_t::matrixTranspose(Np, Np, R5, Np, R5T, Np);
  linAlg_t::matrixTranspose(Np, Np, R6, Np, R6T, Np);

  for (int i = 0; i < Np; ++i) {
  for (int j = 0; j < Np; ++j) {
    dfloat sum1 = 0.0;
    dfloat sum2 = 0.0;
    dfloat sum3 = 0.0;
    dfloat sum4 = 0.0;
    dfloat sum5 = 0.0;
    dfloat sum6 = 0.0;
    for (int k = 0; k < Np; ++k) {
      sum1 += invMM[i*Np + k] * R1T[k*Np + j];
      sum2 += invMM[i*Np + k] * R2T[k*Np + j];
      sum3 += invMM[i*Np + k] * R3T[k*Np + j];
      sum4 += invMM[i*Np + k] * R4T[k*Np + j];
      sum5 += invMM[i*Np + k] * R5T[k*Np + j];
      sum6 += invMM[i*Np + k] * R6T[k*Np + j];

    }
    MI1[i*Np + j] = sum1;
    MI2[i*Np + j] = sum2;
    MI3[i*Np + j] = sum3;
    MI4[i*Np + j] = sum4;
    MI5[i*Np + j] = sum5;
    MI6[i*Np + j] = sum6;
    //printf("MI5[%d]=%f\n",i*Np + j,MI5[i*Np + j]);
    //printf("MI6[%d]=%f\n",i*Np + j,MI6[i*Np + j]);

  }
}

  memory<dfloat>MI1T = MRMT + 0*Np*Np;
  memory<dfloat>MI2T = MRMT + 1*Np*Np;
  memory<dfloat>MI3T = MRMT + 2*Np*Np;
  memory<dfloat>MI4T = MRMT + 3*Np*Np;
  memory<dfloat>MI5T = MRMT + 4*Np*Np;
  memory<dfloat>MI6T = MRMT + 5*Np*Np;

  linAlg_t::matrixTranspose(Np, Np, MI1, Np, MI1T, Np);
  linAlg_t::matrixTranspose(Np, Np, MI2, Np, MI2T, Np);
  linAlg_t::matrixTranspose(Np, Np, MI3, Np, MI3T, Np);
  linAlg_t::matrixTranspose(Np, Np, MI4, Np, MI4T, Np);
  linAlg_t::matrixTranspose(Np, Np, MI5, Np, MI5T, Np);
  linAlg_t::matrixTranspose(Np, Np, MI6, Np, MI6T, Np);

 
//  printf("wj=%f",vgeo[5*3 +  4]);

  o_RM = platform.malloc<dfloat>(MRMT);

}*/

/*void mesh_t::InterpolateToParentTri2D(){

  // Purpose: Constructing 6 different IM to interpolate solution from child to parent,
  
  Np = (N+1)*(N+2)/2;
  const memory<dfloat>r_child(Np);
  const memory<dfloat>s_child(Np);

  memory<dfloat>R1(Np*Np);
  memory<dfloat>R2(Np*Np);
  memory<dfloat>R3(Np*Np);
  memory<dfloat>R4(Np*Np);
  memory<dfloat>R5(Np*Np);
  memory<dfloat>R6(Np*Np); 

  memory<dfloat>MRMT(6*Np*Np);

  // First: r_child = (r-1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-1);
    s_child[i] = s[i];
    printf("r[%d]=%f,s[%d]=%f\n", r[i],s[i]);
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R1);

  //printf("r[%d]=%f,s[%d]=%f\n", r[i],s[i]);
  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-s[i]);
    s_child[i] = s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-1.0);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-r[i]);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R4);

  // Fifth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(s[i]+1)+r[i];
    s_child[i] = 0.5*(s[i]-1);
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R5);

  // Sixth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-1.0);
    s_child[i] = 0.5*(r[i]+1.0)+s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R6);

  linAlg_t::matrixInverse(Np, R1);
  linAlg_t::matrixInverse(Np, R2);
  linAlg_t::matrixInverse(Np, R3);
  linAlg_t::matrixInverse(Np, R4);
  linAlg_t::matrixInverse(Np, R5);
  linAlg_t::matrixInverse(Np, R6);

  memory<dfloat>R1T = MRMT + 0*Np*Np;
  memory<dfloat>R2T = MRMT + 1*Np*Np;
  memory<dfloat>R3T = MRMT + 2*Np*Np;
  memory<dfloat>R4T = MRMT + 3*Np*Np;
  memory<dfloat>R5T = MRMT + 4*Np*Np;
  memory<dfloat>R6T = MRMT + 5*Np*Np;

  linAlg_t::matrixTranspose(Np, Np, R1, Np, R1T, Np);
  linAlg_t::matrixTranspose(Np, Np, R2, Np, R2T, Np);
  linAlg_t::matrixTranspose(Np, Np, R3, Np, R3T, Np);
  linAlg_t::matrixTranspose(Np, Np, R4, Np, R4T, Np);
  linAlg_t::matrixTranspose(Np, Np, R5, Np, R5T, Np);
  linAlg_t::matrixTranspose(Np, Np, R6, Np, R6T, Np);

  o_RM = platform.malloc<dfloat>(MRMT);


  
}*/

void adaptivity_t::InterpolateToParentTri2D(){

  // Purpose: Constructing 3 different RM to interpolate solution from child to parent
  // where r_p = (IM^T*IM)^-1*IM^T*r_c
  
  const dfloat Np = (mesh.N+1)*(mesh.N+2)/2; // Number of points in a triangular element

  const memory<dfloat>r_child(Np);
  const memory<dfloat>s_child(Np);

  memory<dfloat>R12(2*Np*Np);
  memory<dfloat>R34(2*Np*Np);
  memory<dfloat>R56(2*Np*Np);
  
  //memory<dfloat>R1 = R12 + 0*Np*Np;
  //memory<dfloat>R2 = R12 + 1*Np*Np;
  //memory<dfloat>R3 = R34 + 0*Np*Np;
  //memory<dfloat>R4 = R34 + 1*Np*Np;
  //memory<dfloat>R5 = R56 + 0*Np*Np;
  //memory<dfloat>R6 = R56 + 1*Np*Np;

  memory<dfloat>R1;
  memory<dfloat>R2;
  memory<dfloat>R3;
  memory<dfloat>R4;
  memory<dfloat>R5;
  memory<dfloat>R6;

  memory<dfloat>R12T(2*Np*Np);
  memory<dfloat>R34T(2*Np*Np);
  memory<dfloat>R56T(2*Np*Np);

  memory<dfloat>B1(Np*Np);
  memory<dfloat>B2(Np*Np);
  memory<dfloat>B3(Np*Np);


  memory<dfloat>MRMT(6*Np*Np);

  // Calculate IM for 3 different bisection config.

 /*   for (int i = 0; i < Np; ++i)
  {
    printf("r[%d]=%f\n",i,r[i]);
    
  }
    for (int i = 0; i < Np; ++i)
  {
    printf("s[%d]=%f\n",i,s[i]);
  }
*/
    
  // First: r_child = (r-1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-1);
    s_child[i] = mesh.s[i];
    
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R1);


  //printf("r[%d]=%f,s[%d]=%f\n", r[i],s[i]);
  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-mesh.s[i]);
    s_child[i] = mesh.s[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(mesh.s[i]-1.0);
    r_child[i] = mesh.r[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(mesh.s[i]-mesh.r[i]);
    r_child[i] = mesh.r[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R4);

  // Fifth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.s[i]+1)+mesh.r[i];
    s_child[i] = 0.5*(mesh.s[i]-1);
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R5);

  // Sixth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-1.0);
    s_child[i] = 0.5*(mesh.r[i]+1.0)+mesh.s[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,R6);


  /*for (int i = 0; i < Np*Np; ++i)
  {
    printf("R1[%d]=%f\n",i, R1[i]);
  }
    for (int i = 0; i < Np*Np; ++i)
  {
    printf("R2[%d]=%f\n",i, R2[i]);
  }*/

    for (int i = 0; i < Np*Np; ++i)
  {
    R12[i] = R1[i];
    R12[i+Np*Np] = R2[i];
    R34[i] = R3[i];
    R34[i+Np*Np] = R4[i];
    R56[i] = R5[i];
    R56[i+Np*Np] = R6[i];
    //printf("R12[%d]=%f\n",i, R12[i]);
  }

  //  for (int i = 0; i < 2*Np*Np; ++i)
  //{
  //  printf("R12[%d]=%f\n",i, R12[i]);
  //}

  linAlg_t::matrixTranspose(2*Np, Np, R12, Np, R12T, 2*Np);
  linAlg_t::matrixTranspose(2*Np, Np, R34, Np, R34T, 2*Np);
  linAlg_t::matrixTranspose(2*Np, Np, R56, Np, R56T, 2*Np);


// Calculating B=(IM^T*IM)

    for (int i = 0; i < Np; ++i) {
  for (int j = 0; j < Np; ++j) {
    dfloat sum1 = 0.0;
    dfloat sum2 = 0.0;
    dfloat sum3 = 0.0;

    for (int k = 0; k < 2*Np; ++k) {
      sum1 += R12T[i*2*Np + k] * R12[k*Np + j];
      sum2 += R34T[i*2*Np + k] * R34[k*Np + j];
      sum3 += R56T[i*2*Np + k] * R56[k*Np + j];

    }
    B1[i*Np + j] = sum1;
    //printf("B1[%d]=%f\n",i*Np + j, B1[i*Np + j]);
    B2[i*Np + j] = sum2;
    B3[i*Np + j] = sum3;
  
  }
}
  
  linAlg_t::matrixInverse(Np, B1);
  linAlg_t::matrixInverse(Np, B2);
  linAlg_t::matrixInverse(Np, B3);
  /*for (int i = 0; i < Np*Np; ++i)
  {
    printf("B1_inv[%d]=%f\n",i, B1[i]);
  }*/
  
  //memory<dfloat>RM12(2*Np*Np);
  //memory<dfloat>RM34(2*Np*Np);
  //memory<dfloat>RM56(2*Np*Np);

  memory<dfloat>RM12 = MRMT + 0*Np*Np;
  memory<dfloat>RM34 = MRMT + 2*Np*Np;
  memory<dfloat>RM56 = MRMT + 4*Np*Np;

 

    for (int i = 0; i < Np; ++i) {
  for (int j = 0; j < 2*Np; ++j) {
    dfloat sum1 = 0.0;
    dfloat sum2 = 0.0;
    dfloat sum3 = 0.0;

    for (int k = 0; k < Np; ++k) {
      sum1 += B1[i*Np + k] * R12T[2*k*Np + j];
      sum2 += B2[i*Np + k] * R34T[2*k*Np + j];
      sum3 += B3[i*Np + k] * R56T[2*k*Np + j];

    }
    RM12[2*i*Np + j] = sum1;
    RM34[2*i*Np + j] = sum2;
    RM56[2*i*Np + j] = sum3;
    //printf("RM56[%d]=%f\n",2*i*Np + j, RM56[2*i*Np + j]);
  
  }
}

 // memory<dfloat>RMT12 = MRMT + 0*Np*Np;
 // memory<dfloat>RMT34 = MRMT + 2*Np*Np;
 // memory<dfloat>RMT56 = MRMT + 4*Np*Np;
//  linAlg_t::matrixRightSolve( Np, 2*Np, R12T, Np,Np, B1,  RM12);
//  linAlg_t::matrixRightSolve( Np, 2*Np, R34T, Np,Np, B2,  RM34);
//  linAlg_t::matrixRightSolve( Np, 2*Np, R56T, Np,Np, B3,  RM56);




  //linAlg_t::matrixTranspose(Np, 2*Np, RM12, 2*Np, RMT12, Np);
  //linAlg_t::matrixTranspose(Np, 2*Np, RM34, 2*Np, RMT34, Np);
  //linAlg_t::matrixTranspose(Np, 2*Np, RM56, 2*Np, RMT56, Np);


  o_RM = platform.malloc<dfloat>(MRMT);


  
}

} //namespace libp
