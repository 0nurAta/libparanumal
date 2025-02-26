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

/*void mesh_t::InterpolateToParentTri2D(){

  // Purpose: Constructing 6 different IM to interpolate solution from child to parent,
  
  Np = (N+1)*(N+2)/2;
  const memory<dfloat>r_child(Np);
  const memory<dfloat>s_child(Np);

  memory<dfloat>I1(Np*Np);
  memory<dfloat>I2(Np*Np);
  memory<dfloat>I3(Np*Np);
  memory<dfloat>I4(Np*Np);
  memory<dfloat>I5(Np*Np);
  memory<dfloat>I6(Np*Np); 

  memory<dfloat>I1T(Np*Np);
  memory<dfloat>I2T(Np*Np);
  memory<dfloat>I3T(Np*Np);
  memory<dfloat>I4T(Np*Np);
  memory<dfloat>I5T(Np*Np);
  memory<dfloat>I6T(Np*Np); 

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
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I1);
  
  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]+1);
    s_child[i] = s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-1);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]+1);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I4);

  // Fifth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(s[i]+1)+r[i];
    s_child[i] = 0.5*(s[i]-1);
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I5);

  // Sixth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]-1.0);
    s_child[i] = 0.5*(r[i]+1.0)+s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,I6);

  linAlg_t::matrixTranspose(Np, Np, I1, Np, I1T, Np);
  linAlg_t::matrixTranspose(Np, Np, I2, Np, I2T, Np);
  linAlg_t::matrixTranspose(Np, Np, I3, Np, I3T, Np);
  linAlg_t::matrixTranspose(Np, Np, I4, Np, I4T, Np);
  linAlg_t::matrixTranspose(Np, Np, I5, Np, I5T, Np);
  linAlg_t::matrixTranspose(Np, Np, I6, Np, I6T, Np);

  for (int i = 0; i < Np; ++i) {
  for (int j = 0; j < Np; ++j) {
    dfloat sum1 = 0.0;
    dfloat sum2 = 0.0;
    dfloat sum3 = 0.0;
    dfloat sum4 = 0.0;
    dfloat sum5 = 0.0;
    dfloat sum6 = 0.0;
    for (int k = 0; k < Np; ++k) {
      sum1 += invMM[i*Np + k] * I1T[k*Np + j];
      sum2 += invMM[i*Np + k] * I2T[k*Np + j];
      sum3 += invMM[i*Np + k] * I3T[k*Np + j];
      sum4 += invMM[i*Np + k] * I4T[k*Np + j];
      sum5 += invMM[i*Np + k] * I5T[k*Np + j];
      sum6 += invMM[i*Np + k] * I6T[k*Np + j];

    }
    MI1[i*Np + j] = sum1*vgeo[5*3 +  4];
    
    MI2[i*Np + j] = sum2*vgeo[5*3 +  4];
    
    MI3[i*Np + j] = sum3*vgeo[5*3 +  4];
    MI4[i*Np + j] = sum4*vgeo[5*3 +  4];
    MI5[i*Np + j] = sum5*vgeo[5*3 +  4];
    MI6[i*Np + j] = sum6*vgeo[5*3 +  4];
    printf("MI5[%d]=%f\n",i*Np + j,MI5[i*Np + j]);
    printf("MI6[%d]=%f\n",i*Np + j,MI6[i*Np + j]);

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

 
  printf("wj=%f",vgeo[5*3 +  4]);

  o_RM = platform.malloc<dfloat>(MRMT);

}*/

void mesh_t::InterpolateToParentTri2D(){

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
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R1);


  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(r[i]+1);
    s_child[i] = s[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]-1);
    r_child[i] = r[i];
  }
  InterpolationMatrixTri2D(N,r,s,r_child,s_child,R3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(s[i]+1);
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


  
}

} //namespace libp
