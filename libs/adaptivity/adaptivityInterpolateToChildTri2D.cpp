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

void adaptivity_t::InterpolateToChildTri2D(){

  // Purpose: Constructing 6 different IM to interpolate solution from parent to child,
  
  const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
  const memory<dfloat>r_child(Np);
  const memory<dfloat>s_child(Np);

  memory<dfloat>I1(Np*Np);
  memory<dfloat>I2(Np*Np);
  memory<dfloat>I3(Np*Np);
  memory<dfloat>I4(Np*Np);
  memory<dfloat>I5(Np*Np);
  memory<dfloat>I6(Np*Np); 
  memory<dfloat> IMT(6*Np*Np);

  // First: r_child = (r-1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-1);
    s_child[i] = mesh.s[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I1);
  
  // Second: r_child = (r+1)/2
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-mesh.s[i]);
    s_child[i] = mesh.s[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I2);
  
  // Third: s_child = (s-1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(mesh.s[i]-1);
    r_child[i] = mesh.r[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I3);
  
  // Fourth: s_child = (s+1)/2
  for (int i = 0; i < Np; ++i)
  {
    s_child[i] = 0.5*(mesh.s[i]-mesh.r[i]);
    r_child[i] = mesh.r[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I4);

  // Fifth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.s[i]+1)+mesh.r[i];
    s_child[i] = 0.5*(mesh.s[i]-1);
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I5);

  // Sixth: 
  for (int i = 0; i < Np; ++i)
  {
    r_child[i] = 0.5*(mesh.r[i]-1.0);
    s_child[i] = 0.5*(mesh.r[i]+1.0)+mesh.s[i];
  }
  mesh.InterpolationMatrixTri2D(mesh.N,mesh.r,mesh.s,r_child,s_child,I6);


  memory<dfloat>I1T = IMT + 0*Np*Np;
  memory<dfloat>I2T = IMT + 1*Np*Np;
  memory<dfloat>I3T = IMT + 2*Np*Np;
  memory<dfloat>I4T = IMT + 3*Np*Np;
  memory<dfloat>I5T = IMT + 4*Np*Np;
  memory<dfloat>I6T = IMT + 5*Np*Np;

  linAlg_t::matrixTranspose(Np, Np, I1, Np, I1T, Np);
  linAlg_t::matrixTranspose(Np, Np, I2, Np, I2T, Np);
  linAlg_t::matrixTranspose(Np, Np, I3, Np, I3T, Np);
  linAlg_t::matrixTranspose(Np, Np, I4, Np, I4T, Np);
  linAlg_t::matrixTranspose(Np, Np, I5, Np, I5T, Np);
  linAlg_t::matrixTranspose(Np, Np, I6, Np, I6T, Np);
  
  o_IM = platform.malloc<dfloat>(IMT);


}

} //namespace libp
