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

#ifndef ADAPTIVITY_HPP
#define ADAPTIVITY_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "settings.hpp"
#include "ogs.hpp"
#include "mesh.hpp"


namespace libp {

class adaptivitySettings_t: public settings_t {
public:
  adaptivitySettings_t() = default;
  adaptivitySettings_t(comm_t _comm);
  void report();
};

class adaptivity_t {
 public:

  adaptivity_t() = default;

  platform_t platform;
  adaptivitySettings_t settings;
  properties_t props;
  mesh_t mesh;

  dlong Ncoarse=0;            // Coarsened element count
  memory<dlong> EToRefLevel;  // Element Refinement List: size->(Nelements)
  memory<dlong> PToC;         // Parent to Child Connectivity: size->(Nelements*Nchild)
  memory<dlong> IntFlag;      // Interpolation flag for identify different type of configurations: size->(Nelements)
  
  deviceMemory<dlong> o_IntFlag; 
  deviceMemory<dlong> o_PToC;
  deviceMemory<dfloat> o_IM;
  deviceMemory<dfloat> o_RM;

  void adaptivity(deviceMemory<dfloat>& o_q,dlong* _N);
  void Conform(memory<dlong>& RefFlag,  memory<dlong>& FaceFlag, dlong Nrefine);
  void Refine(memory<dfloat>& Q,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong Nrefine);
  void Bisect(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void Red(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn);
  void Blue(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn);

  void Coarse(deviceMemory<dfloat>& o_q,memory<dlong>& RefFlag, dlong Ncoarse);
  void LongestEdge(memory<dlong>& FaceFlag, memory<dlong>& RefFlag);

  // Interpolation for AMR Setup
  void InterpolateToChildTri2D();
  void InterpolateToParentTri2D();

  //  Kernels
  kernel_t indicatorKernel;
  kernel_t combineKernel;
  kernel_t splitKernel;


};


}

#endif