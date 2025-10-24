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

#ifndef ADVECTION_HPP
#define ADVECTION_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DADVECTION LIBP_DIR"/solvers/advection/"

using namespace libp;

class advectionSettings_t: public settings_t {
public:
  advectionSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename);
};

class advection_t: public solver_t {
public:
  mesh_t mesh;
  timeStepper_t timeStepper;

  ogs::halo_t traceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;

  memory<dlong> EToRefLevel;  // Element Refinement List: size->(Nelements)
  memory<dlong> PToC;         // Parent to Child Connectivity: size->(Nelements*Nchild)
  deviceMemory<dlong> o_PToC; 
  //memory<dlong> CToP;         // Child to Parent Connectivity: size->(Nelements)
  memory<dlong> IntFlag;      // Interpolation flag for identify different type of configurations: size->(Nelements)
  deviceMemory<dlong> o_IntFlag; 
  dlong Ncoarse=0;       // Coarsened element count


  kernel_t volumeKernel;
  kernel_t surfaceKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  // AMR Kernels
  kernel_t indicatorKernel;
  kernel_t combineKernel;
  kernel_t splitKernel;

  advection_t() = default;
  advection_t(platform_t &_platform, mesh_t &_mesh,
              advectionSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& platform, mesh_t& mesh,
             advectionSettings_t& settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q,  const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  //// AMR ////
  void Amr(deviceMemory<dfloat>& o_q,dlong* _N);

  void Conform(memory<dlong>& RefFlag,  memory<dlong>& FaceFlag, dlong Nrefine);

  void Refine(memory<dfloat>& Q,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong Nrefine);

  void Bisect(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
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


  void Coarse(memory<dfloat>& Q,memory<dlong>& RefFlag, dlong Ncoarse);

  void LongestEdge(memory<dlong>& FaceFlag, memory<dlong>& RefFlag);

  //// AMR ENDS ////
  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);
};

#endif
