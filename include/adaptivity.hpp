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

  // constructor
  adaptivity_t(platform_t& _platform, 
               mesh_t& _mesh,
               adaptivitySettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  platform_t platform;
  adaptivitySettings_t settings;
  properties_t props;
  mesh_t mesh;

  dlong Ncoarse=0;            // Coarsened element count
  memory<dlong> EToRefLevel;  // Element Refinement List: size->(Nelements)
  memory<dlong> PToC;         // Parent to Child Connectivity: size->(Nelements*Nchild)
  memory<dlong> PCS;         // Parent, Child, Sibling Connectivity: size->(Nelements*(3))
  memory<dlong> IntFlag;      // Interpolation flag for identify different type of configurations: size->(Nelements)
  
  deviceMemory<dlong> o_IntFlag; 
  deviceMemory<dlong> o_EToRefLevel; 
  deviceMemory<dlong> o_PToC;
  deviceMemory<dfloat> o_IM;
  deviceMemory<dfloat> o_RM;

  void adaptivity(deviceMemory<dfloat>& o_q,dlong* _N);
  void Setup(platform_t& _platform, 
             mesh_t& _mesh,
             adaptivitySettings_t& _settings);
  void Conform(memory<dlong>& RefFlag,  memory<dlong>& FaceFlag,dlong& Nrefine);
  void ConformLE(memory<dlong>& RefFlag,  memory<dlong>& FaceFlag, dlong& Nrefine);
  void ConformByVertex(memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag, dlong& Nrefine);
  void ConformByID(memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,memory<dlong>& new_v_id, dlong& Nrefine);
  void ConformByBisect(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                        memory<dlong>& new_v_id, 
                                                                        dlong& Nrefine,
                                                                        memory<dfloat>& EX_new, 
                                                                        memory<dfloat>& EY_new,
                                                                        memory<hlong>& EToV_new,
                                                                        memory<int>& EToB_new,
                                                                        memory<dlong>& SplitFlag,                                                                           
                                                                        hlong* nn,
                                                                        hlong* new_vertex,
                                                                        dlong RefLevel);
  void ConformByBisectMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                        memory<dlong>& new_v_id, 
                                                                        dlong& Nrefine,
                                                                        memory<dfloat>& EX_new, 
                                                                        memory<dfloat>& EY_new,
                                                                        memory<hlong>& EToV_new,
                                                                        memory<int>& EToB_new,
                                                                        memory<dlong>& SplitFlag,                                                                           
                                                                        hlong* nn,
                                                                        hlong* new_vertex,
                                                                        dlong RefLevel);
  void ConformByCoarse(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                        memory<dlong>& new_v_id, 
                                                                        dlong& Nrefine,
                                                                        memory<dfloat>& EX_new, 
                                                                        memory<dfloat>& EY_new,
                                                                        memory<hlong>& EToV_new,
                                                                        memory<int>& EToB_new,
                                                                        memory<dlong>& SplitFlag,                                                                           
                                                                        hlong* nn,
                                                                        hlong* new_vertex,
                                                                        dlong RefLevel);
  
  void Refine(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong Nrefine);
  void RefineLE(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, dlong Nrefine);
  void RefinebyID(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag,memory<dlong>& FaceFlag, dlong Nrefine);
  void RefinebyID2(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag,memory<dlong>& FaceFlag, dlong Nrefine, dlong Nelements_old);
  
  void Bisect(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectNEW(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectbyID(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectbyID2(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectLocal0(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectLocal1(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectLocal2(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
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

  void Coarse(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dlong>& RefFlag, dlong Ncoarse);
  void CoarsebyID(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, dlong Ncoarse,dlong level);
  void LongestEdge(memory<dlong>& FaceFlag, memory<dlong>& RefFlag);
  void LongestEdgeNEW(memory<dlong>& FaceFlag, memory<dlong>& RefFlag);
  // Interpolation for AMR Setup
  void InterpolateToChildTri2D();
  void InterpolateToChildTri2DLE();
  void InterpolateToParentTri2D();
  void InterpolateToParentTri2DLE();

  //  Kernels
  kernel_t indicatorKernel;
  kernel_t combineKernel;
  kernel_t splitKernel;


};


}

#endif