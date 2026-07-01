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

  //dlong Ncoarse=0;            // Coarsened element count
  memory<dlong> EToRefLevel;  // Element Refinement List: size->(Nelements)
  memory<dlong> PToC;         // Parent to Child Connectivity: size->(Nelements*Nchild)
  memory<dlong> PCS;         // Parent, Child, Sibling Connectivity: size->(Nelements*(3))
  memory<dlong> IntFlag;      // Interpolation flag for identify different type of configurations: size->(Nelements)
  memory<dlong> RedFlag;      // 
  
  deviceMemory<dlong> o_IntFlag; 
  deviceMemory<dlong> o_EToRefLevel; 
  deviceMemory<dlong> o_PToC;
  deviceMemory<dfloat> o_IM;
  deviceMemory<dfloat> o_IMRed;
  deviceMemory<dfloat> o_RM;
  deviceMemory<dfloat> o_RMRed;

  void adaptivity(deviceMemory<dfloat>& o_q,dlong* _N);
  void Setup(platform_t& _platform, 
             mesh_t& _mesh,
             adaptivitySettings_t& _settings);
 
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
  void ConformByBisectGB(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& RefFlag2, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                memory<dlong>& new_v_id, 
                                                                dlong& Nrefine,
                                                                memory<dfloat>& EX_new, 
                                                                memory<dfloat>& EY_new,
                                                                memory<hlong>& EToV_new,
                                                                memory<int>& EToB_new,
                                                                memory<dlong>& SplitFlag,
                                                                memory<dlong>& RedFlag,                                                                             
                                                                hlong* nn,
                                                                hlong* new_vertex,
                                                                dlong RefLevel);
                                           
  void ConformByBisectRToR(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
                                                                memory<dlong>& new_v_id, 
                                                                dlong& Nrefine,
                                                                memory<dfloat>& EX_new, 
                                                                memory<dfloat>& EY_new,
                                                                memory<hlong>& EToV_new,
                                                                memory<int>& EToB_new,
                                                                memory<dlong>& SplitFlag,
                                                                memory<dlong>& RedFlag,                                                                             
                                                                hlong* nn,
                                                                hlong* new_vertex,
                                                                dlong RefLevel);                                                              
void ConformByLEMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
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
void ConformByNVBMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
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
 void ConformByRedMultiLvl(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag, memory<dlong>& FaceFlag,
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
  
  void RefinebyBisect(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag,memory<dlong>& FaceFlag, dlong Nrefine,dlong level);
  void RefineRGB(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,
                                                             memory<dlong>& RefFlag, 
                                                             memory<dlong>& ConfFlag,
                                                             memory<dlong>& confGreen,
                                                             memory<dlong>& FaceFlag, 
                                                             memory<hlong>& EToNewV,
                                                             dlong Nrefine,
                                                             dlong level);
  void RefineGreentoRed(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,
                                                             memory<dlong>& RefFlag, 
                                                             memory<dlong>& ConfFlag,
                                                             memory<dlong>& FaceFlag, 
                                                             memory<hlong>& EToNewV,
                                                             dlong Nrefine,
                                                             dlong level);
  void RefineRGB2(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,
                                                              memory<dlong>& RefFlag,
                                                              memory<dlong>& RefFlag2, 
                                                              memory<dlong>& ConfFlag,
                                                              memory<dlong>& FaceFlag, 
                                                              dlong Nrefine, 
                                                              dlong Nelements_old,
                                                              dlong level);
  void RefineRGB3(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,
                                                            memory<dlong>& RefFlag, 
                                                            memory<dlong>& ConfFlag,
                                                            memory<dlong>& FaceFlag, 
                                                            dlong Nrefine, 
                                                            dlong Nelements_old,
                                                            dlong level);  
  void RefinebyID2(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,
                                                               memory<dlong>& RefFlag, 
                                                               memory<dlong>& ConfFlag,
                                                               memory<dlong>& FaceFlag,
                                                               dlong Nrefine, 
                                                               dlong Nelements_old,
                                                               dlong level);
  void RefinebyLE(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, 
                                                                                   memory<dlong>& ConfFlag,
                                                                                   memory<dlong>& FaceFlag, 
                                                                                   dlong Nrefine, 
                                                                                   dlong Nelements_old,
                                                                                   dlong level);
  void RefinebyNV(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& ConfFlag,memory<dlong>& FaceFlag, dlong Nrefine, dlong Nelements_old, dlong level);
  void Bisect(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void BisectNew(memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,
                                                                              memory<hlong>& new_v_id,
                                                                              dlong Nrefine,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel,
                                                                              dlong ConfLevel);
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
  void Red(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,     
                                                                              memory<dlong>& RedFlag,  
                                                                              memory<hlong>& new_v_id, 
                                                                              dlong Nrefine,                                                                          
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);

    void RedOld(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,     
                                                                              memory<dlong>& RedFlag,  
                                                                              memory<hlong>& new_v_id,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
void RedLocal(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, 
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,     
                                                                              memory<dlong>& RedFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
  void Blue(memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dfloat>& EX_new, memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn);
   void Green0(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
   void Green1(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);
   void Green2(dlong e,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& FaceFlag, memory<dlong>& ConfFlag, memory<dlong>& new_v_id,
                                                                              memory<dfloat>& EX_new, 
                                                                              memory<dfloat>& EY_new,
                                                                              memory<hlong>& EToV_new,
                                                                              memory<int>& EToB_new,
                                                                              memory<dlong>& SplitFlag,                                                                           
                                                                              hlong* nn,
                                                                              hlong* new_vertex,
                                                                              dlong RefLevel);                                                                                                                                                                                                                                    

  void CoarsebyID(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, dlong Ncoarse,dlong level);
  void Coarse(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, dlong Ncoarse,dlong level);
  void CoarseRed(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag, memory<dlong>& RedFlag,  memory<dlong>& ConfFlag,dlong Ncoarse,dlong level);
  
  void CoarseGreentoRed(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& RedFlag,  memory<dlong>& ConfFlag, memory<dlong>& ConfGreen,
    memory<hlong>& EToNewV,dlong Ncoarse,dlong level);
  void CoarseBluetoRed(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& RedFlag,  memory<dlong>& ConfFlag, dlong Ncoarse,dlong level);
 
  void CoarseGreen(deviceMemory<dfloat>& o_q,memory<dfloat>& Q,memory<dfloat>& Qold,memory<dlong>& RefFlag,memory<dlong>& RedFlag,  memory<dlong>& ConfFlag, memory<hlong>& EToNewV, dlong Ncoarse,dlong level);

  // Decision for split edge
  void LongestEdge(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, dlong level);
  void LongestEdge2(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, dlong level, dlong Nrefine,memory<hlong>& new_v_id,hlong* new_vertex);
  void NewestVertex2(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, dlong level, dlong Nrefine,memory<hlong>& new_v_id,hlong* new_vertex);
  void RGB_flag(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, dlong level, dlong Nrefine,memory<hlong>& new_v_id,memory<dlong>& ConfRed,memory<dlong>& ConfGreen,memory<hlong>& EtoNewV,hlong* new_vertex); 
  void LongestEdgeConform(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, memory<dlong>& ConfFlag, dlong level, dlong Nrefine,memory<hlong>& new_v_id,hlong* new_vertex);
  void NewestVertexConform(memory<dlong>& FaceFlag, memory<dlong>& RefFlag, memory<dlong>& ConfFlag, dlong level, dlong Nrefine,memory<hlong>& new_v_id,hlong* new_vertex);
  
  void NewestVertex(memory<dlong>& FaceFlag, memory<dlong>& RefFlag);
  // Interpolation for AMR Setup
  void InterpolateToChildTri2D();
  void InterpolateToChildTri2DLE();

  void InterpolateToParentTri2D();
  void InterpolateToParentTri2DLE();

  void InterpolateToParentTri2DRed();
  void InterpolateToChildTri2DRed();

  //  Kernels
  kernel_t indicatorKernel;
  kernel_t combineKernel;
  kernel_t splitKernel;


};


}

#endif
