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
// Conduct Adaptive Mesh Refinement
void adaptivity_t::Setup(platform_t& _platform, 
                         mesh_t& _mesh,
                         adaptivitySettings_t& _settings){

  platform = _platform;
  settings = _settings;
  mesh = _mesh;
  props = platform.props();


  //comm = _comm.Dup();
  //rank = comm.rank();
  //size = comm.size();
  
  // compute interpolation matrices for amr
  //InterpolateToChildTri2D();
  InterpolateToChildTri2DLE();
  InterpolateToParentTri2DLE();
  // For RGB
  InterpolateToChildTri2DRed();
  InterpolateToParentTri2DRed();

  dlong const MAX_LEVEL = 12;
  
  // Set Lists related to AMR
  EToRefLevel.calloc(128*mesh.Nelements); //
  PToC.calloc(128*mesh.Nelements*(MAX_LEVEL+3)); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  PCS.calloc(128*mesh.Nelements*3); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  IntFlag.calloc(128*mesh.Nelements*(MAX_LEVEL+3)); //
  RedFlag.malloc(128*mesh.Nelements*(MAX_LEVEL+3),-1);
  
  for (int i = 0; i < 128*mesh.Nelements; ++i)
  {
    for (int n = 0; n < (MAX_LEVEL+3); ++n)
    {
      dlong id = i*(MAX_LEVEL+3)+n;
       PToC[id] = -1;
    }
    
  }

// OCCA build stuff

  properties_t kernelInfo = mesh.props; //copy base occa properties

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"] = NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"] = NblockS;

  kernelInfo["defines/" "p_Nfields"] = 1;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = ADAPTIVITY_DIR "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

    // indicator kernel
  fileName   = oklFilePrefix + "Indicator" + suffix + oklFileSuffix;
  kernelName = "IndicatorTest" + suffix;
  indicatorKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  // combine solution kernel
  fileName   = oklFilePrefix + "Combine" + suffix + oklFileSuffix;
  kernelName = "Combine" + suffix;
  combineKernel = platform.buildKernel(fileName, kernelName,
                                       kernelInfo);
  // split solution kernel
  fileName   = oklFilePrefix + "Split" + suffix + oklFileSuffix;
  kernelName = "Split" + suffix;
  splitKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);
  
  // candidate kernel
  fileName   = oklFilePrefix + "Candidate" + suffix + oklFileSuffix;
  kernelName = "Candidate" + suffix;
  candidateKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  // assign kernel
  fileName   = oklFilePrefix + "Assign" + suffix + oklFileSuffix;
  kernelName = "AssignSimple" + suffix;
  assignKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  // assign conform kernel
  fileName   = oklFilePrefix + "AssignConform" + suffix + oklFileSuffix;
  kernelName = "AssignConformSimple" + suffix;
  assignconformKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
                                         
  // conform kernel
  fileName   = oklFilePrefix + "Conform" + suffix + oklFileSuffix;
  kernelName = "ConformCandidate" + suffix;
  conformKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  // bisect kernel
  fileName   = oklFilePrefix + "Bisect" + suffix + oklFileSuffix;
  kernelName = "Bisect" + suffix;
  bisectKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);                                       



}
    
}
