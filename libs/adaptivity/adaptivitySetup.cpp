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
void adaptivity_t::Setup(platform_t& _platform, adaptivitySettings_t& _settings,
                   comm_t _comm){

  platform = _platform;
  settings = _settings;
  mesh = _mesh;
  props = platform.props();


  comm = _comm.Dup();
  rank = comm.rank();
  size = comm.size();
  
    // compute interpolation matrices for amr
InterpolateToChildTri2D();
InterpolateToParentTri2D();

// OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

    // indicator kernel
fileName   = oklFilePrefix + "bnsIndicator" + suffix + oklFileSuffix;
kernelName = "bnsIndicatorTest" + suffix;
indicatorKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);
  // combine solution kernel
fileName   = oklFilePrefix + "bnsCombine" + suffix + oklFileSuffix;
kernelName = "bnsCombine" + suffix;
combineKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);
  // split solution kernel
fileName   = oklFilePrefix + "bnsSplit" + suffix + oklFileSuffix;
kernelName = "bnsSplit" + suffix;
splitKernel = platform.buildKernel(fileName, kernelName,
                                   kernelInfo);
}
    
}
