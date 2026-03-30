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
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace libp {

// =============================================================================
// Driver: preprocess all probes on Quad2D mesh
// =============================================================================
void mesh_t::ProbeQuad2D() {

  const dlong Nprobes = probeX.length();

    LIBP_ABORT("ProbeQuad2D: probeX and probeY must have same length",probeY.length() != Nprobes);
  

  probeElement.malloc(Nprobes);
  probeR.malloc(Nprobes);
  probeS.malloc(Nprobes);

  for (dlong p = 0; p < Nprobes; ++p) {
    probeElement[p] = -1;
    probeR[p] = 0.0;
    probeS[p] = 0.0;

    ProbeLocateQuad2D(probeX[p], probeY[p],
                      probeElement[p],
                      probeR[p],
                      probeS[p]);

    if (probeElement[p] < 0) {
      printf("WARNING: Probe %lld at (% .8e,% .8e) not found in any element\n",
             (long long) p, probeX[p], probeY[p]);
    }
  }

  // flattened matrix of size Nprobes x Np
  probeInterp.malloc(Nprobes*Np);
  for (dlong i = 0; i < Nprobes*Np; ++i) probeInterp[i] = 0.0;

  BuildProbeInterpQuad2D(probeR, probeS, probeElement, probeInterp);
}

// =============================================================================
// Locate element and compute reference coordinates (r,s)
// =============================================================================
void mesh_t::ProbeLocateQuad2D(const dfloat pX,
                               const dfloat pY,
                               dlong &elem,
                               dfloat &rp,
                               dfloat &sp) {

  const dfloat bboxTol = 1e-10;
  const dfloat rstol   = 1e-08;

  elem = -1;
  rp   = 0.0;
  sp   = 0.0;

  for (dlong e = 0; e < Nelements; ++e) {

    const dlong id = e*Nverts;

    const dfloat x1 = EX[id+0];
    const dfloat x2 = EX[id+1];
    const dfloat x3 = EX[id+2];
    const dfloat x4 = EX[id+3];

    const dfloat y1 = EY[id+0];
    const dfloat y2 = EY[id+1];
    const dfloat y3 = EY[id+2];
    const dfloat y4 = EY[id+3];

    // cheap bounding-box prefilter
    const dfloat xmin = std::min(std::min(x1,x2), std::min(x3,x4));
    const dfloat xmax = std::max(std::max(x1,x2), std::max(x3,x4));
    const dfloat ymin = std::min(std::min(y1,y2), std::min(y3,y4));
    const dfloat ymax = std::max(std::max(y1,y2), std::max(y3,y4));

    const bool inBBox =
      (pX >= xmin - bboxTol && pX <= xmax + bboxTol &&
       pY >= ymin - bboxTol && pY <= ymax + bboxTol);

    if (!inBBox) continue;

    dfloat rtmp = 0.0;
    dfloat stmp = 0.0;

    const bool inside = PhysicalToReferenceQuad2D(pX, pY, e, rtmp, stmp);

    if (inside &&
        rtmp >= -1.0-rstol && rtmp <= 1.0+rstol &&
        stmp >= -1.0-rstol && stmp <= 1.0+rstol) {
      elem = e;
      rp   = rtmp;
      sp   = stmp;
      return;
    }
  }
}

// =============================================================================
// Invert bilinear mapping for one Quad4 element using Newton
// Returns true if the final (r,s) is inside reference element
// =============================================================================
bool mesh_t::PhysicalToReferenceQuad2D(const dfloat pX,
                                       const dfloat pY,
                                       const dlong e,
                                       dfloat &rp,
                                       dfloat &sp) {

  const dlong id = e*Nverts;

  const dfloat x1 = EX[id+0];
  const dfloat x2 = EX[id+1];
  const dfloat x3 = EX[id+2];
  const dfloat x4 = EX[id+3];

  const dfloat y1 = EY[id+0];
  const dfloat y2 = EY[id+1];
  const dfloat y3 = EY[id+2];
  const dfloat y4 = EY[id+3];

  dfloat r_m = 0.0;
  dfloat s_m = 0.0;

  const int maxIter = 20;
  const dfloat tol = 1e-12;

  for (int iter = 0; iter < maxIter; ++iter) {

    // bilinear shape functions
    const dfloat N1 = 0.25*(1.0-r_m)*(1.0-s_m);
    const dfloat N2 = 0.25*(1.0+r_m)*(1.0-s_m);
    const dfloat N3 = 0.25*(1.0+r_m)*(1.0+s_m);
    const dfloat N4 = 0.25*(1.0-r_m)*(1.0+s_m);

    // mapped point
    const dfloat x_m =
      N1*x1 + N2*x2 + N3*x3 + N4*x4;

    const dfloat y_m =
      N1*y1 + N2*y2 + N3*y3 + N4*y4;

    // residual
    const dfloat Rx = x_m - pX;
    const dfloat Ry = y_m - pY;

    if (std::abs(Rx) < tol && std::abs(Ry) < tol) {
      rp = r_m;
      sp = s_m;

      const dfloat eps = 1e-8;
      return (rp >= -1.0-eps && rp <= 1.0+eps &&
              sp >= -1.0-eps && sp <= 1.0+eps);
    }

    // shape function derivatives
    const dfloat dN1dr = -0.25*(1.0-s_m);
    const dfloat dN2dr =  0.25*(1.0-s_m);
    const dfloat dN3dr =  0.25*(1.0+s_m);
    const dfloat dN4dr = -0.25*(1.0+s_m);

    const dfloat dN1ds = -0.25*(1.0-r_m);
    const dfloat dN2ds = -0.25*(1.0+r_m);
    const dfloat dN3ds =  0.25*(1.0+r_m);
    const dfloat dN4ds =  0.25*(1.0-r_m);

    // Jacobian
    const dfloat dxdr =
      dN1dr*x1 + dN2dr*x2 + dN3dr*x3 + dN4dr*x4;
    const dfloat dxds =
      dN1ds*x1 + dN2ds*x2 + dN3ds*x3 + dN4ds*x4;

    const dfloat dydr =
      dN1dr*y1 + dN2dr*y2 + dN3dr*y3 + dN4dr*y4;
    const dfloat dyds =
      dN1ds*y1 + dN2ds*y2 + dN3ds*y3 + dN4ds*y4;

    const dfloat detJ = dxdr*dyds - dxds*dydr;

    if (std::abs(detJ) < 1e-16) {
      rp = r_m;
      sp = s_m;
      return false;
    }

    // Newton update
    const dfloat dr = (-dyds*Rx + dxds*Ry)/detJ;
    const dfloat ds = ( dydr*Rx - dxdr*Ry)/detJ;

    r_m += dr;
    s_m += ds;

    if (std::abs(dr) < tol && std::abs(ds) < tol) {
      rp = r_m;
      sp = s_m;

      const dfloat eps = 1e-8;
      return (rp >= -1.0-eps && rp <= 1.0+eps &&
              sp >= -1.0-eps && sp <= 1.0+eps);
    }
  }

  rp = r_m;
  sp = s_m;

  const dfloat eps = 1e-6;
  return (rp >= -1.0-eps && rp <= 1.0+eps &&
          sp >= -1.0-eps && sp <= 1.0+eps);
}

// =============================================================================
// Build flattened interpolation matrix:
//   probeInterp[p*Np + n] = value of basis/node n at (probeR[p], probeS[p])
//
// For probes not found in mesh (probeElement[p] < 0), row is left zero.
// =============================================================================
void mesh_t::BuildProbeInterpQuad2D(const memory<dfloat> &probeR_,
                                    const memory<dfloat> &probeS_,
                                    const memory<dlong>  &probeElement_,
                                    memory<dfloat> &probeInterp_) {

  const dlong Nprobes = probeR_.length();

  LIBP_ABORT("BuildProbeInterpQuad2D: inconsistent probe array lengths",
              probeS_.length() != Nprobes || probeElement_.length() != Nprobes);

  LIBP_ABORT("BuildProbeInterpQuad2D: probeInterp has wrong size",
              probeInterp_.length() != Nprobes*Np);

  // temporary 1-point arrays for calling interpolation routine
  memory<dfloat> rProbe(1), sProbe(1), I1(Np);

  for (dlong p = 0; p < Nprobes; ++p) {

    // zero row by default
    for (int n = 0; n < Np; ++n) {
      probeInterp_[p*Np + n] = 0.0;
    }

    if (probeElement_[p] < 0) continue;

    rProbe[0] = probeR_[p];
    sProbe[0] = probeS_[p];

    // Expected behavior:
    // InterpolationMatrixQuad2D(N, r, s, rProbe, sProbe, I1)
    // fills I1[n] = basis/node n evaluated at (rProbe[0], sProbe[0])
    InterpolationMatrixQuad2D(N, r, s, rProbe, sProbe, I1);

    for (int n = 0; n < Np; ++n) {
      probeInterp_[p*Np + n] = I1[n];
    }
  }
}

} // namespace libp