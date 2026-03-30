#include "mesh.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace libp {

// =============================================================================
// Driver: preprocess all probes on Hex3D mesh
// =============================================================================
void mesh_t::ProbeHex3D() {

  const dlong Nprobes = probeX.length();

  LIBP_ABORT("ProbeHex3D: probeX, probeY and probeZ must have same length",
             probeY.length() != Nprobes || probeZ.length() != Nprobes);

  probeElement.malloc(Nprobes);
  probeR.malloc(Nprobes);
  probeS.malloc(Nprobes);
  probeT.malloc(Nprobes);

  for (dlong p = 0; p < Nprobes; ++p) {
    probeElement[p] = -1;
    probeR[p] = 0.0;
    probeS[p] = 0.0;
    probeT[p] = 0.0;

    ProbeLocateHex3D(probeX[p], probeY[p], probeZ[p],
                     probeElement[p],
                     probeR[p],
                     probeS[p],
                     probeT[p]);

    if (probeElement[p] < 0) {
      printf("WARNING: Probe %lld at (% .8e,% .8e,% .8e) not found in any element\n",
             (long long) p, probeX[p], probeY[p], probeZ[p]);
    }
  }

  probeInterp.malloc(Nprobes*Np);
  for (dlong i = 0; i < Nprobes*Np; ++i) probeInterp[i] = 0.0;

  BuildProbeInterpHex3D(probeR, probeS, probeT, probeElement, probeInterp);
}

// =============================================================================
// Locate element and compute reference coordinates (r,s,t)
// =============================================================================
void mesh_t::ProbeLocateHex3D(const dfloat pX,
                              const dfloat pY,
                              const dfloat pZ,
                              dlong &elem,
                              dfloat &rp,
                              dfloat &sp,
                              dfloat &tp) {

  const dfloat bboxTol = 1e-10;
  const dfloat rstol   = 1e-08;

  elem = -1;
  rp   = 0.0;
  sp   = 0.0;
  tp   = 0.0;

  for (dlong e = 0; e < Nelements; ++e) {

    const dlong id = e*Nverts;

    const dfloat x1 = EX[id+0];
    const dfloat x2 = EX[id+1];
    const dfloat x3 = EX[id+2];
    const dfloat x4 = EX[id+3];
    const dfloat x5 = EX[id+4];
    const dfloat x6 = EX[id+5];
    const dfloat x7 = EX[id+6];
    const dfloat x8 = EX[id+7];

    const dfloat y1 = EY[id+0];
    const dfloat y2 = EY[id+1];
    const dfloat y3 = EY[id+2];
    const dfloat y4 = EY[id+3];
    const dfloat y5 = EY[id+4];
    const dfloat y6 = EY[id+5];
    const dfloat y7 = EY[id+6];
    const dfloat y8 = EY[id+7];

    const dfloat z1 = EZ[id+0];
    const dfloat z2 = EZ[id+1];
    const dfloat z3 = EZ[id+2];
    const dfloat z4 = EZ[id+3];
    const dfloat z5 = EZ[id+4];
    const dfloat z6 = EZ[id+5];
    const dfloat z7 = EZ[id+6];
    const dfloat z8 = EZ[id+7];

    const dfloat xmin = std::min(std::min(std::min(x1,x2), std::min(x3,x4)),
                                 std::min(std::min(x5,x6), std::min(x7,x8)));
    const dfloat xmax = std::max(std::max(std::max(x1,x2), std::max(x3,x4)),
                                 std::max(std::max(x5,x6), std::max(x7,x8)));

    const dfloat ymin = std::min(std::min(std::min(y1,y2), std::min(y3,y4)),
                                 std::min(std::min(y5,y6), std::min(y7,y8)));
    const dfloat ymax = std::max(std::max(std::max(y1,y2), std::max(y3,y4)),
                                 std::max(std::max(y5,y6), std::max(y7,y8)));

    const dfloat zmin = std::min(std::min(std::min(z1,z2), std::min(z3,z4)),
                                 std::min(std::min(z5,z6), std::min(z7,z8)));
    const dfloat zmax = std::max(std::max(std::max(z1,z2), std::max(z3,z4)),
                                 std::max(std::max(z5,z6), std::max(z7,z8)));

    const bool inBBox =
      (pX >= xmin - bboxTol && pX <= xmax + bboxTol &&
       pY >= ymin - bboxTol && pY <= ymax + bboxTol &&
       pZ >= zmin - bboxTol && pZ <= zmax + bboxTol);

    if (!inBBox) continue;

    dfloat rtmp = 0.0;
    dfloat stmp = 0.0;
    dfloat ttmp = 0.0;

    const bool inside = PhysicalToReferenceHex3D(pX, pY, pZ, e, rtmp, stmp, ttmp);

    if (inside &&
        rtmp >= -1.0-rstol && rtmp <= 1.0+rstol &&
        stmp >= -1.0-rstol && stmp <= 1.0+rstol &&
        ttmp >= -1.0-rstol && ttmp <= 1.0+rstol) {
      elem = e;
      rp   = rtmp;
      sp   = stmp;
      tp   = ttmp;
      return;
    }
  }
}

// =============================================================================
// Invert trilinear mapping for one Hex8 element using Newton
// Returns true if the final (r,s,t) is inside reference element
// =============================================================================
bool mesh_t::PhysicalToReferenceHex3D(const dfloat pX,
                                      const dfloat pY,
                                      const dfloat pZ,
                                      const dlong e,
                                      dfloat &rp,
                                      dfloat &sp,
                                      dfloat &tp) {

  const dlong id = e*Nverts;

  const dfloat x1 = EX[id+0];
  const dfloat x2 = EX[id+1];
  const dfloat x3 = EX[id+2];
  const dfloat x4 = EX[id+3];
  const dfloat x5 = EX[id+4];
  const dfloat x6 = EX[id+5];
  const dfloat x7 = EX[id+6];
  const dfloat x8 = EX[id+7];

  const dfloat y1 = EY[id+0];
  const dfloat y2 = EY[id+1];
  const dfloat y3 = EY[id+2];
  const dfloat y4 = EY[id+3];
  const dfloat y5 = EY[id+4];
  const dfloat y6 = EY[id+5];
  const dfloat y7 = EY[id+6];
  const dfloat y8 = EY[id+7];

  const dfloat z1 = EZ[id+0];
  const dfloat z2 = EZ[id+1];
  const dfloat z3 = EZ[id+2];
  const dfloat z4 = EZ[id+3];
  const dfloat z5 = EZ[id+4];
  const dfloat z6 = EZ[id+5];
  const dfloat z7 = EZ[id+6];
  const dfloat z8 = EZ[id+7];

  dfloat r_m = 0.0;
  dfloat s_m = 0.0;
  dfloat t_m = 0.0;

  const int maxIter = 25;
  const dfloat tol = 1e-12;

  for (int iter = 0; iter < maxIter; ++iter) {

    // trilinear Hex8 shape functions
    const dfloat N1 = 0.125*(1.0-r_m)*(1.0-s_m)*(1.0-t_m);
    const dfloat N2 = 0.125*(1.0+r_m)*(1.0-s_m)*(1.0-t_m);
    const dfloat N3 = 0.125*(1.0+r_m)*(1.0+s_m)*(1.0-t_m);
    const dfloat N4 = 0.125*(1.0-r_m)*(1.0+s_m)*(1.0-t_m);
    const dfloat N5 = 0.125*(1.0-r_m)*(1.0-s_m)*(1.0+t_m);
    const dfloat N6 = 0.125*(1.0+r_m)*(1.0-s_m)*(1.0+t_m);
    const dfloat N7 = 0.125*(1.0+r_m)*(1.0+s_m)*(1.0+t_m);
    const dfloat N8 = 0.125*(1.0-r_m)*(1.0+s_m)*(1.0+t_m);

    const dfloat x_m =
      N1*x1 + N2*x2 + N3*x3 + N4*x4 +
      N5*x5 + N6*x6 + N7*x7 + N8*x8;

    const dfloat y_m =
      N1*y1 + N2*y2 + N3*y3 + N4*y4 +
      N5*y5 + N6*y6 + N7*y7 + N8*y8;

    const dfloat z_m =
      N1*z1 + N2*z2 + N3*z3 + N4*z4 +
      N5*z5 + N6*z6 + N7*z7 + N8*z8;

    const dfloat Rx = x_m - pX;
    const dfloat Ry = y_m - pY;
    const dfloat Rz = z_m - pZ;

    if (std::abs(Rx) < tol && std::abs(Ry) < tol && std::abs(Rz) < tol) {
      rp = r_m;
      sp = s_m;
      tp = t_m;

      const dfloat eps = 1e-8;
      return (rp >= -1.0-eps && rp <= 1.0+eps &&
              sp >= -1.0-eps && sp <= 1.0+eps &&
              tp >= -1.0-eps && tp <= 1.0+eps);
    }

    // derivatives of shape functions
    const dfloat dN1dr = -0.125*(1.0-s_m)*(1.0-t_m);
    const dfloat dN2dr =  0.125*(1.0-s_m)*(1.0-t_m);
    const dfloat dN3dr =  0.125*(1.0+s_m)*(1.0-t_m);
    const dfloat dN4dr = -0.125*(1.0+s_m)*(1.0-t_m);
    const dfloat dN5dr = -0.125*(1.0-s_m)*(1.0+t_m);
    const dfloat dN6dr =  0.125*(1.0-s_m)*(1.0+t_m);
    const dfloat dN7dr =  0.125*(1.0+s_m)*(1.0+t_m);
    const dfloat dN8dr = -0.125*(1.0+s_m)*(1.0+t_m);

    const dfloat dN1ds = -0.125*(1.0-r_m)*(1.0-t_m);
    const dfloat dN2ds = -0.125*(1.0+r_m)*(1.0-t_m);
    const dfloat dN3ds =  0.125*(1.0+r_m)*(1.0-t_m);
    const dfloat dN4ds =  0.125*(1.0-r_m)*(1.0-t_m);
    const dfloat dN5ds = -0.125*(1.0-r_m)*(1.0+t_m);
    const dfloat dN6ds = -0.125*(1.0+r_m)*(1.0+t_m);
    const dfloat dN7ds =  0.125*(1.0+r_m)*(1.0+t_m);
    const dfloat dN8ds =  0.125*(1.0-r_m)*(1.0+t_m);

    const dfloat dN1dt = -0.125*(1.0-r_m)*(1.0-s_m);
    const dfloat dN2dt = -0.125*(1.0+r_m)*(1.0-s_m);
    const dfloat dN3dt = -0.125*(1.0+r_m)*(1.0+s_m);
    const dfloat dN4dt = -0.125*(1.0-r_m)*(1.0+s_m);
    const dfloat dN5dt =  0.125*(1.0-r_m)*(1.0-s_m);
    const dfloat dN6dt =  0.125*(1.0+r_m)*(1.0-s_m);
    const dfloat dN7dt =  0.125*(1.0+r_m)*(1.0+s_m);
    const dfloat dN8dt =  0.125*(1.0-r_m)*(1.0+s_m);

    // Jacobian entries
    const dfloat dxdr =
      dN1dr*x1 + dN2dr*x2 + dN3dr*x3 + dN4dr*x4 +
      dN5dr*x5 + dN6dr*x6 + dN7dr*x7 + dN8dr*x8;
    const dfloat dxds =
      dN1ds*x1 + dN2ds*x2 + dN3ds*x3 + dN4ds*x4 +
      dN5ds*x5 + dN6ds*x6 + dN7ds*x7 + dN8ds*x8;
    const dfloat dxdt =
      dN1dt*x1 + dN2dt*x2 + dN3dt*x3 + dN4dt*x4 +
      dN5dt*x5 + dN6dt*x6 + dN7dt*x7 + dN8dt*x8;

    const dfloat dydr =
      dN1dr*y1 + dN2dr*y2 + dN3dr*y3 + dN4dr*y4 +
      dN5dr*y5 + dN6dr*y6 + dN7dr*y7 + dN8dr*y8;
    const dfloat dyds =
      dN1ds*y1 + dN2ds*y2 + dN3ds*y3 + dN4ds*y4 +
      dN5ds*y5 + dN6ds*y6 + dN7ds*y7 + dN8ds*y8;
    const dfloat dydt =
      dN1dt*y1 + dN2dt*y2 + dN3dt*y3 + dN4dt*y4 +
      dN5dt*y5 + dN6dt*y6 + dN7dt*y7 + dN8dt*y8;

    const dfloat dzdr =
      dN1dr*z1 + dN2dr*z2 + dN3dr*z3 + dN4dr*z4 +
      dN5dr*z5 + dN6dr*z6 + dN7dr*z7 + dN8dr*z8;
    const dfloat dzds =
      dN1ds*z1 + dN2ds*z2 + dN3ds*z3 + dN4ds*z4 +
      dN5ds*z5 + dN6ds*z6 + dN7ds*z7 + dN8ds*z8;
    const dfloat dzdt =
      dN1dt*z1 + dN2dt*z2 + dN3dt*z3 + dN4dt*z4 +
      dN5dt*z5 + dN6dt*z6 + dN7dt*z7 + dN8dt*z8;

    // Solve J * [dr ds dt]^T = -R
    const dfloat detJ =
        dxdr*(dyds*dzdt - dydt*dzds)
      - dxds*(dydr*dzdt - dydt*dzdr)
      + dxdt*(dydr*dzds - dyds*dzdr);

    if (std::abs(detJ) < 1e-16) {
      rp = r_m;
      sp = s_m;
      tp = t_m;
      return false;
    }

    const dfloat b1 = -Rx;
    const dfloat b2 = -Ry;
    const dfloat b3 = -Rz;

    const dfloat dr =
        ( b1*(dyds*dzdt - dydt*dzds)
        - dxds*(b2*dzdt - dydt*b3)
        + dxdt*(b2*dzds - dyds*b3) ) / detJ;

    const dfloat ds =
        ( dxdr*(b2*dzdt - dydt*b3)
        - b1*(dydr*dzdt - dydt*dzdr)
        + dxdt*(dydr*b3 - b2*dzdr) ) / detJ;

    const dfloat dt =
        ( dxdr*(dyds*b3 - b2*dzds)
        - dxds*(dydr*b3 - b2*dzdr)
        + b1*(dydr*dzds - dyds*dzdr) ) / detJ;

    r_m += dr;
    s_m += ds;
    t_m += dt;

    if (std::abs(dr) < tol && std::abs(ds) < tol && std::abs(dt) < tol) {
      rp = r_m;
      sp = s_m;
      tp = t_m;

      const dfloat eps = 1e-8;
      return (rp >= -1.0-eps && rp <= 1.0+eps &&
              sp >= -1.0-eps && sp <= 1.0+eps &&
              tp >= -1.0-eps && tp <= 1.0+eps);
    }
  }

  rp = r_m;
  sp = s_m;
  tp = t_m;

  const dfloat eps = 1e-6;
  return (rp >= -1.0-eps && rp <= 1.0+eps &&
          sp >= -1.0-eps && sp <= 1.0+eps &&
          tp >= -1.0-eps && tp <= 1.0+eps);
}

// =============================================================================
// Build flattened interpolation matrix for Hex3D
//   probeInterp[p*Np + n] = value of basis/node n at (probeR[p], probeS[p], probeT[p])
// =============================================================================
void mesh_t::BuildProbeInterpHex3D(const memory<dfloat> &probeR_,
                                   const memory<dfloat> &probeS_,
                                   const memory<dfloat> &probeT_,
                                   const memory<dlong>  &probeElement_,
                                   memory<dfloat> &probeInterp_) {

  const dlong Nprobes = probeR_.length();

  LIBP_ABORT("BuildProbeInterpHex3D: inconsistent probe array lengths",
             probeS_.length() != Nprobes ||
             probeT_.length() != Nprobes ||
             probeElement_.length() != Nprobes);

  LIBP_ABORT("BuildProbeInterpHex3D: probeInterp has wrong size",
             probeInterp_.length() != Nprobes*Np);

  memory<dfloat> rProbe(1), sProbe(1), tProbe(1), I1(Np);

  for (dlong p = 0; p < Nprobes; ++p) {

    for (int n = 0; n < Np; ++n) {
      probeInterp_[p*Np + n] = 0.0;
    }

    if (probeElement_[p] < 0) continue;

    rProbe[0] = probeR_[p];
    sProbe[0] = probeS_[p];
    tProbe[0] = probeT_[p];

    // Expected behavior:
    // InterpolationMatrixHex3D(N, r, s, t, rProbe, sProbe, tProbe, I1)
    InterpolationMatrixHex3D(N, r, s, t, rProbe, sProbe, tProbe, I1);

    for (int n = 0; n < Np; ++n) {
      probeInterp_[p*Np + n] = I1[n];
    }
  }
}

} // namespace libp