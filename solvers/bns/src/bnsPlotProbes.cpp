#include "bns.hpp"
#include <cstdio>
#include <cmath>

void bns_t::ReportProbes(memory<dfloat>& Q,
                         memory<dfloat>& V,
                         const dfloat time,
                         const std::string& fileName,
                         const bool writeHeader) {

  const dlong Nprobes = mesh.probeElement.length();

  if (mesh.probeX.length() != Nprobes || mesh.probeY.length() != Nprobes) {
    printf("ERROR: ReportProbes: inconsistent probe coordinate array sizes\n");
    return;
  }

  if (mesh.probeInterp.length() != Nprobes*mesh.Np) {
    printf("ERROR: ReportProbes: interpolation array has wrong size\n");
    return;
  }

  if (mesh.dim == 3 && mesh.probeZ.length() != Nprobes) {
    printf("ERROR: ReportProbes: probeZ has wrong size for 3D case\n");
    return;
  }

  FILE* fp = fopen(fileName.c_str(), writeHeader ? "w" : "a");
  if (!fp) {
    printf("ERROR: could not open probe output file: %s\n", fileName.c_str());
    return;
  }

  if (writeHeader) {
    fprintf(fp, "# time probeId x y z rho p u v w vortx vorty vortz\n");
  }

  for (dlong pId = 0; pId < Nprobes; ++pId) {
    const dlong e = mesh.probeElement[pId];

    if (e < 0 || e >= mesh.Nelements) continue;

    dfloat rhoP  = 0.0;
    dfloat rhouP = 0.0;
    dfloat rhovP = 0.0;
    dfloat rhowP = 0.0;

    dfloat vortX = 0.0;
    dfloat vortY = 0.0;
    dfloat vortZ = 0.0;

    for (int n = 0; n < mesh.Np; ++n) {
      const dfloat wgt  = mesh.probeInterp[pId*mesh.Np + n];
      const dlong baseQ = e*mesh.Np*Nfields + n;

      rhoP  += wgt * Q[baseQ + 0*mesh.Np];
      rhouP += wgt * Q[baseQ + 1*mesh.Np];
      rhovP += wgt * Q[baseQ + 2*mesh.Np];

      if (mesh.dim == 3) {
        rhowP += wgt * Q[baseQ + 3*mesh.Np];
      }
    }

    dfloat uP = 0.0;
    dfloat vP = 0.0;
    dfloat wP = 0.0;
    dfloat pP = 0.0;

    if (rhoP > 1e-14) {
      uP = c * rhouP / rhoP;
      vP = c * rhovP / rhoP;
      if (mesh.dim == 3) wP = c * rhowP / rhoP;
      pP = RT * rhoP;
    }

    if (V.length() != 0) {
      if (mesh.dim == 2) {
        dfloat vortScalar = 0.0;
        for (int n = 0; n < mesh.Np; ++n) {
          const dfloat wgt = mesh.probeInterp[pId*mesh.Np + n];
          vortScalar += wgt * V[e*mesh.Np + n];
        }
        vortX = 0.0;
        vortY = 0.0;
        vortZ = vortScalar;
      } else {
        for (int n = 0; n < mesh.Np; ++n) {
          const dfloat wgt  = mesh.probeInterp[pId*mesh.Np + n];
          const dlong baseV = e*mesh.Np*3 + n;

          vortX += wgt * V[baseV + 0*mesh.Np];
          vortY += wgt * V[baseV + 1*mesh.Np];
          vortZ += wgt * V[baseV + 2*mesh.Np];
        }
      }
    }

    const dfloat zP = (mesh.dim == 3) ? mesh.probeZ[pId] : 0.0;

    fprintf(fp,
            "%.16e %lld %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            time,
            (long long) pId,
            mesh.probeX[pId],
            mesh.probeY[pId],
            zP,
            rhoP,
            pP,
            uP,
            vP,
            wP,
            vortX,
            vortY,
            vortZ);
  }

  fclose(fp);
}