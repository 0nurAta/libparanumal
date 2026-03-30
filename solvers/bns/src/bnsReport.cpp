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

#include "bns.hpp"
#include <cmath>

void bns_t::Report(dfloat time, int tstep){

  static int frame = 0;
  static int forceFrame = 0;
  static int surfaceFrame = 0;
  static bool probeHeaderWritten = false;

  // write local surface coefficients every N report calls
  static const int surfaceOutputInterval = 50;

  // --------------------------------------------------------------------------
  // compute ||q||_M
  // --------------------------------------------------------------------------
  const dlong Nentries = mesh.Nelements * mesh.Np * Nfields;

  deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
  mesh.MassMatrixApply(o_q, o_Mq);

  const dfloat norm2 =
    std::sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

  o_Mq.free();

  if (mesh.rank == 0) {
    printf("%5.4f (%d), %5.4f (time, timestep, norm)\n", time, tstep, norm2);
  }

  // --------------------------------------------------------------------------
  // reporting choices
  // --------------------------------------------------------------------------
  const bool writeFields = false;              // no regular field output
  const bool writeProbes = mesh.probesEnabled; // probes every report step
  const bool writeForcesNow = true;            // forces every report step
  const bool writeSurfaceNow = (surfaceFrame % surfaceOutputInterval == 0);

  const bool writeQCritSnapshot =
      std::abs(time - 2000.0) < 1e-6 ||
      std::abs(time - 4000.0) < 1e-6 ||
      std::abs(time - 9000.0) < 1e-6;

  const bool needVorticity = writeFields || writeProbes;
  const bool needQCrit     = writeQCritSnapshot;
  const bool needHostQ     = writeFields || writeProbes || writeForcesNow || writeSurfaceNow || needQCrit;

  if (!needHostQ) return;

  // --------------------------------------------------------------------------
  // copy solution to host
  // --------------------------------------------------------------------------
  o_q.copyTo(q);

  // --------------------------------------------------------------------------
  // compute vorticity only if needed
  // --------------------------------------------------------------------------
  memory<dfloat> Vort;
  if (needVorticity) {
    deviceMemory<dfloat> o_Vort =
      platform.reserve<dfloat>(mesh.dim * mesh.Nelements * mesh.Np);

    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    Vort.malloc(mesh.dim * mesh.Nelements * mesh.Np);
    o_Vort.copyTo(Vort);
    o_Vort.free();
  }

  // --------------------------------------------------------------------------
  // compute Q-criterion only at selected snapshot times
  // --------------------------------------------------------------------------
  memory<dfloat> QCrit;
  if (needQCrit) {
    deviceMemory<dfloat> o_QCrit =
      platform.reserve<dfloat>(mesh.Nelements * mesh.Np);

    qcriterionKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_QCrit);

    QCrit.malloc(mesh.Nelements * mesh.Np);
    o_QCrit.copyTo(QCrit);
    o_QCrit.free();
  }

  // --------------------------------------------------------------------------
  // file base name
  // --------------------------------------------------------------------------
  std::string outName;
  settings.getSetting("OUTPUT FILE NAME", outName);

  // --------------------------------------------------------------------------
  // regular field output disabled
  // --------------------------------------------------------------------------
  if (writeFields) {
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", outName.c_str(), mesh.rank, frame++);
    PlotFields(q, Vort, std::string(fname));
  }

  // --------------------------------------------------------------------------
  // Q-criterion snapshots only at selected times
  // --------------------------------------------------------------------------
  if (writeQCritSnapshot) {
    char fname[BUFSIZ];
    sprintf(fname, "%s_qcrit_%04d_%04d.vtu", outName.c_str(), mesh.rank, frame++);
    PlotFields(q, QCrit, std::string(fname));
  }

  // --------------------------------------------------------------------------
  // integrated coefficients every report step
  // --------------------------------------------------------------------------
  if (writeForcesNow) {
    writeForces(time, tstep, forceFrame);
  }

  // --------------------------------------------------------------------------
  // local surface coefficients every surfaceOutputInterval report calls
  // --------------------------------------------------------------------------
  if (writeSurfaceNow) {
    writeSurfaceCoeffs(time, surfaceFrame);
  }

  forceFrame++;
  surfaceFrame++;

  // --------------------------------------------------------------------------
  // probes every report step
  // --------------------------------------------------------------------------
  if (writeProbes) {
    char pfname[BUFSIZ];
    sprintf(pfname, "%s_%04d.dat", mesh.probeOutputFile.c_str(), mesh.rank);

    ReportProbes(q, Vort, time, std::string(pfname), !probeHeaderWritten);
    probeHeaderWritten = true;
  }



 


  /*
  if(bns->dim==3){
    if(options.compareArgs("OUTPUT FILE FORMAT","ISO")){

      for (int gr=0; gr<bns->isoGNgroups; gr++){

        bns->isoNtris[0] = 0;
        bns->o_isoNtris.copyFrom(bns->isoNtris);
        if(mesh->nonPmlNelements){
        bns->isoSurfaceKernel(mesh->nonPmlNelements,    // Numner of elements
                              mesh->o_nonPmlElementIds,    // Element Ids
                              bns->isoField,               // which field to use for isosurfacing
                              bns->isoColorField,          // which field to use for isosurfacing
                              bns->isoGNlevels[gr],        // number of isosurface levels
                              bns->o_isoGLvalues[gr],      // array of isosurface levels
                              bns->isoMaxNtris,            // maximum number of generated triangles
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              bns->o_q,
                              bns->o_Vort,
                              bns->o_VortMag,
                              bns->o_plotInterp,
                              bns->o_plotEToV,
                              bns->o_isoNtris,             // output: number of generated triangles
                              bns->o_isoq);                // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)

      }
        // find number of generated triangles
        bns->o_isoNtris.copyTo(bns->isoNtris);
        bns->isoNtris[0] = mymin(bns->isoNtris[0], bns->isoMaxNtris);

        //
        printf("Rank:%2d Group:%2d Triangles:%8d\n", mesh->rank, bns->isoNtris[0], gr);
        //
        int offset = 0;
        bns->o_isoq.copyTo(bns->isoq, bns->isoNtris[0]*(mesh->dim+bns->isoNfields)*3*sizeof(dfloat), offset);

        char fname[BUFSIZ];
        string outName;
        options.getArgs("OUTPUT FILE NAME", outName);


        if(options.compareArgs("OUTPUT FILE FORMAT", "WELD"))
        {
          int Ntris1 = bns->isoNtris[0];
          int Ntris2 = bnsWeldTriVerts(bns, Ntris1, bns->isoq);

          printf("Welding triangles:%8d to:%8d\n", Ntris1, Ntris2);
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoWeldPlotVTU(bns,  fname);
        }
        else
        {
          sprintf(fname, "%s_%d_%d_ %04d_%04d.vtu",(char*)outName.c_str(), bns->isoField, gr, mesh->rank, bns->frame);
          bnsIsoPlotVTU(bns, bns->isoNtris[0], bns->isoq, fname);
        }
      }
      bns->frame++;
    }
  }
  */
}

void bns_t::writeForces(dfloat time, int tstep, int frame){

  dfloat rref = 1.0;
  dfloat uref = 0.1;
  dfloat vref = 0.0;
  dfloat wref = 0.0;

  const dfloat velRef = (mesh.dim == 2)
      ? std::sqrt(uref*uref + vref*vref)
      : std::sqrt(uref*uref + vref*vref + wref*wref);

  dfloat Aref = 1.0; // reference length in 2D, reference area in 3D

  LIBP_ABORT("writeForces: reference velocity magnitude is zero", velRef <= 0);
  LIBP_ABORT("writeForces: reference area/length is zero", Aref <= 0);
  LIBP_ABORT("writeForces: reference density is zero", rref <= 0);

  const dfloat rcp_dynp = 1.0/(0.5*rref*velRef*velRef*Aref);

  if (mesh.rank == 0) {
    printf("----------------------------------------------------------------------\n");
    if (mesh.dim == 2) printf("Drag\t\tLift\n");
    else               printf("X-Forces\tY-Forces\tZ-Forces\n");
  }

  const dlong NforceFields = (mesh.dim == 2) ? 4 : 6;
  const dlong Nentries = mesh.Nelements * mesh.Np * NforceFields;

  deviceMemory<dfloat> o_F = platform.reserve<dfloat>(Nentries);

  forceKernel(mesh.Nelements, c,
              mesh.o_sgeo, mesh.o_sM, mesh.o_vmapM, mesh.o_EToB,
              mesh.o_x, mesh.o_y, mesh.o_z,
              o_q, o_F);

  const dlong shift = mesh.Nelements * mesh.Np;

  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name += "_analysis.dat";

  if (mesh.dim == 2){
    const dfloat vFx = rcp_dynp*platform.linAlg().sum(shift, o_F + 0*shift, mesh.comm);
    const dfloat vFy = rcp_dynp*platform.linAlg().sum(shift, o_F + 1*shift, mesh.comm);
    const dfloat pFx = rcp_dynp*platform.linAlg().sum(shift, o_F + 2*shift, mesh.comm);
    const dfloat pFy = rcp_dynp*platform.linAlg().sum(shift, o_F + 3*shift, mesh.comm);

    if (mesh.rank == 0){
      FILE* fp = fopen(name.c_str(), (frame == 0) ? "w" : "a");

      if (frame == 0) fprintf(fp, "time Cd Cf_x Cp_x Cl Cf_y Cp_y\n");

      const dfloat Fx = vFx + pFx;
      const dfloat Fy = vFy + pFy;

      printf("%.6e \t %.6e\n", Fx, Fy);
      fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
              time, Fx, vFx, pFx, Fy, vFy, pFy);
      fclose(fp);

      printf("----------------------------------------------------------------------\n");
    }

  } else {
    const dfloat vFx = rcp_dynp*platform.linAlg().sum(shift, o_F + 0*shift, mesh.comm);
    const dfloat vFy = rcp_dynp*platform.linAlg().sum(shift, o_F + 1*shift, mesh.comm);
    const dfloat vFz = rcp_dynp*platform.linAlg().sum(shift, o_F + 2*shift, mesh.comm);
    const dfloat pFx = rcp_dynp*platform.linAlg().sum(shift, o_F + 3*shift, mesh.comm);
    const dfloat pFy = rcp_dynp*platform.linAlg().sum(shift, o_F + 4*shift, mesh.comm);
    const dfloat pFz = rcp_dynp*platform.linAlg().sum(shift, o_F + 5*shift, mesh.comm);

    if (mesh.rank == 0){
      FILE* fp = fopen(name.c_str(), (frame == 0) ? "w" : "a");

      if (frame == 0) fprintf(fp, "time Fx Fvx Fpx Fy Fvy Fpy Fz Fvz Fpz\n");

      const dfloat Fx = vFx + pFx;
      const dfloat Fy = vFy + pFy;
      const dfloat Fz = vFz + pFz;

      printf("%.6e \t %.6e \t %.6e\n", Fx, Fy, Fz);
      fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
              time, Fx, vFx, pFx, Fy, vFy, pFy, Fz, vFz, pFz);
      fclose(fp);

      printf("----------------------------------------------------------------------\n");
    }
  }
}


void bns_t::writeSurfaceCoeffs(dfloat time, int frame){
  if (mesh.dim == 2) {
    writeSurfaceCoeffsQuad2D(time, frame);
  } else if (mesh.dim == 3) {
    writeSurfaceCoeffsHex3D(time, frame);
  } else {
    LIBP_ABORT("writeSurfaceCoeffs: unsupported mesh dimension", 1);
  }
}

void bns_t::writeSurfaceCoeffsQuad2D(dfloat time, int frame){

  const dfloat rref = 1.0;
  const dfloat uref = 0.1;
  const dfloat vref = 0.0;

  const dfloat velRef = std::sqrt(uref*uref + vref*vref);
  LIBP_ABORT("writeSurfaceCoeffsQuad2D: reference velocity magnitude is zero",
             velRef <= 0.0);

  const dfloat pInf = rref*c*c;
  const dfloat dynp = 0.5*rref*velRef*velRef;

  o_q.copyTo(q);

  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name += "_surfcoeffs_rank";
  name += std::to_string(mesh.rank);
  name += "_";
  name += std::to_string(frame);
  name += ".dat";

  FILE* fp = fopen(name.c_str(), (frame == 0) ? "w" : "a");
  if (!fp) {
    printf("ERROR: writeSurfaceCoeffsQuad2D: could not open file %s\n", name.c_str());
    return;
  }

  if (frame == 0) {
    fprintf(fp, "# time element face node x y nx ny tx ty rho p Cp tauw Cf\n");
  }
   
   const int NXID  = 0;
   const int NYID  = 1;
   
  for (dlong e = 0; e < mesh.Nelements; ++e) {
    for (int f = 0; f < mesh.Nfaces; ++f) {

      const int bc = mesh.EToB[e*mesh.Nfaces + f];
      if (bc != 1) continue;

      for (int n = 0; n < mesh.Nfp; ++n) {

        const dlong sk = e*mesh.Nfaces*mesh.Nfp + f*mesh.Nfp + n;

        const dfloat nx = mesh.sgeo[sk*mesh.Nsgeo + NXID];
        const dfloat ny = mesh.sgeo[sk*mesh.Nsgeo + NYID];

        const dfloat tx = -ny;
        const dfloat ty =  nx;

        const dlong id   = mesh.vmapM[sk];
        const int   vid  = id % mesh.Np;
        const dlong base = e*mesh.Np*Nfields + vid;

        const dfloat q1 = q[base + 0*mesh.Np];
        const dfloat q2 = q[base + 1*mesh.Np];
        const dfloat q3 = q[base + 2*mesh.Np];
        const dfloat q4 = q[base + 3*mesh.Np];
        const dfloat q5 = q[base + 4*mesh.Np];
        const dfloat q6 = q[base + 5*mesh.Np];

        if (q1 <= 1e-14) continue;

        const dfloat p  = q1*c*c;
        const dfloat Cp = (p - pInf)/dynp;

        const dfloat s11 = -c*c*(std::sqrt(2.0)*q5 - q2*q2/q1);
        const dfloat s22 = -c*c*(std::sqrt(2.0)*q6 - q3*q3/q1);
        const dfloat s12 = -c*c*(q4 - q2*q3/q1);

        const dfloat tvx = -(s11*nx + s12*ny);
        const dfloat tvy = -(s12*nx + s22*ny);

        const dfloat tauw = tvx*tx + tvy*ty;
        const dfloat Cf   = tauw/dynp;

        const dlong nodeId = e*mesh.Np + vid;
        const dfloat x = mesh.x[nodeId];
        const dfloat y = mesh.y[nodeId];

        fprintf(fp,
                "%.16e %lld %d %d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                time,
                (long long)e,
                f,
                n,
                x,
                y,
                nx,
                ny,
                tx,
                ty,
                q1,
                p,
                Cp,
                tauw,
                Cf);
      }
    }
  }

  fclose(fp);
}

void bns_t::writeSurfaceCoeffsHex3D(dfloat time, int frame){

  const dfloat rref = 1.0;
  const dfloat uref = 0.1;
  const dfloat vref = 0.0;
  const dfloat wref = 0.0;

  const dfloat velRef = std::sqrt(uref*uref + vref*vref + wref*wref);
  LIBP_ABORT("writeSurfaceCoeffsHex3D: reference velocity magnitude is zero",
             velRef <= 0.0);

  const dfloat pInf = rref*c*c;
  const dfloat dynp = 0.5*rref*velRef*velRef;

  o_q.copyTo(q);

  std::string name;
  settings.getSetting("OUTPUT FILE NAME", name);
  name += "_surfcoeffs_rank";
  name += std::to_string(mesh.rank);
  name += "_";
  name += std::to_string(frame);
  name += ".dat";

  FILE* fp = fopen(name.c_str(), "w");
  if (!fp) {
    printf("ERROR: writeSurfaceCoeffsHex3D: could not open file %s\n", name.c_str());
    return;
  }

  fprintf(fp,
          "# time element face node x y z nx ny nz rho p Cp "
          "tvx tvy tvz tauwx tauwy tauwz tauw Cf\n");

  const int NXID = 0;
  const int NYID = 1;
  const int NZID = 2;

  for (dlong e = 0; e < mesh.Nelements; ++e) {
    for (int f = 0; f < mesh.Nfaces; ++f) {

      const int bc = mesh.EToB[e*mesh.Nfaces + f];
      if (bc != 1) continue;

      for (int n = 0; n < mesh.Nfp; ++n) {

        const dlong sk = e*mesh.Nfaces*mesh.Nfp + f*mesh.Nfp + n;

        const dfloat nx = mesh.sgeo[sk*mesh.Nsgeo + NXID];
        const dfloat ny = mesh.sgeo[sk*mesh.Nsgeo + NYID];
        const dfloat nz = mesh.sgeo[sk*mesh.Nsgeo + NZID];

        const dlong id   = mesh.vmapM[sk];
        const int   vid  = id % mesh.Np;
        const dlong base = e*mesh.Np*Nfields + vid;

        const dfloat q1  = q[base + 0*mesh.Np];
        const dfloat q2  = q[base + 1*mesh.Np];
        const dfloat q3  = q[base + 2*mesh.Np];
        const dfloat q4  = q[base + 3*mesh.Np];
        const dfloat q5  = q[base + 4*mesh.Np];
        const dfloat q6  = q[base + 5*mesh.Np];
        const dfloat q7  = q[base + 6*mesh.Np];
        const dfloat q8  = q[base + 7*mesh.Np];
        const dfloat q9  = q[base + 8*mesh.Np];
        const dfloat q10 = q[base + 9*mesh.Np];

        if (q1 <= 1e-14) continue;

        const dfloat p  = q1*c*c;
        const dfloat Cp = (p - pInf)/dynp;

        const dfloat s12 = -c*c*(          q5  - q2*q3/q1);
        const dfloat s13 = -c*c*(          q6  - q2*q4/q1);
        const dfloat s23 = -c*c*(          q7  - q3*q4/q1);
        const dfloat s11 = -c*c*(std::sqrt(2.0)*q8  - q2*q2/q1);
        const dfloat s22 = -c*c*(std::sqrt(2.0)*q9  - q3*q3/q1);
        const dfloat s33 = -c*c*(std::sqrt(2.0)*q10 - q4*q4/q1);

        // viscous traction on body: t_v = -sigma*n
        const dfloat tvx = -(s11*nx + s12*ny + s13*nz);
        const dfloat tvy = -(s12*nx + s22*ny + s23*nz);
        const dfloat tvz = -(s13*nx + s23*ny + s33*nz);

        // tangential wall-shear vector
        const dfloat tv_dot_n = tvx*nx + tvy*ny + tvz*nz;

        const dfloat tauwx = tvx - tv_dot_n*nx;
        const dfloat tauwy = tvy - tv_dot_n*ny;
        const dfloat tauwz = tvz - tv_dot_n*nz;

        const dfloat tauw = std::sqrt(tauwx*tauwx + tauwy*tauwy + tauwz*tauwz);
        const dfloat Cf   = tauw/dynp;

        const dlong nodeId = e*mesh.Np + vid;
        const dfloat x = mesh.x[nodeId];
        const dfloat y = mesh.y[nodeId];
        const dfloat z = mesh.z[nodeId];

        fprintf(fp,
                "%.16e %lld %d %d "
                "%.16e %.16e %.16e "
                "%.16e %.16e %.16e "
                "%.16e %.16e %.16e "
                "%.16e %.16e %.16e "
                "%.16e %.16e %.16e %.16e %.16e\n",
                time,
                (long long)e,
                f,
                n,
                x, y, z,
                nx, ny, nz,
                q1,
                p,
                Cp,
                tvx, tvy, tvz,
                tauwx, tauwy, tauwz,
                tauw,
                Cf);
      }
    }
  }

  fclose(fp);
}