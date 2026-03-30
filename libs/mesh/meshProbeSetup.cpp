#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "mesh.hpp"

namespace libp {

void mesh_t::ProbeSetup() {
  probesEnabled = settings.compareSetting("PROBES", "TRUE");
  if (!probesEnabled) return;

  std::string probeFile;

  settings.getSetting("PROBE FILE", probeFile);
  settings.getSetting("PROBE OUTPUT FILE", probeOutputFile);
  settings.getSetting("PROBE OUTPUT INTERVAL", probeOutputInterval);
  Nprobes = 0;
  if (dim == 2) {
    ReadProbeFile2D(probeFile, probeX, probeY);
    ProbeQuad2D(); // replace with ProbeTri2D() later if needed
    Nprobes = probeElement.length();
  } else if (dim == 3) {
    ReadProbeFile3D(probeFile, probeX, probeY, probeZ);
    ProbeHex3D();
    Nprobes = probeElement.length();
  } else {
    LIBP_ABORT("ProbeSetup: unsupported mesh dimension", 1);
  }

  if (rank == 0) {
    std::cout << "Loaded " << Nprobes << " probes\n";
  }
}

void mesh_t::ReadProbeFile2D(const std::string& fileName,
                             memory<dfloat>& pX,
                             memory<dfloat>& pY) {

  std::ifstream fp(fileName);
  if (!fp.is_open()) {
    printf("ReadProbeFile2D: could not open probe file: %s\n", fileName.c_str());
    LIBP_ABORT("ReadProbeFile2D: could not open probe file", 1);
  }

  std::vector<dfloat> xlist, ylist;
  std::string line;

  while (std::getline(fp, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    dfloat x_p, y_p;

    if (!(iss >> x_p >> y_p)) continue;

    xlist.push_back(x_p);
    ylist.push_back(y_p);
  }

  const dlong N_probe = static_cast<dlong>(xlist.size());

  pX.malloc(N_probe);
  pY.malloc(N_probe);

  for (dlong n = 0; n < N_probe; ++n) {
    pX[n] = xlist[n];
    pY[n] = ylist[n];
  }
}

void mesh_t::ReadProbeFile3D(const std::string& fileName,
                             memory<dfloat>& pX,
                             memory<dfloat>& pY,
                             memory<dfloat>& pZ) {

  std::ifstream fp(fileName);
  if (!fp.is_open()) {
    printf("ReadProbeFile3D: could not open probe file: %s\n", fileName.c_str());
    LIBP_ABORT("ReadProbeFile3D: could not open probe file", 1);
  }

  std::vector<dfloat> xlist, ylist, zlist;
  std::string line;

  while (std::getline(fp, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    dfloat x_p, y_p, z_p;

    if (!(iss >> x_p >> y_p >> z_p)) continue;

    xlist.push_back(x_p);
    ylist.push_back(y_p);
    zlist.push_back(z_p);
  }

  const dlong N_probe = static_cast<dlong>(xlist.size());

  pX.malloc(N_probe);
  pY.malloc(N_probe);
  pZ.malloc(N_probe);

  for (dlong n = 0; n < N_probe; ++n) {
    pX[n] = xlist[n];
    pY[n] = ylist[n];
    pZ[n] = zlist[n];
  }
}

} // namespace libp