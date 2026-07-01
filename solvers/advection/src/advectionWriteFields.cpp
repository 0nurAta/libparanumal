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

#include "advection.hpp"

// interpolate data to plot nodes and save to file (one per process)
void advection_t::PlotAdaptivity(std::string fileName, dfloat time, dfloat norm){

  FILE *fp;


  if(mesh.rank==0)
    fp = fopen(fileName.c_str(), "a+");

  //fprintf(fp, "\n Kinetic Energy>\n");
    const dlong TotalElem = mesh.Nelements;
    memory<dfloat>  Reg_Array(mesh.Nelements,0);
    dfloat WorstReg = 0.0;
    // write out velocity
    
    WorstReg = 0.0;
    
    for(dlong e = 0; e < mesh.Nelements; ++e){
    
        const dlong id = e*mesh.Nfaces;
    
        const dfloat x0 = mesh.EX[id+0], y0 = mesh.EY[id+0];
        const dfloat x1 = mesh.EX[id+1], y1 = mesh.EY[id+1];
        const dfloat x2 = mesh.EX[id+2], y2 = mesh.EY[id+2];
    
        // Edge lengths
        const dfloat L0 = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)); // edge (0,1)
        const dfloat L1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); // edge (1,2)
        const dfloat L2 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0)); // edge (2,0)
    
        // Triangle area using cross product
        const dfloat Area =
            0.5*fabs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
    
        // Diameter = longest edge
        dfloat Diameter = L0;
    
        if(L1 > Diameter){
            Diameter = L1;
        }
      
        if(L2 > Diameter){
            Diameter = L2;
        }
      
        // Shape regularity factor
        const dfloat Reg = Diameter / sqrt(Area);
      
        Reg_Array[e] = Reg;
      
        if(Reg > WorstReg){
            WorstReg = Reg;
        }
    }

    dfloat Reg_Sum = 0.0 ;
    // Averaging Regularity
    for(dlong e=0;e<mesh.Nelements;++e){
    
   	 Reg_Sum += Reg_Array[e];
    }

    dfloat AvgReg = Reg_Sum / TotalElem; 

if(mesh.rank==0){
  fprintf(fp,"%5.4f %5.8f %d %13.12f %13.12f\n",time, norm, TotalElem, AvgReg, WorstReg);
  fclose(fp);
}

}
