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

namespace libp {

void mesh_t::CurvedQuad2D(){

  //printf("HELLO\n");
  //printf("Nbf = %lld\n", NboundaryFaces);

  /* Origin coordinates and radius of the cylinder */
  dfloat x0 = 0.0;
  dfloat y0 = 0.0;
  dfloat r_cyl =1.0;
  
    /* count number of boundary faces (i.e. not yet connected) */
  hlong bcnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToB[e*Nfaces+f]==1) { // || EToE[e*Nfaces+f]==e)
        //printf("bcnt = %d\n", e);

        dlong id = e*Nverts;

        dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
        dfloat xe4 = EX[id+3];

        dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
        dfloat ye4 = EY[id+3];
        
         //printf("nodes of boundary faces x4= %f\n, y4= %f\n", xe4, ye4);

        dfloat theta1 = atan2(ye1-y0,xe1-x0); 
        dfloat theta2 = atan2(ye4-y0,xe4-x0);

        dfloat newx1 = x0 + r_cyl*cos(theta1); dfloat newy1 = y0 + r_cyl*sin(theta1);
        dfloat newx4 = x0 + r_cyl*cos(theta2); dfloat newy4 = y0 + r_cyl*sin(theta2);
        
        EX[id+0] = newx1;
        EX[id+3] = newx4;
        EY[id+0] = newy1;
        EY[id+3] = newy4;
        bcnt = bcnt+1;  
        //printf("nodes of boundary faces x4= %f\n, y4= %f\n", newx4, newy4);
      }
    }
  }
        
      

  //printf("bcnt = %lld\n", bcnt);
}

} //namespace libp
