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

void mesh_t::CurvedHex3D(){

  
  //printf("Nbf = %lld\n", NboundaryFaces);

  /* Origin coordinates and radius of the sphere */
  dfloat x0 = 0.0;
  dfloat y0 = 0.0;
  dfloat z0 = 0.0;
  dfloat r_cyl =0.5;
  
  /* Count number of wall boundary faces  */
  hlong bcnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToB[e*Nfaces+f]==1) { // || EToE[e*Nfaces+f]==e)
        //printf("bcnt = %d\n", e);

        dlong id = e*Nverts;

    

        dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = EX[id+1];
        dfloat xe3 = EX[id+2];
        dfloat xe4 = EX[id+3];
        dfloat xe5 = EX[id+4];
        dfloat xe6 = EX[id+5];
        dfloat xe7 = EX[id+6];
        dfloat xe8 = EX[id+7];
    
        dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = EY[id+1];
        dfloat ye3 = EY[id+2];
        dfloat ye4 = EY[id+3];
        dfloat ye5 = EY[id+4];
        dfloat ye6 = EY[id+5];
        dfloat ye7 = EY[id+6];
        dfloat ye8 = EY[id+7];
    
        dfloat ze1 = EZ[id+0]; /* z-coordinates of vertices */
        dfloat ze2 = EZ[id+1];
        dfloat ze3 = EZ[id+2];
        dfloat ze4 = EZ[id+3];
        dfloat ze5 = EZ[id+4];
        dfloat ze6 = EZ[id+5];
        dfloat ze7 = EZ[id+6];
        dfloat ze8 = EZ[id+7];
        

        // printf("nodes of boundary faces xe1= %f, ye1= %f, ze1=%f\n", xe1, ye1,ze1);
        // printf("nodes of boundary faces xe2= %f, ye2= %f, ze2=%f\n", xe2, ye2,ze2);
        // printf("nodes of boundary faces xe3= %f, ye3= %f, ze3=%f\n", xe3, ye3,ze3);
        // printf("nodes of boundary faces xe4= %f, ye4= %f, ze4=%f\n", xe4, ye4,ze4);
      /*  if (e==33621)
       {
         printf("nodes of boundary faces xe1= %f, ye1= %f, ze1=%f\n", xe1, ye1,ze1);
        printf("nodes of boundary faces xe2= %f, ye2= %f, ze2=%f\n", xe2, ye2,ze2);
        printf("nodes of boundary faces xe3= %f, ye3= %f, ze3=%f\n", xe3, ye3,ze3);
        printf("nodes of boundary faces xe4= %f, ye4= %f, ze4=%f\n", xe4, ye4,ze4);
        }*/

        
         dfloat theta1 = atan2(ye1-y0,xe1-x0); 
         dfloat theta2 = atan2(ye2-y0,xe2-x0);
         dfloat theta3 = atan2(ye3-y0,xe3-x0);
         dfloat theta4 = atan2(ye4-y0,xe4-x0);

         dfloat phi1 = atan2(sqrt((ye1-y0)*(ye1-y0)+(xe1-x0)*(xe1-x0)),ze1-z0); 
         dfloat phi2 = atan2(sqrt((ye2-y0)*(ye2-y0)+(xe2-x0)*(xe2-x0)),ze2-z0);
         dfloat phi3 = atan2(sqrt((ye3-y0)*(ye3-y0)+(xe3-x0)*(xe3-x0)),ze3-z0);
         dfloat phi4 = atan2(sqrt((ye4-y0)*(ye4-y0)+(xe4-x0)*(xe4-x0)),ze4-z0);

        dfloat newx1 = x0 + r_cyl*cos(theta1)*sin(phi1); dfloat newy1 = y0 + r_cyl*sin(theta1)*sin(phi1); dfloat newz1 = r_cyl*cos(phi1);
        dfloat newx2 = x0 + r_cyl*cos(theta2)*sin(phi2); dfloat newy2 = y0 + r_cyl*sin(theta2)*sin(phi2); dfloat newz2 = r_cyl*cos(phi2);
        dfloat newx3 = x0 + r_cyl*cos(theta3)*sin(phi3); dfloat newy3 = y0 + r_cyl*sin(theta3)*sin(phi3); dfloat newz3 = r_cyl*cos(phi3);
        dfloat newx4 = x0 + r_cyl*cos(theta4)*sin(phi4); dfloat newy4 = y0 + r_cyl*sin(theta4)*sin(phi4); dfloat newz4 = r_cyl*cos(phi4);
        
        EX[id+0] = newx1; EY[id+0] = newy1; EZ[id+0] = newz1;
        EX[id+1] = newx2; EY[id+1] = newy2; EZ[id+1] = newz2;
        EX[id+2] = newx3; EY[id+2] = newy3; EZ[id+2] = newz3;
        EX[id+3] = newx4; EY[id+3] = newy4; EZ[id+3] = newz4;
        
        
        bcnt = bcnt+1;  

        /*if (e==33621)
        {
         printf("nodes of boundary faces xe1= %f, ye1= %f, ze1=%f, t1=%f\n", xe1, ye1,ze1,theta1);
        printf("nodes of boundary faces xe2= %f, ye2= %f, ze2=%f, t2=%f\n", xe2, ye2,ze2,theta2);
        printf("nodes of boundary faces xe3= %f, ye3= %f, ze3=%f, t3=%f\n", xe3, ye3,ze3,theta3);
        printf("nodes of boundary faces xe4= %f, ye4= %f, ze4=%f, t4=%f\n", xe4, ye4,ze4,theta4);
        }*/
        //printf("nodes of boundary faces x4= %f\n, y4= %f\n", newx4, newy4);
        //LIBP_FORCE_ABORT("Stopped");
      }
    }
  }
        
      

  printf("bcnt = %lld\n", bcnt);
}

} //namespace libp
