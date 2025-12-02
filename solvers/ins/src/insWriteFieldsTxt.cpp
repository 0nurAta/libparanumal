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

#include "ins.hpp"

// interpolate data to plot nodes and save to file (one per process)
void ins_t::WriteFieldsTxt(memory<dfloat>& U, memory<dfloat>& P, std::string fileName, dfloat time){

	FILE *fp;


	  fp = fopen(fileName.c_str(), "w");

	for(dlong e=0;e<mesh.Nelements;++e){
		for(int n=0;n<mesh.Np;++n){
			const dfloat x = mesh.x[e*mesh.Np+n];
			const dfloat y = mesh.y[e*mesh.Np+n];
			const dfloat p = P[e*mesh.Np+n];
			const dfloat u  = U[e*mesh.Np*NVfields+n+mesh.Np*0];
			const dfloat v  = U[e*mesh.Np*NVfields+n+mesh.Np*1];
			if(mesh.dim==3){
				const dfloat z = mesh.z[e*mesh.Np+n];
				const dfloat w  = U[e*mesh.Np*NVfields+n+mesh.Np*1];
				fprintf(fp, "%g,%g,%g,%g,%g,%g,%g\n", x, y, z, p, u, v, w);
			}
			else
				fprintf(fp, "%g,%g, %g,%g,%g\n", x, y, p, u, v);
			
		}
	}

	fclose(fp);

}