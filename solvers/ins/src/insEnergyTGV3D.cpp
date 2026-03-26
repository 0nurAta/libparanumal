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

// Calculate the kinetic energy for TGV test cases in HEXes
void ins_t::EnergyTGV(memory<dfloat>& U, std::string fileName, dfloat time){

	FILE *fp;
	dfloat KE = 0.0;

	if(mesh.rank==0)
		fp = fopen(fileName.c_str(), "a+");

	memory<dfloat> u_tgv(mesh.Np);
	memory<dfloat> v_tgv(mesh.Np);
	memory<dfloat> w_tgv(mesh.Np);

	for(dlong e=0;e<mesh.Nelements;++e){
		for(int n=0;n<mesh.Np;++n){

			int id = n + e*mesh.Np*NVfields;

			u_tgv[n] = U[id + 0*mesh.Np];
			v_tgv[n] = U[id + 1*mesh.Np];
			w_tgv[n] = U[id + 2*mesh.Np];

			KE += mesh.wJ[n+e*mesh.Np]*(u_tgv[n]*u_tgv[n]+v_tgv[n]*v_tgv[n]+w_tgv[n]*w_tgv[n]);
		}
	}

	comm.Allreduce(KE, Comm::Sum);
	const dfloat volScale = 0.5/(6.28318530718*6.28318530718*6.28318530718); 

	KE *= volScale;

	if(mesh.rank==0){
		fprintf(fp, "%5.4f %13.12f\n", time, KE);
		fclose(fp);
	}
}
