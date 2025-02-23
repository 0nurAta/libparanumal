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

// 
void advection_t::Coarse(memory<dfloat>& Q,
                         memory<dlong>& RefFlag,
                         dlong Nrefine){

  

  // extract q halo on DEVICE
  //traceHalo.ExchangeStart(o_Q, 1);

  // Store old info & Allocate new arrays
  
  // Solution info
  //memory <dfloat> qold = Q;
  
  // Element to vertex & Element to boundary connectivity 
  //memory<hlong>EToV_old(mesh.Nelements*mesh.Nverts);
  memory<hlong>EToV_new(mesh.Nelements*mesh.Nverts+2*Nrefine*mesh.Nverts);
  memory<int>EToB_new(mesh.Nelements*mesh.Nverts+2*Nrefine*mesh.Nverts);

  // Vertex physical coordinates
  //memory<dfloat>EX_old(mesh.Nelements*mesh.Nverts);
  //memory<dfloat>EY_old(mesh.Nelements*mesh.Nverts);

  memory<dfloat>EX_new(mesh.Nelements*mesh.Nverts+2*Nrefine*mesh.Nverts);
  memory<dfloat>EY_new(mesh.Nelements*mesh.Nverts+2*Nrefine*mesh.Nverts);
   //printf("EToVsize=%ld\n",EToV_new.size());
  
  // A flag to ommit informations of deleted element
  hlong Delete_Flag[mesh.Nelements]={};
  // Copy old Element to Vertex Connectivity to the New One
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    for (int n = 0; n < 3; ++n)
    {
    const dlong id = e*mesh.Nverts+n;
    EToB_new[id] = mesh.EToB[id];
    EToV_new[id] = mesh.EToV[id];
    EX_new[id] = mesh.EX[id]; 
    EY_new[id] = mesh.EY[id]; 
    //printf("EToE=%lld\n",mesh.EToE[id]);
    //printf("EX=%g\n",mesh.EX[id]);
    //printf("EY=%f\n",mesh.EY[id]);
    }
  }
  
    // Determine Triangles To be Coarsened
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      
      const dlong id = e*mesh.Nverts; 
      if (RefFlag[e]==-1 && EToRefLevel[e]!=0 && e==PToC[2*e+0])
      {

        hlong sib_e ;
        hlong sib_id;
        hlong v0 ;
        hlong v1 ;
        hlong v2 ;
        hlong v0_sib;
        hlong v1_sib;
        hlong v2_sib;     
           
        // Extract Vertex Number of Elements to be Combined
        sib_e = PToC[2*e+1];
        sib_id = sib_e*mesh.Nverts;
        v0 = mesh.EToV[id+0]; 
        v1 = mesh.EToV[id+1]; 
        v2 = mesh.EToV[id+2];
        v0_sib = mesh.EToV[sib_id+0]; 
        v1_sib = mesh.EToV[sib_id+1]; 
        v2_sib = mesh.EToV[sib_id+2];
        

        
         dlong id_del; // Id of the node that will be unnecessary
         if (v0 >= v1 && v0 >= v2) {
         id_del = id+0;
         } else if (v1 >= v2) {
         id_del = id+1;
         } else {
         id_del = id+2;
         }

         dlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
         if (v0_sib >= v1_sib && v0_sib >= v2_sib) {
         id_del2 = sib_id+0;
         } else if (v1_sib >= v2_sib) {
         id_del2 = sib_id+1;
         } else {
         id_del2 = sib_id+2;
         }
 
        if (id_del==id+1&&id_del2==sib_id+0)
        {
          hlong const neighbor_id = mesh.EToE[id+0];

          if(RefFlag[neighbor_id]==0 && mesh.EToE[id+1]!=-1){
             RefFlag[neighbor_id]=-1;
             Ncoarse++; 
             printf("hello1\n");
          }
        }

        if (id_del==id+2&&id_del2==sib_id+1)
        {         
          hlong const neighbor_id = mesh.EToE[id+1];

          if(RefFlag[neighbor_id]==0 && mesh.EToE[id+1]!=-1){
             RefFlag[neighbor_id]=-1;
             Ncoarse++; 
             printf("neighbor_id=%d\n",sib_e);
             printf("neighbor_id=%lld\n",neighbor_id);
             printf("Level=%lld\n",EToRefLevel[neighbor_id]);
           
             printf("neighbor_id%lld\n",e);
          }
        }

        if (id_del==id+2&&id_del2==sib_id+2)
        {          
          hlong const neighbor_id = mesh.EToE[id+2];

          if(RefFlag[neighbor_id]==0 && mesh.EToE[id+1]!=-1){
             RefFlag[neighbor_id]=-1; 
             Ncoarse++;
             //printf("hello3\n");
          }        
        } 

      }             
  }
  // Coarsement Loop
  // Determine ids of new vertices and EToV
  
  for (int e = 0; e < mesh.Nelements; ++e)
  {
      
    dlong id = e*mesh.Nverts; 

    if (RefFlag[e]==-1 && EToRefLevel[e]!=0)
    {
      
      hlong sib_e=0 ;
      hlong sib_id=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      hlong v0_sib=0;
      hlong v1_sib=0;
      hlong v2_sib=0;
      
      if (e==PToC[2*e+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[2*e+1];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
      }

      if (e==PToC[2*e+1]) // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[2*e+0];
         sib_id = sib_e*mesh.Nverts;
         // swap id's 
         hlong id_swap = id;
         id = sib_id;
         sib_id = id_swap;
         // swap id's 
         hlong e_swap = e;
         e = sib_e;
         sib_e = e_swap;

         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
      }

         dlong id_del; // Id of the node that will be unnecessary
         if (v0 >= v1 && v0 >= v2) {
         id_del = id+0;
         } else if (v1 >= v2) {
         id_del = id+1;
         } else {
         id_del = id+2;
         }

         dlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
         if (v0_sib >= v1_sib && v0_sib >= v2_sib) {
         id_del2 = sib_id+0;
         } else if (v1_sib >= v2_sib) {
         id_del2 = sib_id+1;
         } else {
         id_del2 = sib_id+2;
         }
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
        if (id_del==id+1&&id_del2==sib_id+0)
        {

          EToV_new[id+1] = v1_sib;
          EX_new[id+1] = EX_new[sib_id+1];
          EY_new[id+1] = EY_new[sib_id+1];
          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          
          //PToC[e*2]   = e;
          //PToC[e*2+1] = id_new;
        }

        if (id_del==id+2&&id_del2==sib_id+1)
        {
          
          EToV_new[id+2] = v2_sib;
          EX_new[id+2] = EX_new[sib_id+2];
          EY_new[id+2] = EY_new[sib_id+2];
          //printf("e=%d\n",e);
          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
            
          //PToC[e*2]   = e;
          //PToC[e*2+1] = id_new;
        }

        if (id_del==id+2&&id_del2==sib_id+2)
        {

          EToV_new[id+2] = v1_sib;
          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];
          //printf("e=%d\n",e);
          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;         
        }

        
    }             
  }   

  // Coarsening Loop

      if (Ncoarse!=0)
      { 
        printf("inside coarsement\n"); 
        
        //mesh.EToV.free();
        //mesh.EToB.free();
        //mesh.EX.free();
        //mesh.EY.free();
        
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              //printf("id=%d\n",e);
            } else{
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = e_new*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id]; 
                    //printf("EToV=%lld\n",mesh.EToV[id]);
                    //printf("EX=%g\n",mesh.EX[id]);
                    //printf("EY=%f\n",mesh.EY[id]);
                    
                    }
                    e_new++;
                  }
        }        
        mesh.Nelements = mesh.Nelements - Ncoarse;
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
        //mesh.SetupUpdate(Ncoarse);
        mesh_t updatedMeshPtr = mesh.SetupUpdate(Ncoarse);
        mesh = updatedMeshPtr;
      }



         
}
