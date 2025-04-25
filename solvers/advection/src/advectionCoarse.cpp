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
                         dlong Ncoarse){

  // Store old info & Allocate new arrays
  

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);

  
  // A flag to ommit informations of deleted element
  hlong Delete_Flag[2*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(2*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  #pragma omp parallel for
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    for (int n = 0; n < 3; ++n)
    {
    const dlong id = e*mesh.Nverts+n;
    EToB_new[id] = mesh.EToB[id];
    EToV_new[id] = mesh.EToV[id];
    EX_new[id] = mesh.EX[id]; 
    EY_new[id] = mesh.EY[id]; 
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*2,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming
    dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==-1)
    {
      Coar[ii] = e;
      ii = ii + 1;
    }
  }
    // Determine Triangles To be Coarsened
  hlong del_vertex = 0; // Counts each new_vertex that will be deleted
  #pragma omp parallel for
  for (int i = 0; i < ii; ++i)
  {
      int e = Coar[i];
      const hlong id = e*mesh.Nverts; 
      if (RefFlag[e]==-1 && EToRefLevel[e]>0 /*&& e==PToC[2*e+0]*/)
      {
        //printf("e_conf=%d\n",PToC[2*e+0]*mesh.Np);
        hlong sib_e ;
        hlong sib0_id;
        hlong sib_id;
        hlong v0 ;
        hlong v1 ;
        hlong v2 ;
        hlong v0_sib;
        hlong v1_sib;
        hlong v2_sib;     
  
        // Extract Vertex Number of Elements to be Combined
        sib0_id = PToC[2*e+0]*mesh.Nverts; 
        sib_e = PToC[2*e+1];
        sib_id = sib_e*mesh.Nverts;
        v0 = mesh.EToV[sib0_id+0]; 
        v1 = mesh.EToV[sib0_id+1]; 
        v2 = mesh.EToV[sib0_id+2];
        v0_sib = mesh.EToV[sib_id+0]; 
        v1_sib = mesh.EToV[sib_id+1]; 
        v2_sib = mesh.EToV[sib_id+2];
        

        
        hlong id_del; // Id of the node that will be unnecessary
        if (v0 >= v1 && v0 >= v2) {
        id_del = sib0_id+0;
        } else if (v1 >= v2) {
        id_del = sib0_id+1;
        } else {
        id_del = sib0_id+2;
        }

        hlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
        if (v0_sib >= v1_sib && v0_sib >= v2_sib) {
        id_del2 = sib_id+0;
        } else if (v1_sib >= v2_sib) {
        id_del2 = sib_id+1;
        } else {
        id_del2 = sib_id+2;
        }
        
        if (id_del==sib0_id+1&&id_del2==sib_id+0)
        {
          hlong const neighbor_id = mesh.EToE[sib0_id+0];
           
          if(RefFlag[neighbor_id]!=-1 && mesh.EToE[sib0_id+0]!=-1){
             RefFlag[neighbor_id]=-1;
             Ncoarse++; 
             del_vertex--;
          }
        }

        if (id_del==sib0_id+2&&id_del2==sib_id+1)
        {         
          hlong const neighbor_id = mesh.EToE[sib0_id+1];
          //printf("mesh.EToE[id+1]=%d\n",mesh.EToE[id+1]);
        
          if(RefFlag[neighbor_id]!=-1 && mesh.EToE[id+1]!=-1){

             RefFlag[neighbor_id]=-1;
             Ncoarse++; 
             del_vertex--;
          }
        }

        if (id_del==sib0_id+2&&id_del2==sib_id+0)
        {          
          hlong const neighbor_id = mesh.EToE[sib0_id+2];
     
          if(RefFlag[neighbor_id]!=-1 && mesh.EToE[id+2]!=-1){
             RefFlag[neighbor_id]=-1; 
             Ncoarse++;
             del_vertex--;
          }        
        } 

      }         
  }
  

  ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==-1)
    {
      Coar[ii] = e;
      ii = ii + 1;
    }
  }

  // Coarsement Loop
  // Determine ids of new vertices and EToV
  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 

    if (RefFlag[e]==-1 && EToRefLevel[e]>0)
    {
      
      hlong sib_e=0 ;
      hlong sib_id=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      hlong v0_sib=0;
      hlong v1_sib=0;
      hlong v2_sib=0;
      //printf("e=%d\n",PToC[2*e+0]*mesh.Np);
      
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
          
          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[sib_id+1];
          EToB_new[id+2] = mesh.EToB[id+2];

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;
        }

        if (id_del==id+2&&id_del2==sib_id+1)
        {
          
          EToV_new[id+2] = v2_sib;
          EX_new[id+2] = EX_new[sib_id+2];
          EY_new[id+2] = EY_new[sib_id+2];
       
          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
               
          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[id+1];
          EToB_new[id+2] = mesh.EToB[sib_id+2];

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;  
          
        }

        if (id_del==id+2&&id_del2==sib_id+0)
        {

          EToV_new[id+2] = v2_sib;
          EX_new[id+2] = EX_new[sib_id+2];
          EY_new[id+2] = EY_new[sib_id+2];
    
          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
           
          EToB_new[id+0] = mesh.EToB[id+0];
          EToB_new[id+1] = mesh.EToB[sib_id+1];
          EToB_new[id+2] = mesh.EToB[id+2];

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;         
        } 
    }             
  }   

  // Loop for updating mesh information if coarsening done.

      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              //printf("id=%d\n",e);
            } else{
                    EToRefLevel[e_new] = EToRefLevel[e]; 
                    //IntFlag[e_new] = IntFlag[e];
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


        deviceMemory<dlong> o_combineFlag = platform.malloc<dlong>(CombineFlag);

        deviceMemory<dfloat> o_Q = platform.malloc<dfloat>(Q);
        combineKernel(mesh.Nelements,Ncoarse,o_Q ,o_q, o_combineFlag,o_IntFlag,o_PToC,mesh.o_RM);        
        
        o_q.copyTo(Q);
        
        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
                            if (e==PToC[e*2+0]){
              dlong  sib_e = PToC[e*2+1];
              PToC[e*2+0] = 0; PToC[e*2+1] = 0;
              PToC[sib_e*2+0] = 0; PToC[sib_e*2+1] = 0;
            } else{
              dlong sib_e = PToC[e*2+0];
              PToC[e*2+0] = 0; PToC[e*2+1] = 0;
              PToC[sib_e*2+0] = 0; PToC[sib_e*2+1] = 0;
              //printf("e=%d,sib_e=%d\n",e*mesh.Np,sib_e*mesh.Np);
            }
              //printf("id=%d\n",e);
            } else{ 
                    IntFlag[e_new] = IntFlag[e];
                    if (e==PToC[e*2+1])
                                       {
                    

                    PToC[e_new*2+0] = PToC[e*2+0];
                    PToC[e_new*2+1] = e_new;

                    dlong p_id=PToC[e*2+0];
                    PToC[p_id*2+1] = e_new;

                    //printf("e=%d,e_new=%d,p_id=%d\n",e*mesh.Np,e_new*mesh.Np ,p_id*mesh.Np);
                   // printf("PToC[e_new*2+1]=%d,PToC[p_id*2+1]=%d\n",PToC[p_id*2+1]*mesh.Np,PToC[e_new*2+1]*mesh.Np);
                                       }                   

                   
                    
                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = e_new*mesh.Np+n; 
          
                    Q[id_new] = Q[id];
                    
                    //printf("EToV=%lld\n",mesh.EToV[id]);
                    //printf("EX=%g\n",mesh.EX[id]);
                    //printf("EY=%f\n",mesh.EY[id]);
                    
                    }
                    e_new++;
                  }
        }

        o_PToC = platform.malloc<dlong>(PToC);
        
        for (int e = e_new; e < mesh.Nelements+nn; ++e)
        {
          EToRefLevel[e] = 0; 
          for (int n = 0; n < 2; ++n)
          {
          dlong id = e*2+n;

          PToC[id] = 0;
    
          }         
        }

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%d\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        
        deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        o_oldq = platform.malloc<dfloat>(Q);
        o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
      }         
}
