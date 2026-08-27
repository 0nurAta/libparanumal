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

#include "adaptivity.hpp"
#include "mesh.hpp"

namespace libp {

void adaptivity_t::Coarse(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(2*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(2*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  
  // Permutation map
  dlong perm[2*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[2*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(2*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < level+3; ++i)
    {
      const dlong id = e*(level+3)+i;
      PToC_new[id] = PToC[id];
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
  //#pragma omp parallel for
  for (int i = 0; i < ii; ++i)
  {

      int e = Coar[i];
      const hlong id = e*mesh.Nverts; 
       //printf("I am Here!!!!\n");  
      //printf("refflag=%d\n",RefFlag[e]); 
      if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1 && EToRefLevel[e]>0)
      {

        
        //printf("e_coarse=%d,e=%d,level[%d]=%d\n",PToC[4*e+0],e,PToC[4*e+0],EToRefLevel[e]);
        //printf("x1=%f,x2=%f,x3=%f\n",EX_new[id+0],EX_new[id+1],EX_new[id+2]);
        //printf("y1=%f,y2=%f,y3=%f\n",EY_new[id+0],EY_new[id+1],EY_new[id+2]);
        hlong sib_e ;
        hlong e_parent;
        hlong id;
        hlong sib_id;
        hlong v0 ;
        hlong v1 ;
        hlong v2 ;
        hlong v0_sib;
        hlong v1_sib;
        hlong v2_sib;    
        dlong swap = 0; 
        
        if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         id = PToC[e*(level+3)+0]*mesh.Nverts; 
         sib_e = PToC[e*(level+3)+EToRefLevel[e]];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         

      }

      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         id = PToC[e*(level+3)+0]*mesh.Nverts; 
         sib_e = e;
         e = PToC[e*(level+3)+0];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         swap = 1;


      }
                      if (EToRefLevel[e]==EToRefLevel[sib_e] )
        {
       

        hlong id_del; // Id of the node that will be unnecessary
        id_del = id+2;
        printf("RefFlag[%d]=%d\n",e,RefFlag[e] );
        hlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
        id_del2 = sib_id+2;
        // For Rule 0 & 1
        if (mesh.EToE[id+2]==sib_e)
        {
           
          hlong const neighbor_id = mesh.EToE[id+1];
          hlong const neighbor_sibid = mesh.EToE[sib_id+2];
          printf("e_coarse01=%d,e_coarse01=%d,e_coarse01=%d,EToRefLevel[neighbor_id]=%d,EToRefLevel[neighbor_sibid]=%d\n",e,neighbor_id,neighbor_sibid,PToC[neighbor_id*(level+3)+0],PToC[neighbor_sibid*(level+3)+1]);
          printf("e_coarse01=%d,e_coarse01=%d,e_coarse01=%d,True=%d,True=%d\n",e,neighbor_id,
          neighbor_sibid,EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid] ,(PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id));

          if (EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid] && (PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id))
          {
          
          if(mesh.EToE[id+1]!=-1){
            
             RefFlag[PToC[neighbor_id*(level+3)+0]]=-1;
             RefFlag[PToC[neighbor_sibid*(level+3)+0]]=-1;
             Ncoarse++; 
             del_vertex--;
          } 
          } else if(mesh.EToE[id+1]==-1){}
          else{printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}

        }

        // For Rule 2
        else if (mesh.EToE[id+1]==sib_e)
        {
          hlong const neighbor_id = mesh.EToE[id+2];
          hlong const neighbor_sibid = mesh.EToE[sib_id+1];
          printf("e_coarse02=%d,e_coarse02=%d,e_coarse02=%d\n",e,neighbor_id,neighbor_sibid);
           if (EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid]&& (PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id))
          {
            
          if(mesh.EToE[id+2]!=-1){
             RefFlag[PToC[neighbor_id*(level+3)+0]]=-1;
             RefFlag[PToC[neighbor_sibid*(level+3)+0]]=-1;
             Ncoarse++; 
             del_vertex--;
           }
          }  else if(mesh.EToE[id+1]==-1){}
          else {printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}
        }

      }  else {printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}
    } else {RefFlag[e]=0;}
       
  }
  
  printf("RefFlag[93]=%d,EToRefLevel[93]=%d,sib93=%d,sib93=%d,,sib93=%d\n",RefFlag[93],EToRefLevel[93],mesh.EToE[93*mesh.Nverts+0],mesh.EToE[93*mesh.Nverts+1],mesh.EToE[93*mesh.Nverts+2]);
  ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==-1)
    {
      Coar[ii] = e;
      ii = ii + 1;
    }
  }
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

        const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(6*Np*Np,0);
      o_RM.copyTo(RM);


  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 

    if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1 && EToRefLevel[e]>0)
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
      
      if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+EToRefLevel[e]];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+0];
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

                if (EToRefLevel[e]==EToRefLevel[sib_e])
        {

          // Rule 0, IntFlag[e]==1
          if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          {    
          EToV_new[id+0] = v1;
          EToV_new[id+1] = v0_sib;  
          EToV_new[id+2] = v1_sib;

          EX_new[id+0] = EX_new[id+1];
          EX_new[id+1] = EX_new[sib_id+0];
          EX_new[id+2] = EX_new[sib_id+1];
          
          EY_new[id+0] = EY_new[id+1];
          EY_new[id+1] = EY_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id+2];
          EToB_new[id+1] = mesh.EToB[sib_id+0];
          EToB_new[id+2] = mesh.EToB[id+0];

          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 0;
          const dlong id_int_c = 1;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
          printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          }

          // Rule 1, IntFlag[e]==5
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==5)
          {    
          EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
                    EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                          for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 4;
          const dlong id_int_c = 5;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          }

          // Rule 2, IntFlag[e]==3
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==3)
          {    
          EToV_new[id+2] = v1_sib;

          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

           for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 2;
          const dlong id_int_c = 3;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);
          } 
        }
    }             
  }   

  // Loop for updating mesh information if coarsening done.
                printf("PToC[52]_Before=%d,%d,%d,%d,%d\n",PToC[52*(level+3)+0],PToC[52*(level+3)+1],PToC[52*(level+3)+2],PToC[52*(level+3)+3],PToC[52*(level+3)+4] );
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        // Map ids of elements to new ones
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
          if (Delete_Flag[e]==1)
          {
            perm[e] = -1;
          }else{
            perm[e] = e_new;
            e_new++;
          }
        }

        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              //printf("id=%d\n",e);
            } else{
                   // EToRefLevel[e_new] = EToRefLevel[e]; 
                    //IntFlag[e_new] = IntFlag[e];
                    if (perm[e]>= 0)
                   {
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = perm[e]*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id];
                    }

                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = perm[e]*mesh.Np+n; 
                    
                    Q_new[id_new] = Q[id];
                    
                    //printf("EToV=%lld\n",mesh.EToV[id]);
                    //printf("EX=%g\n",mesh.EX[id]);
                    //printf("EY=%f\n",mesh.EY[id]);
                    
                    }
                    
                    for (int i = 0; i < level+3; ++i) {
                    dlong oldVal = PToC[e*(level+3)+i];
                    dlong newVal = 0;
                  
                    if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
                      newVal = perm[oldVal];     // valid surviving element
                    else
                      newVal = -1;                // no link / deleted
                  
                    PToC_new[perm[e]*(level+3)+i] = newVal;
                    IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
                  }

                     EToRefLevel_new[perm[e]] = EToRefLevel[e];
                     RefFlag_new[perm[e]] = RefFlag[e];


                    e_new++;  

                    }
                    }
 
        }


        for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
        {
          EToRefLevel_new[e] = 0; 
          RefFlag_new[e] = 0;

          for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    Q_new[id] = 0.f;
                    }

          for (int n = 0; n < (level+3); ++n)
          {
          dlong id = e*(level+3)+n;

          PToC_new[id] = -1;
          IntFlag_new[id] = 0;
          }         
        }
        RefFlag= RefFlag_new;
        PToC = PToC_new;
        EToRefLevel = EToRefLevel_new;
        IntFlag = IntFlag_new;
        Q_new.copyTo(Q); 
        o_PToC = platform.malloc<dlong>(PToC);

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        //mesh.PmlSetup();
        //deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        //o_oldq = platform.malloc<dfloat>(Q);
        //o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
        o_q.copyFrom(Q);
      }         
}

void adaptivity_t::CoarsebyID(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(16*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(16*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(16*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(16*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(16*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(16*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(16*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(16*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  
  // Permutation map
  dlong perm[16*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[16*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(16*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < level+3; ++i)
    {
      const dlong id = e*(level+3)+i;
      PToC_new[id] = PToC[id];
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*64,0); // Array holds element ids for refining. Holds some extra mem. for 
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
  //#pragma omp parallel for
  for (int i = 0; i < ii; ++i)
  {

      int e = Coar[i];
      const hlong id = e*mesh.Nverts; 
       //printf("I am Here!!!!\n");  
      //printf("refflag=%d\n",RefFlag[e]); 
      if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1 && EToRefLevel[e]>0)
      {

        
        //printf("e_coarse=%d,e=%d,level[%d]=%d\n",PToC[4*e+0],e,PToC[4*e+0],EToRefLevel[e]);
        //printf("x1=%f,x2=%f,x3=%f\n",EX_new[id+0],EX_new[id+1],EX_new[id+2]);
        //printf("y1=%f,y2=%f,y3=%f\n",EY_new[id+0],EY_new[id+1],EY_new[id+2]);
        hlong sib_e ;
        hlong e_parent;
        hlong id;
        hlong sib_id;
        hlong v0 ;
        hlong v1 ;
        hlong v2 ;
        hlong v0_sib;
        hlong v1_sib;
        hlong v2_sib;    
        dlong swap = 0; 
        
        if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         id = PToC[e*(level+3)+0]*mesh.Nverts; 
         sib_e = PToC[e*(level+3)+EToRefLevel[e]];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         

      }

      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         id = PToC[e*(level+3)+0]*mesh.Nverts; 
         sib_e = e;
         e = PToC[e*(level+3)+0];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         swap = 1;


      }
                      if (EToRefLevel[e]==EToRefLevel[sib_e] )
        {
       

        hlong id_del; // Id of the node that will be unnecessary
        id_del = id+2;
        printf("RefFlag[%d]=%d\n",e,RefFlag[e] );
        hlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
        id_del2 = sib_id+2;
        // For Rule 0 & 1
        if (mesh.EToE[id+2]==sib_e)
        {
           
          hlong const neighbor_id = mesh.EToE[id+1];
          hlong const neighbor_sibid = mesh.EToE[sib_id+2];
          printf("e_coarse01=%d,e_coarse01=%d,e_coarse01=%d,EToRefLevel[neighbor_id]=%d,EToRefLevel[neighbor_sibid]=%d\n",e,neighbor_id,neighbor_sibid,PToC[neighbor_id*(level+3)+0],PToC[neighbor_sibid*(level+3)+1]);
          printf("e_coarse01=%d,e_coarse01=%d,e_coarse01=%d,True=%d,True=%d\n",e,neighbor_id,
          neighbor_sibid,EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid] ,(PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id));

          if (EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid] && (PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id))
          {
          
          if(mesh.EToE[id+1]!=-1){
            
             RefFlag[PToC[neighbor_id*(level+3)+0]]=-1;
             RefFlag[PToC[neighbor_sibid*(level+3)+0]]=-1;
             Ncoarse++; 
             del_vertex--;
          } 
          } else if(mesh.EToE[id+1]==-1){}
          else{printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}

        }

        // For Rule 2
        else if (mesh.EToE[id+1]==sib_e)
        {
          hlong const neighbor_id = mesh.EToE[id+2];
          hlong const neighbor_sibid = mesh.EToE[sib_id+1];
          printf("e_coarse02=%d,e_coarse02=%d,e_coarse02=%d\n",e,neighbor_id,neighbor_sibid);
           if (EToRefLevel[neighbor_id]==EToRefLevel[neighbor_sibid]&& (PToC[neighbor_id*(level+3)+EToRefLevel[neighbor_id]]==neighbor_sibid || PToC[neighbor_sibid*(level+3)+EToRefLevel[neighbor_sibid]]==neighbor_id))
          {
            
          if(mesh.EToE[id+2]!=-1){
             RefFlag[PToC[neighbor_id*(level+3)+0]]=-1;
             RefFlag[PToC[neighbor_sibid*(level+3)+0]]=-1;
             Ncoarse++; 
             del_vertex--;
           }
          }  else if(mesh.EToE[id+1]==-1){}
          else {printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}
        }

      }  else {printf("unflagged =%d\n",e);(swap==0)? RefFlag[e]=0:RefFlag[sib_e]=0;}
    } else {RefFlag[e]=0;}
       
  }
  
  printf("RefFlag[93]=%d,EToRefLevel[93]=%d,sib93=%d,sib93=%d,,sib93=%d\n",RefFlag[93],EToRefLevel[93],mesh.EToE[93*mesh.Nverts+0],mesh.EToE[93*mesh.Nverts+1],mesh.EToE[93*mesh.Nverts+2]);
  ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==-1)
    {
      Coar[ii] = e;
      ii = ii + 1;
    }
  }
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

        const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(6*Np*Np,0);
      o_RM.copyTo(RM);


  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 

    if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1 && EToRefLevel[e]>0)
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
      
      if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+EToRefLevel[e]];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+0];
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

         //dlong id_del; // Id of the node that will be unnecessary
         //id_del = id+2;
      

         //dlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
         //id_del2 = sib_id+2;
         
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e])
        {

           // Rule 01, IntFlag[e]==5
      /*     if (mesh.EToE[id+2]==sib_e)
        {
            EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          
          for(int n=0;n<mesh.Np;++n){
          const dlong id_p = e*mesh.Np;
          Q[id_p+n] = 0.0;
          }                  

          for(int n=0;n<mesh.Np;++n){
            dfloat qn1=0.;
            dfloat qn2=0.;
            const dlong id_p = e*mesh.Np;
            const dlong id_c = sib_e*mesh.Np;
            const dlong id_int_p = 4;
            const dlong id_int_c = 5;
              for(int ip=0;ip<mesh.Np;++ip){
              const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
      
                   qn1 += Ri*Qold[id_p+ip];       
              }
             for(int ip=0;ip<mesh.Np;++ip){
              const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                   qn2 += Ri*Qold[id_c+ip];  
              }        
            Q[id_p+n] = qn1+qn2;
            }     

            for(int n=0;n<mesh.Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }

          if (PToC[e*(level+3)+EToRefLevel[e]]==mesh.EToE[id+0])
          {
            EToV_new[id+2] = EToV_new[id+0];
            EToV_new[id+0] = EToV_new[id+1];
            EToV_new[id+1] = v0_sib;

            EX_new[id+2] = EX_new[id+0];
            EX_new[id+0] = EX_new[id+1];
            EX_new[id+1] = EX_new[sib_id+0];           
            

            EY_new[id+2] = EY_new[id+0];
            EY_new[id+0] = EY_new[id+1]; 
            EY_new[id+1] = EY_new[sib_id+0];

            EToB_new[id+2] = mesh.EToB[id+0];      
            EToB_new[id+0] = mesh.EToB[id+1];
            EToB_new[id+1] = mesh.EToB[id+2];
            
            printf("INSIDE CORRECT, e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
             printf("RefFlag[%d]=%d,EToRefLevel[53]=%d,sib53=%d,sib53=%d,,sib53=%d\n",e,RefFlag[53],EToRefLevel[53],mesh.EToE[53*mesh.Nverts+0],mesh.EToE[53*mesh.Nverts+1],mesh.EToE[53*mesh.Nverts+2]);
          }

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          IntFlag[e]=5;
          IntFlag[sib_e]=6;
          del_vertex++;
          nn++;
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          CombineFlag[e]=0;
          CombineFlag[sib_e]=0;
        }  

        // Rule 2, IntFlag[e]==3
           if (mesh.EToE[id+1]==sib_e)
        {
          EToV_new[id+2] = v1_sib;
          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          IntFlag[e]=3;
          IntFlag[sib_e]=4;
          del_vertex++;
          nn++;
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);



      for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

    for(int n=0;n<mesh.Np;++n){
      dfloat qn1=0.;
      dfloat qn2=0.;
      const dlong id_p = e*mesh.Np;
      const dlong id_c = sib_e*mesh.Np;
      const dlong id_int_p = 2;
      const dlong id_int_c = 3;
        for(int ip=0;ip<mesh.Np;++ip){
        const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];

             qn1 += Ri*Qold[id_p+ip];       
        }
       for(int ip=0;ip<mesh.Np;++ip){
        const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
             qn2 += Ri*Qold[id_c+ip];  
        }        
      Q[id_p+n] = qn1+qn2;
      }                  
          CombineFlag[e]=0;
          CombineFlag[sib_e]=0;

        } */
          // Rule 0, IntFlag[e]==1
          if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          {    
          EToV_new[id+0] = v1;
          EToV_new[id+1] = v0_sib;  
          EToV_new[id+2] = v1_sib;

          EX_new[id+0] = EX_new[id+1];
          EX_new[id+1] = EX_new[sib_id+0];
          EX_new[id+2] = EX_new[sib_id+1];
          
          EY_new[id+0] = EY_new[id+1];
          EY_new[id+1] = EY_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id+2];
          EToB_new[id+1] = mesh.EToB[sib_id+0];
          EToB_new[id+2] = mesh.EToB[id+0];

          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 0;
          const dlong id_int_c = 1;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
          printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          }

          // Rule 1, IntFlag[e]==5
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==5)
          {    
          EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
                    EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                          for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 4;
          const dlong id_int_c = 5;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          }

          // Rule 2, IntFlag[e]==3
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==3)
          {    
          EToV_new[id+2] = v1_sib;

          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          RefFlag[sib_e]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

           for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 2;
          const dlong id_int_c = 3;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);
          } 
        }
    }             
  }   

  // Loop for updating mesh information if coarsening done.
                printf("PToC[52]_Before=%d,%d,%d,%d,%d\n",PToC[52*(level+3)+0],PToC[52*(level+3)+1],PToC[52*(level+3)+2],PToC[52*(level+3)+3],PToC[52*(level+3)+4] );
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        // Map ids of elements to new ones
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
          if (Delete_Flag[e]==1)
          {
            perm[e] = -1;
          }else{
            perm[e] = e_new;
            e_new++;
          }
        }

        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              //printf("id=%d\n",e);
            } else{
                   // EToRefLevel[e_new] = EToRefLevel[e]; 
                    //IntFlag[e_new] = IntFlag[e];
                    if (perm[e]>= 0)
                   {
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = perm[e]*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id];
                    }

                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = perm[e]*mesh.Np+n; 
                    
                    Q_new[id_new] = Q[id];
                    
                    //printf("EToV=%lld\n",mesh.EToV[id]);
                    //printf("EX=%g\n",mesh.EX[id]);
                    //printf("EY=%f\n",mesh.EY[id]);
                    
                    }
                    
                    for (int i = 0; i < level+3; ++i) {
                    dlong oldVal = PToC[e*(level+3)+i];
                    dlong newVal = 0;
                  
                    if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
                      newVal = perm[oldVal];     // valid surviving element
                    else
                      newVal = -1;                // no link / deleted
                  
                    PToC_new[perm[e]*(level+3)+i] = newVal;
                    IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
                  }

                     EToRefLevel_new[perm[e]] = EToRefLevel[e];
                     RefFlag_new[perm[e]] = RefFlag[e];


                    e_new++;  

                    }
                    }
 
            //        e_new++;
              //    }
        }


        //deviceMemory<dlong> o_combineFlag = platform.malloc<dlong>(CombineFlag);
        //o_q.copyFrom(Q);
        //deviceMemory<dfloat> o_qOld = platform.malloc<dfloat>(Q);
        //o_IntFlag = platform.malloc<dlong>(IntFlag);  
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 
        //combineKernel(mesh.Nelements,Ncoarse,o_qOld ,o_q, o_combineFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_RM,level);        
        
        //o_q.copyTo(Q);

        for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
        {
          EToRefLevel_new[e] = 0; 
          RefFlag_new[e] = 0;

          for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    Q_new[id] = 0.f;
                    }

          for (int n = 0; n < (level+3); ++n)
          {
          dlong id = e*(level+3)+n;

          PToC_new[id] = -1;
          IntFlag_new[id] = 0;
          }         
        }
        RefFlag= RefFlag_new;
        PToC = PToC_new;
        EToRefLevel = EToRefLevel_new;
        IntFlag = IntFlag_new;
        Q_new.copyTo(Q); 
        o_PToC = platform.malloc<dlong>(PToC);

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        //mesh.PmlSetup();
        //deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        //o_oldq = platform.malloc<dfloat>(Q);
        //o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
        o_q.copyFrom(Q);
      }         
}

void adaptivity_t::CoarseGPU(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         deviceMemory<dlong>& o_RefFlag,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(16*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(16*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(16*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(16*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(16*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(16*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(16*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(16*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  
  // Permutation map
  dlong perm[16*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[16*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(16*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < level+3; ++i)
    {
      const dlong id = e*(level+3)+i;
      PToC_new[id] = PToC[id];
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*64,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming
//    dlong ii = 0;
//  for (int e = 0; e < mesh.Nelements; ++e)
//  {
//    if (RefFlag[e]==-1)
//    {
//      Coar[ii] = e;
//      ii = ii + 1;
//    }
//  }

  const dlong hierarchyStride = level+3;
  deviceMemory<dlong> o_initialCandidate = platform.reserve<dlong>(mesh.Nelements);
  deviceMemory<dlong> o_coarseList = platform.reserve<dlong>(4*mesh.Nelements);
  deviceMemory<long long int> o_DeleteFlag = platform.reserve<long long int>(mesh.Nelements);
  deviceMemory<dlong> o_NcoarseOut = platform.reserve<dlong>(1);

  deviceMemory<dfloat> o_qold = platform.malloc<dfloat>(Q);
  deviceMemory<long long int> o_EToV_new = platform.reserve<long long int>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dlong> o_EToB_new = platform.reserve<int>(16*mesh.Nelements*mesh.Nverts);
  o_EToE = platform.malloc<long long int>(mesh.EToE);
  o_EToV_new = platform.malloc<long long int>(EToV_new);
  o_EToV = platform.malloc<long long int>(EToV_new);

  o_EToB_new = platform.malloc<int>(EToB_new);
  o_EToRefLevel = platform.malloc<dlong>(EToRefLevel);
  o_EX = platform.malloc<dfloat>(mesh.EX);
  o_EY = platform.malloc<dfloat>(mesh.EY);
  deviceMemory<dfloat> o_EX_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  deviceMemory<dfloat> o_EY_new = platform.malloc<dfloat>(16*mesh.Nelements*mesh.Nverts);
  o_EX_new = platform.malloc<dfloat>(EX_new);
  o_EY_new = platform.malloc<dfloat>(EY_new);

  memory<dlong> NcoarseOut(1, 0);
  //deviceMemory<dlong> o_NcoarseOut =platform.malloc<dlong>(NcoarseOut);
  
  memory<dlong> checkRef(mesh.Nelements);
memory<dlong> checkLevel(mesh.Nelements);
memory<dlong> checkPToC(mesh.Nelements*hierarchyStride);
//memory<dlong> checkInt(128*mesh.Nelements*(5+3));

o_RefFlag.copyTo(checkRef, mesh.Nelements);
o_EToRefLevel.copyTo(checkLevel, mesh.Nelements);
o_PToC.copyTo(checkPToC, mesh.Nelements*hierarchyStride);
o_IntFlag.copyTo(IntFlag);

dlong marked = 0;
dlong validHierarchy = 0;
dlong validRule = 0;

for(dlong e=0; e<mesh.Nelements; ++e) {
  if(checkRef[e] != -1)
    continue;

  marked++;

  const dlong L = checkLevel[e];
  const dlong parent =
      (L > 0) ? checkPToC[e*hierarchyStride] : -1;

  const dlong sibling =
      (parent >= 0 && parent < mesh.Nelements)
      ? checkPToC[parent*hierarchyStride+L] : -1;

  dlong rule = -1;
  if(parent >= 0 && parent < mesh.Nelements && L > 0)
    rule = IntFlag[parent*hierarchyStride+L-1];

  const bool hierarchyOK =
      L > 0 &&
      parent >= 0 && parent < mesh.Nelements &&
      sibling >= 0 && sibling < mesh.Nelements &&
      checkLevel[parent] == L &&
      checkLevel[sibling] == L;

  if(hierarchyOK)
    validHierarchy++;

  if(rule == 1 || rule == 3 || rule == 5)
    validRule++;

  printf("candidate e=%d L=%d parent=%d sibling=%d "
         "parentL=%d siblingL=%d rule=%d\n",
         e, L, parent, sibling,
         parent >= 0 ? checkLevel[parent] : -1,
         sibling >= 0 ? checkLevel[sibling] : -1,
         rule);
}

printf("marked=%d validHierarchy=%d validRule=%d\n",
       marked, validHierarchy, validRule);

   coarseCandidateKernel(mesh.Nelements,hierarchyStride,o_RefFlag,o_EToRefLevel,
                  o_PToC,o_IntFlag,o_EToE,o_initialCandidate,
                  o_coarseList,o_NcoarseOut);

   o_NcoarseOut.copyTo(NcoarseOut);

   const dlong Npairs = NcoarseOut[0];               
    const dlong nn = Npairs;

    printf("Elements to be coarsened %d \n", nn);
    // Determine Triangles To be Coarsened
  hlong del_vertex = 0; // Counts each new_vertex that will be deleted
  //#pragma omp parallel for
 // ii = 0;
 // for (int e = 0; e < mesh.Nelements; ++e)
 // {
 //   if (RefFlag[e]==-1)
 //   {
 //     Coar[ii] = e;
 //     ii = ii + 1;
 //   }
 // }
  //printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

  coarseKernel(Npairs,hierarchyStride,o_coarseList,o_q,o_qold,
         o_RefFlag,o_EToRefLevel,o_PToC,o_IntFlag,
         o_EX,o_EY,o_EToV,mesh.o_EToB,
         o_EX_new,o_EY_new,o_EToV_new,o_EToB_new,
         o_DeleteFlag,o_RM);

// Loop for updating mesh information if coarsening done.
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement GPU\n"); 
        
        // Map ids of elements to new ones
        
        //memory<dlong> DeleteFlag(mesh.Nelements, 0);
        //deviceMemory<dlong> o_DeleteFlag = platform.malloc<dlong>(DeleteFlag);
        deviceMemory<dlong> o_perm = platform.reserve<dlong>(mesh.Nelements);
        memory<dlong> NelementsOut(1, 0);
        deviceMemory<dlong> o_NelementsOut = platform.malloc<dlong>(NelementsOut);

        deviceMemory<dfloat> o_qCompact = platform.malloc<dfloat>(Q);
        deviceMemory<dfloat> o_EXCompact = platform.malloc<dfloat>(mesh.EX);
        deviceMemory<dfloat> o_EYCompact = platform.malloc<dfloat>(mesh.EY);
        deviceMemory<long long int> o_EToVCompact = platform.malloc<long long int>(EToV_new);
        deviceMemory<int> o_EToBCompact = platform.malloc<int>(EToB_new);
        deviceMemory<dlong> o_RefFlagCompact = platform.reserve<dlong>(mesh.Nelements);
        deviceMemory<dlong> o_EToRefLevelCompact = platform.malloc<dlong>(EToRefLevel);
        deviceMemory<dlong> o_PToCCompact = platform.malloc<dlong>(PToC);
        deviceMemory<dlong> o_IntFlagCompact = platform.malloc<dlong>(IntFlag);
        
        permKernel(mesh.Nelements,o_DeleteFlag,o_perm,
                   o_NelementsOut);  

        compactKernel(mesh.Nelements,hierarchyStride,o_perm,
                      o_q,o_qold,o_EX_new,o_EY_new,
                      o_EToV_new,o_EToB_new,o_RefFlag,
                      o_EToRefLevel,o_PToC,o_IntFlag,
                      o_qCompact,o_EXCompact,o_EYCompact,
                      o_EToVCompact,o_EToBCompact,o_RefFlagCompact,o_EToRefLevelCompact,
                      o_PToCCompact,o_IntFlagCompact);

        o_NelementsOut.copyTo(NelementsOut);

        const dlong newNelements = NelementsOut[0];
        const dlong nn = mesh.Nelements-newNelements;

        o_EToVCompact.copyTo(EToV_new);
        o_EToBCompact.copyTo(EToB_new);
        o_EXCompact.copyTo(mesh.EX);
        o_EYCompact.copyTo(mesh.EY);
        o_EToRefLevelCompact.copyTo(EToRefLevel);
        
        mesh.EToV = EToV_new;
        mesh.EToB = EToB_new;
        o_qCompact.copyTo(Q); 
       
        o_RefFlagCompact.copyTo(o_RefFlag);
        o_PToCCompact.copyTo(o_PToC);
        o_IntFlagCompact.copyTo(o_IntFlag);
        
        mesh.Nelements = newNelements;
        mesh = mesh.SetupUpdate(nn);
                            
        //mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",newNelements,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        //mesh.PmlSetup();
        o_q.copyFrom(Q);
      }         
}


void adaptivity_t::CoarseGreentoRed(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         memory<dlong>& RedFlag,
                         memory<dlong>& ConfFlag,
                         memory<dlong>& ConfGreen,
                         memory<hlong>& EToNewV,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  
  printf("Inside CoarseGreentoRed loop=%d\n",Ncoarse);
  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);

  const dlong stride = 16*level;

  memory<dlong>RefFlag_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>RedFlag_new(8*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>ConfFlag_new(8*mesh.Nelements,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>ConfGreen_new(6*mesh.Nelements*8,-1); //                       
  memory<dlong>PToC_new(8*mesh.Nelements*stride,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(8*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<hlong>EToNewV_new(8*mesh.Nelements*3,-1);
  // Permutation map
  dlong perm[8*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[8*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(8*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < stride; ++i)
    {
      const dlong id = e*stride+i;
      PToC_new[id] = PToC[id];
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*2,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming

   for (int e = 0; e < mesh.Nelements; ++e)
  {
     dlong const L = EToRefLevel[e];
     dlong id = e*mesh.Nverts; 
     dlong marker = 0;

    if (RefFlag[e]==1 && (L<level+1 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1)  )
    {
      
      printf("redflag=%d",RedFlag[e*(level+3)+(EToRefLevel[e])-1]);  
    if (EToRefLevel[e]>0)
    {
     dlong rep = PToC[e*(stride)+(L-1)*4+0];
     dlong sib = PToC[e*(stride)+(L-1)*4+1];
      if (rep == e && sib!=-1)
      {
        dlong n0 = mesh.EToE[id+0];
        dlong n1 = mesh.EToE[id+1];
        dlong n2 = mesh.EToE[id+2];

        if(n0>-1 && PToC[n0*(stride)+(L-1)*4+1]==e) {RefFlag[n0]=1; marker = 1;}
        else if(n1>-1 && PToC[n1*(stride)+(L-1)*4+1]==e) {RefFlag[n1]=1; marker = 1;}
        else if(n2>-1 && PToC[n2*(stride)+(L-1)*4+1]==e) {RefFlag[n2]=1; marker = 1;}
      }
      else if (rep==PToC[rep*(stride)+(L-1)*4+0] && sib==-1) {RefFlag[rep]=1; marker = 1;}
    }
    }

  } 

  // Checking older parents for double blue refinement
    for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==1 && (EToRefLevel[e]<level+3 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1) )
    { 
          if (EToRefLevel[e]>1 && e== PToC[e*(stride)+(EToRefLevel[e]-1)*4+0])
    { 
      dlong old_parent = PToC[e*(stride)+(EToRefLevel[e]-2)*4+0];
      dlong old_sibling = PToC[old_parent*(stride)+(EToRefLevel[e]-2)*4+1];
      if(old_sibling!=-1 && EToRefLevel[old_sibling]==EToRefLevel[e]) {RefFlag[old_sibling]=1; 
      printf("old_sibling %d of element %d is marked\n",old_sibling,e);}
    if(old_parent!=-1 && EToRefLevel[old_parent]==EToRefLevel[e]) {RefFlag[old_parent]=1;
     printf("old_parent %d of element %d is marked\n",old_parent,e);} 
  } 
  if(EToRefLevel[e]>1 && e!=PToC[e*(stride)+(EToRefLevel[e]-1)*4+0] && PToC[e*(stride)+(EToRefLevel[e]-1)*4+0]!=-1) {
    dlong parent = PToC[e*(stride)+(EToRefLevel[e]-1)*4+0];
    dlong old_parent = PToC[parent*(stride)+(EToRefLevel[parent]-2)*4+0];
    dlong old_sibling = PToC[old_parent*(stride)+(EToRefLevel[parent]-2)*4+1];
      
       if(old_sibling!=-1 && EToRefLevel[old_sibling]==EToRefLevel[parent]) {RefFlag[old_sibling]=1; 
      printf("old_sibling %d of element %d is marked\n",old_sibling,parent);}
     if(old_parent!=-1 && EToRefLevel[old_parent]==EToRefLevel[parent]) {RefFlag[old_parent]=1;
     printf("old_parent %d of element %d is marked\n",old_parent,parent);} 

  }

    }
  }

// Checking double green blue refine elements 
//NEED THIS LATER IT LL PROVIDE ONLY LAST LEVEL IS GREEN/BLUE
//for (int e = 0; e < mesh.Nelements; ++e)
//  {
//    if (EToRefLevel[e]>1 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1
//    && RedFlag[e*(level+3)+(EToRefLevel[e])-2]!=1)
//    { 
//          if (e== PToC[e*(stride)+(EToRefLevel[e]-1)*4+0])
//    { dlong sibling = PToC[e*(stride)+(EToRefLevel[e]-2)*4+1];
//      dlong old_parent = PToC[e*(stride)+(EToRefLevel[e]-2)*4+0];
//      dlong old_sibling = PToC[old_parent*(stride)+(EToRefLevel[e]-2)*4+1];
//      if(old_sibling!=-1 && EToRefLevel[old_sibling]==EToRefLevel[e] && sibling!=old_sibling) 
//      {RefFlag[old_sibling]=1;RefFlag[e]=1;}
//    if(old_parent!=-1 && EToRefLevel[old_parent]==EToRefLevel[e]) {RefFlag[old_parent]=1;
//     RefFlag[e]=1;} 
//  } 
//    }
//  }
    for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==1 && (EToRefLevel[e]==level+1 && RedFlag[e*(level+3)+(EToRefLevel[e])-2]==1) )
    {
          RefFlag[e]=0;
    }
  }
  
    dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==1 && (EToRefLevel[e]<level+3 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1) )
    {
          Coar[ii] = e;
      ii = ii + 1;
      printf("elements to be green to red e=%d \n",e);
    }
  }


    // Determine Triangles To be Coarsened
  hlong del_vertex = 0; // Counts each new_vertex that will be deleted
  //#pragma omp parallel for
 
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

        const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(6*Np*Np,0);
      o_RM.copyTo(RM);



  
  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 
    dlong const L = EToRefLevel[e];
    if(L<=0) continue;
    dlong rep = PToC[e*(stride)+(L-1)*4+0];
    dlong sib = PToC[e*(stride)+(L-1)*4+1];
    
    printf("element %d in the loop for GToR with rep=%d and sib=%d \n",e,rep,sib);
    if (rep!=-1 && sib!=-1 && L>0)
    {
      
      hlong sib_e=0 ;
      hlong sib_id=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      hlong v0_sib=0;
      hlong v1_sib=0;
      hlong v2_sib=0;
      
      if (e==rep) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = sib;
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = rep;
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
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e])
        {

          // Rule 0, IntFlag[e]==1
          if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          {    
          
          EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+3]= EToNewV[e*4+1]; EToNewV[e*4+1]=mesh.EToV[id+2]; printf("Rule 0 Double Coarsening applied to e=%d \n",e);} 
          if (EToNewV[sib_e*4+1]!=-1 ) 
          {EToNewV[e*4+2]= EToNewV[sib_e*4+1]; EToNewV[e*4+1]=mesh.EToV[id+2]; printf("Rule 0 Double Coarsening applied from sib to e=%d \n",e);}
          if (EToNewV[e*4+1]==-1 ) 
          {EToNewV[e*4+1] = mesh.EToV[id+2] ; printf("Rule 0 First Coarsening applied to e=%d \n",e);}  
          
           
         
          
          if (ConfGreen[sib_e*6+0]!=-1)  
          {ConfGreen[e*6+2]=ConfGreen[sib_e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+0]);}
          if (ConfGreen[sib_e*6+1]!=-1)  
          {ConfGreen[e*6+3]=ConfGreen[sib_e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+1]);}    

          if (ConfGreen[e*6+0]!=-1)  
          {ConfGreen[e*6+4]=ConfGreen[e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+0]);}
          if (ConfGreen[e*6+1]!=-1)  
          {ConfGreen[e*6+5]=ConfGreen[e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+1]);}    

          if (mesh.EToE[id+1]!=-1)
          { hlong n0 = mesh.EToE[id+1];
            // if the neigh is green/blue refined parent should be stored
            if (RedFlag[n0*(level+3)+EToRefLevel[n0]-1]!=1) n0 = PToC[n0*(stride)+(EToRefLevel[n0]-1)*4+0]; 
            ConfGreen[e*6+0] = n0;   
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[id+1]); 
          }
          if (mesh.EToE[sib_id+2]!=-1)
          { hlong n1 = mesh.EToE[sib_id+2];
             if (RedFlag[n1*(level+3)+EToRefLevel[n1]-1]!=1) n1 = PToC[n1*(stride)+(EToRefLevel[n1]-1)*4+0];
            ConfGreen[e*6+1] = n1;   
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[sib_id+2]); 
          }
          

          EToV_new[id+0] = v1;
          EToV_new[id+1] = v0_sib;  
          EToV_new[id+2] = v1_sib;

          EX_new[id+0] = EX_new[id+1];
          EX_new[id+1] = EX_new[sib_id+0];
          EX_new[id+2] = EX_new[sib_id+1];
          
          EY_new[id+0] = EY_new[id+1];
          EY_new[id+1] = EY_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id+2];
          EToB_new[id+1] = mesh.EToB[sib_id+0];
          EToB_new[id+2] = mesh.EToB[id+0];

          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          
       
 
          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 0;
          const dlong id_int_c = 1;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
          printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          }

          // Rule 1, IntFlag[e]==5
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==5)
          {    
          
          dlong marker = 0;  
          EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+1]= EToNewV[e*4+1]; printf("Rule 1 Double Coarsening applied to e=%d \n",e);
          marker = 1;} 
          if (EToNewV[sib_e*4+1]!=-1) 
          {EToNewV[e*4+3]= EToNewV[sib_e*4+1];  printf("Rule 1 Double Coarsening applied from sib to e=%d \n",e);
          marker = 1;}
          if (EToNewV[e*4+2]==-1)
          { 
            EToNewV[e*4+2] = mesh.EToV[id+2] ;printf("Rule 1 First Coarsening applied to e=%d \n",e);
          if((EToRefLevel[e]>1 && RedFlag[e*(level+3)+EToRefLevel[e]-2]!=1)&& marker!=1) {EToNewV[e*4+2] = -1;}
          }  
          

          if (ConfGreen[sib_e*6+0]!=-1)  
          {ConfGreen[e*6+2]=ConfGreen[sib_e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+0]);}
          if (ConfGreen[sib_e*6+1]!=-1)  
          {ConfGreen[e*6+3]=ConfGreen[sib_e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+1]);}    

          if (ConfGreen[e*6+0]!=-1)  
          {ConfGreen[e*6+4]=ConfGreen[e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+0]);}
          if (ConfGreen[e*6+1]!=-1)  
          {ConfGreen[e*6+5]=ConfGreen[e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+1]);}    

          if (mesh.EToE[id+1]!=-1)
          { hlong n0 = mesh.EToE[id+1];
            // if the neigh is green/blue refined parent should be stored
            if (RedFlag[n0*(level+3)+EToRefLevel[n0]-1]!=1) n0 = PToC[n0*(stride)+(EToRefLevel[n0]-1)*4+0]; 
          ConfGreen[e*6+0] = n0;   
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[id+1]); 
          }
          if (mesh.EToE[sib_id+2]!=-1)
          {hlong n1 = mesh.EToE[sib_id+2];
             if (RedFlag[n1*(level+3)+EToRefLevel[n1]-1]!=1) n1 = PToC[n1*(stride)+(EToRefLevel[n1]-1)*4+0];
            ConfGreen[e*6+1] = n1;     
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[sib_id+2]); 
          }

          EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
                    EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          
        
   
          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                          for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 4;
          const dlong id_int_c = 5;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          }

          // Rule 2, IntFlag[e]==3
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==3)
          {    
          
          dlong marker = 0;   
          EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+1]= EToNewV[e*4+1]; printf("Rule 2 Double Coarsening applied to e=%d \n",e);
          marker=1;} 
          if (EToNewV[sib_e*4+1]!=-1) 
          {EToNewV[e*4+2]= EToNewV[sib_e*4+1];  printf("Rule 2 Double Coarsening applied from sib to e=%d \n",e);
          marker=1;}
          if (EToNewV[e*4+3]==-1)
          {EToNewV[e*4+3] = mesh.EToV[id+2] ;printf("Rule 2 First Coarsening applied to e=%d \n",e);
            if((EToRefLevel[e]>1 && RedFlag[e*(level+3)+EToRefLevel[e]-2]!=1) && marker!=1) {EToNewV[e*4+3] = -1;}}

          if (ConfGreen[sib_e*6+0]!=-1)  
          {ConfGreen[e*6+2]=ConfGreen[sib_e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+0]);}
          if (ConfGreen[sib_e*6+1]!=-1)  
          {ConfGreen[e*6+3]=ConfGreen[sib_e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[sib_e*6+1]);}    

          if (ConfGreen[e*6+0]!=-1)  
          {ConfGreen[e*6+4]=ConfGreen[e*6+0];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+0]);}
          if (ConfGreen[e*6+1]!=-1)  
          {ConfGreen[e*6+5]=ConfGreen[e*6+1];
          printf("Rule 0 second conform call applied to e=%d by  %d\n",e,ConfGreen[e*6+1]);}    

          if (mesh.EToE[id+2]!=-1)
          { hlong n0 = mesh.EToE[id+2];
            // if the neigh is green/blue refined parent should be stored
            if (RedFlag[n0*(level+3)+EToRefLevel[n0]-1]!=1) n0 = PToC[n0*(stride)+(EToRefLevel[n0]-1)*4+0]; 
            ConfGreen[e*6+0] = n0;   
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[id+1]); 
          }
          if (mesh.EToE[sib_id+1]!=-1)
          {hlong n1 = mesh.EToE[sib_id+1];
            // if the neigh is green/blue refined parent should be stored
            if (RedFlag[n1*(level+3)+EToRefLevel[n1]-1]!=1) n1 = PToC[n1*(stride)+(EToRefLevel[n1]-1)*4+0]; 
            ConfGreen[e*6+1] = n1;   
          printf("Rule 0 first conform call applied to e=%d by  %d\n",e,mesh.EToE[sib_id+2]); 
          }
          
          EToV_new[id+2] = v1_sib;

          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          
       
         
          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

           for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 2;
          const dlong id_int_c = 3;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);
          } 
        }
    }             
  }   

 // Loop for updating mesh information if coarsening done.

if (Ncoarse!=0 && nn!=0)
{ 
  printf("inside coarsement\n"); 
 
  // Map ids of elements to new ones
  int e_new = 0; 
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (Delete_Flag[e]==1)
    {
      perm[e] = -1;
    }else{
      perm[e] = e_new;
      e_new++;
    }
  }
       
  e_new = 0; 
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (Delete_Flag[e]==1)
    {
      // deleted
    }
    else
    {
      if (perm[e] >= 0)
      {
        for (int n = 0; n < 3; ++n)
        {
          dlong id = e*mesh.Nverts+n;
          dlong id_new = perm[e]*mesh.Nverts+n; 
          
          mesh.EToB[id_new] = EToB_new[id];
          mesh.EToV[id_new] = EToV_new[id];
          mesh.EX[id_new]   = EX_new[id]; 
          mesh.EY[id_new]   = EY_new[id];
        }

        for (int n = 0; n < mesh.Np; ++n)
        {
          dlong id = e*mesh.Np+n;
          dlong id_new = perm[e]*mesh.Np+n; 
          Q_new[id_new] = Q[id];
        }

        for (int i = 0; i < stride; ++i)
        {
          dlong oldVal = PToC[e*stride+i];
          dlong newVal = -1;
        
          if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
            newVal = perm[oldVal];
        
          PToC_new[perm[e]*stride+i] = newVal;
        }

        for (int i = 0; i < (level+3); ++i)
        {
          IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
          RedFlag_new[perm[e]*(level+3)+i] = RedFlag[e*(level+3)+i];
        }

        EToRefLevel_new[perm[e]] = EToRefLevel[e];
        RefFlag_new[perm[e]]     = RefFlag[e];
        
        dlong oldConf = ConfFlag[e];
        dlong newConf = -1;

        if (oldConf >= 0 && oldConf < mesh.Nelements && perm[oldConf] >= 0)
          newConf = perm[oldConf];

        ConfFlag_new[perm[e]] = newConf;       
        

        // ------------------------------------------
        // Added: carry EToNewV through permutation
        // Layout: [elem, mid01, mid12, mid20]
        // ------------------------------------------
        EToNewV_new[perm[e]*4+0] = perm[e];
        EToNewV_new[perm[e]*4+1] = EToNewV[e*4+1];
        EToNewV_new[perm[e]*4+2] = EToNewV[e*4+2];
        EToNewV_new[perm[e]*4+3] = EToNewV[e*4+3];

        // ------------------------------------------
        // Added: carry ConfGreen through permutation
        // Layout: [n0, n1, n2, n3, n4, n5]
        // ------------------------------------------
for (int k=0; k<6; ++k) {
  dlong oldVal = ConfGreen[e*6+k];
  dlong newVal = -1;
  if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
    newVal = perm[oldVal];
  ConfGreen_new[perm[e]*6+k] = newVal;
}

        e_new++;  
      }
    }
  }

  for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
  {
    EToRefLevel_new[e] = 0; 
    RefFlag_new[e] = 0;
    ConfFlag_new[e] = -1;
    

    // ------------------------------------------
    // Added: clear compacted tail of EToNewV_new
    // ------------------------------------------
    EToNewV_new[e*4+0] = -1;
    EToNewV_new[e*4+1] = -1;
    EToNewV_new[e*4+2] = -1;
    EToNewV_new[e*4+3] = -1;

    ConfGreen_new[e*6+0] = -1;
    ConfGreen_new[e*6+1] = -1;
    ConfGreen_new[e*6+2] = -1;
    ConfGreen_new[e*6+3] = -1;
    ConfGreen_new[e*6+4] = -1;
    ConfGreen_new[e*6+5] = -1;

    for (int n = 0; n < mesh.Np; ++n)
    {
      dlong id = e*mesh.Np+n;
      Q_new[id] = 0.f;
    }

    for (int n = 0; n < (level+3); ++n)
    {
      dlong id = e*(level+3)+n;
      IntFlag_new[id] = 0;
      RedFlag_new[id] = -1;
    }   
    
        for (int n = 0; n < stride; ++n)
    {
      dlong id = e*stride+n;
      PToC_new[id] = -1; 
    }   
  }

  RefFlag    = RefFlag_new;
  ConfFlag   = ConfFlag_new;
  RedFlag    = RedFlag_new;
  PToC       = PToC_new;
  EToRefLevel= EToRefLevel_new;
  IntFlag    = IntFlag_new;
  EToNewV    = EToNewV_new;   // Added
  ConfGreen  = ConfGreen_new;

  Q_new.copyTo(Q); 
  o_PToC = platform.malloc<dlong>(PToC);

  mesh.Nelements = mesh.Nelements-nn;
  //mesh.Nnodes = mesh.Nnodes-del_vertex;
  mesh.o_EToB = platform.malloc<int>(mesh.EToB);

  printf("first coarsement done!!\n");
  printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
  printf("Ncoarse inside coarsening loop=%lld\n",nn);
  printf("del_vertex_count=%lld\n",del_vertex);
        
  mesh = mesh.SetupUpdate(Ncoarse);
  o_q.copyFrom(Q);
  
}
      
}


void adaptivity_t::CoarseBluetoRed(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         memory<dlong>& RedFlag,
                         memory<dlong>& ConfFlag,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(2*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(2*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(2*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(2*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>RedFlag_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>ConfFlag_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(2*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(2*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(2*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  
  // Permutation map
  dlong perm[2*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[2*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(2*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < level+3; ++i)
    {
      const dlong id = e*(level+3)+i;
      PToC_new[id] = PToC[id];
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*2,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming

   for (int e = 0; e < mesh.Nelements; ++e)
  {
     dlong id = e*mesh.Nverts; 
    if (RefFlag[e]==1 && (EToRefLevel[e]<level+3 && RedFlag[e]!=1) )
    printf("redflag=%d",RedFlag[e]);     
    {
    if (EToRefLevel[e]>0)
    {
      if (PToC[e*(level+3)+0]==e && PToC[e*(level+3)+EToRefLevel[e]]==-1)
      {
        printf("INSIDE NEW LOOP e=%d\n",PToC[2*e+0]*mesh.Np);
        printf("BEFORE PToC[e*(level+3)+0]=%d,PToC[e*(level+3)+EToRefLevel[e]]=%d, EToRefLevel[e]=%d\n", PToC[e*(level+3)+0],PToC[e*(level+3)+EToRefLevel[e]], EToRefLevel[e]);
        if(PToC[mesh.EToE[id+0]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+0]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+0]]=1;}
        else if(PToC[mesh.EToE[id+1]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+1]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+1]]=1;}
        else if(PToC[mesh.EToE[id+2]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+2]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+2]]=1;}
        printf("AFTER PToC[e*(level+3)+0]=%d,PToC[e*(level+3)+EToRefLevel[e]]=%d, EToRefLevel[e]=%d\n", PToC[e*(level+3)+0],PToC[e*(level+3)+EToRefLevel[e]], EToRefLevel[e]);
      }
      else if (PToC[e*(level+3)+0]==PToC[PToC[e*(level+3)+0]*(level+3)+0] && PToC[e*(level+3)+EToRefLevel[e]]==-1) {RefFlag[PToC[e*(level+3)+0]]=1;}
    }
    }
  } 

    dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==1 && (EToRefLevel[e]<level+3 && RedFlag[e]==3) )
    {
      Coar[ii] = e;
      ii = ii + 1;
      printf("PToC[e*(level+3)+0]=%d,PToC[e*(level+3)+EToRefLevel[e]]=%d, EToRefLevel[e]=%d\n", PToC[e*(level+3)+0],PToC[e*(level+3)+EToRefLevel[e]], EToRefLevel[e]);
    }
  }


    // Determine Triangles To be Coarsened
  hlong del_vertex = 0; // Counts each new_vertex that will be deleted
  //#pragma omp parallel for
 
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

        const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(6*Np*Np,0);
      o_RM.copyTo(RM);



  
  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 


    if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && EToRefLevel[e]>0)
    {
      
      hlong sib_e=0 ;
      hlong sib_id=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      hlong v0_sib=0;
      hlong v1_sib=0;
      hlong v2_sib=0;
      printf("e=%d\n",PToC[2*e+0]*mesh.Np);
      
      if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+EToRefLevel[e]];
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = PToC[e*(level+3)+0];
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
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e])
        {

          // Rule 0, IntFlag[e]==1
          if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          {    
          EToV_new[id+0] = v1;
          EToV_new[id+1] = v0_sib;  
          EToV_new[id+2] = v1_sib;

          EX_new[id+0] = EX_new[id+1];
          EX_new[id+1] = EX_new[sib_id+0];
          EX_new[id+2] = EX_new[sib_id+1];
          
          EY_new[id+0] = EY_new[id+1];
          EY_new[id+1] = EY_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id+2];
          EToB_new[id+1] = mesh.EToB[sib_id+0];
          EToB_new[id+2] = mesh.EToB[id+0];

          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 0;
          const dlong id_int_c = 1;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
          printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          }

          // Rule 1, IntFlag[e]==5
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==5)
          {    
          EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
                    EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                          for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 4;
          const dlong id_int_c = 5;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          }

          // Rule 2, IntFlag[e]==3
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==3)
          {    
          EToV_new[id+2] = v1_sib;

          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=1;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

           for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 2;
          const dlong id_int_c = 3;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);
          } 
        }
    }             
  }   

  // Loop for updating mesh information if coarsening done.
                printf("PToC[52]_Before=%d,%d,%d,%d,%d\n",PToC[52*(level+3)+0],PToC[52*(level+3)+1],PToC[52*(level+3)+2],PToC[52*(level+3)+3],PToC[52*(level+3)+4] );
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        // Map ids of elements to new ones
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
          if (Delete_Flag[e]==1)
          {
            perm[e] = -1;
          }else{
            perm[e] = e_new;
            e_new++;
          }
        }

        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              //printf("id=%d\n",e);
            } else{
                   // EToRefLevel[e_new] = EToRefLevel[e]; 
                    //IntFlag[e_new] = IntFlag[e];
                    if (perm[e]>= 0)
                   {
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = perm[e]*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id];
                    }

                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = perm[e]*mesh.Np+n; 
                    
                    Q_new[id_new] = Q[id];
                    
                    //printf("EToV=%lld\n",mesh.EToV[id]);
                    //printf("EX=%g\n",mesh.EX[id]);
                    //printf("EY=%f\n",mesh.EY[id]);
                    
                    }
                    
                    for (int i = 0; i < level+3; ++i) {
                    dlong oldVal = PToC[e*(level+3)+i];
                    dlong newVal = 0;
                  
                    if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
                      newVal = perm[oldVal];     // valid surviving element
                    else
                      newVal = -1;                // no link / deleted
                  
                    PToC_new[perm[e]*(level+3)+i] = newVal;
                    IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
                  }

                     EToRefLevel_new[perm[e]] = EToRefLevel[e];
                     RefFlag_new[perm[e]] = RefFlag[e];
                     ConfFlag_new[perm[e]] = ConfFlag[e];
                     RedFlag_new[perm[e]]  = RedFlag[e];

                    e_new++;  

                    }
                    }
 
            //        e_new++;
              //    }
        }


        //deviceMemory<dlong> o_combineFlag = platform.malloc<dlong>(CombineFlag);
        //o_q.copyFrom(Q);
        //deviceMemory<dfloat> o_qOld = platform.malloc<dfloat>(Q);
        //o_IntFlag = platform.malloc<dlong>(IntFlag);  
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 
        //combineKernel(mesh.Nelements,Ncoarse,o_qOld ,o_q, o_combineFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_RM,level);        
        
        //o_q.copyTo(Q);

        for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
        {
          EToRefLevel_new[e] = 0; 
          RefFlag_new[e] = 0;
          ConfFlag_new[e] = -1;
          RedFlag_new[e] = 0;

          for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    Q_new[id] = 0.f;
                    }

          for (int n = 0; n < (level+3); ++n)
          {
          dlong id = e*(level+3)+n;

          PToC_new[id] = -1;
          IntFlag_new[id] = 0;
          }         
        }
        RefFlag= RefFlag_new;
        ConfFlag = ConfFlag_new;
        RedFlag = RedFlag_new;
        PToC = PToC_new;
        EToRefLevel = EToRefLevel_new;
        IntFlag = IntFlag_new;
        Q_new.copyTo(Q); 
        o_PToC = platform.malloc<dlong>(PToC);

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        //mesh.PmlSetup();
        //deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        //o_oldq = platform.malloc<dfloat>(Q);
        //o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
        o_q.copyFrom(Q);
      }
      
}

void adaptivity_t::CoarseGreen(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         memory<dlong>& RedFlag,
                         memory<dlong>& ConfFlag,
                         memory<hlong>& EToNewV,
                         dlong Ncoarse,dlong level){

  // Store old info & Allocate new arrays
  

  printf("CoarseGreen inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 
  const dlong stride =level*16;
  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>GreenFlag(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>RedFlag_new(8*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>ConfFlag_new(8*mesh.Nelements,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(8*mesh.Nelements*stride,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(8*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<hlong>EToNewV_new(8*mesh.Nelements*3,-1);
  
  // Permutation map
  dlong perm[8*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[8*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(8*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < stride; ++i)
    {
      const dlong id = e*stride+i;
      PToC_new[id] = PToC[id];
    }
  }

  // Determine elements to be coarsened by using Coarse Flag
  
  memory<dlong> Coar(Ncoarse*2,0); // Array holds element ids for refining. Holds some extra mem. for 
                                  // conforming

   for (int e = 0; e < mesh.Nelements; ++e)
  {
     dlong id = e*mesh.Nverts; 
     dlong const L = EToRefLevel[e]; if (L<=0) continue;
     dlong rep = PToC[e*(stride)+(L-1)*4+0];
     dlong sib = PToC[e*(stride)+(L-1)*4+1];
      dlong marker = 0;
      dlong n0 = mesh.EToE[id+0];
      dlong n1 = mesh.EToE[id+1];
      dlong n2 = mesh.EToE[id+2];
            if (e==111)
      {
        printf("111 in the loop with n0=%d, n1=%d, n2=%d \n",n0,n1,n2);
        printf("111 in the loop with rep=%d, sib=%d \n",rep,sib);
      }

                  if (e==475)
      {
        printf("475 in the loop with n0=%d, n1=%d, n2=%d \n",n0,n1,n2);
        printf("475 in the loop with rep=%d, sib=%d \n",rep,sib);
      }
    if (RefFlag[e]==-1 && (L<level+3 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1) )
    {

      
    if (L>0)
    {
 
        if (n0==-1 || n1==-1 || n2==-1) {
      if (rep==e && sib==-1)
      {

        if(n0>-1 && PToC[n0*(stride)+(L-1)*4+1]==e) {GreenFlag[n0]=-1; marker = 1;
        printf("1.element %d marked for green coarsening with n0=%d, n1=%d, n2=%d \n",e,n0,n1,n2);}
        else if(n1>-1 && PToC[n1*(stride)+(L-1)*4+1]==e) {GreenFlag[n1]=-1; marker = 1;
        printf("2.element %d marked for green coarsening with n0=%d, n1=%d, n2=%d \n",e,n0,n1,n2);}
        else if(n2>-1 && PToC[n2*(stride)+(L-1)*4+1]==e) {GreenFlag[n2]=-1; marker = 1;
        printf("3.element %d marked for green coarsening with n0=%d, n1=%d, n2=%d \n",e,n0,n1,n2);} 
      }
      else if (EToRefLevel[rep] == L && rep==PToC[rep*(stride)+(L-1)*4+0]  ) {GreenFlag[rep]=-1; marker = 1;
      printf("4.element %d marked for green coarsening with rep=%d,n0=%d, n1=%d, n2=%d \n",e,rep ,n0,n1,n2);} 
    }
       
    }
    if (GreenFlag[e]==-1) printf("5.element %d marked for green coarsening with n0=%d, n1=%d, n2=%d \n",e,n0,n1,n2);
   
    }
    if (marker ==0) GreenFlag[e] = 0; 
  }

    dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (GreenFlag[e]==-1 && (EToRefLevel[e]<level+3 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]!=1 && EToRefLevel[e]>0) )
    {
      Coar[ii] = e;
      ii = ii + 1;
      printf("element %d saved for green coarsening \n",e);
    }
  }


    // Determine Triangles To be Coarsened
  hlong del_vertex = 0; // Counts each new_vertex that will be deleted
  //#pragma omp parallel for
 
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

        const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(6*Np*Np,0);
      o_RM.copyTo(RM);

  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 
    dlong const L = EToRefLevel[e]; if (L<=0) continue;
    dlong rep = PToC[e*(stride)+(L-1)*4+0];
    dlong sib = PToC[e*(stride)+(L-1)*4+1];

    if (rep!=-1 && sib!=-1 && L>0)
    {
      
      hlong sib_e=0 ;
      hlong sib_id=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      hlong v0_sib=0;
      hlong v1_sib=0;
      hlong v2_sib=0;
      printf("e=%d\n",rep);
      
      if (e==rep) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e = sib;
         sib_id = sib_e*mesh.Nverts;
         v0 = mesh.EToV[id+0]; 
         v1 = mesh.EToV[id+1]; 
         v2 = mesh.EToV[id+2];
         v0_sib = mesh.EToV[sib_id+0]; 
         v1_sib = mesh.EToV[sib_id+1]; 
         v2_sib = mesh.EToV[sib_id+2];
         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
         // Extract Vertex Number of Elements to be Combined
         sib_e = rep;
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
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e])
        {

          // Rule 0, IntFlag[e]==1
          if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          {    
          EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+3]= EToNewV[e*4+1]; EToNewV[e*4+1]=mesh.EToV[id+2]; printf("Rule 0 Double Coarsening applied to e=%d \n",e);} 
          else if (EToNewV[sib_e*4+1]!=-1 ) 
          {EToNewV[e*4+2]= EToNewV[sib_e*4+1]; EToNewV[e*4+1]=mesh.EToV[id+2]; printf("Rule 0 Double Coarsening applied from sib to e=%d \n",e);}
          else 
          {EToNewV[e*4+1] = mesh.EToV[id+2] ; printf("Rule 0 First Coarsening applied to e=%d \n",e);}

          EToV_new[id+0] = v1;
          EToV_new[id+1] = v0_sib;  
          EToV_new[id+2] = v1_sib;

          EX_new[id+0] = EX_new[id+1];
          EX_new[id+1] = EX_new[sib_id+0];
          EX_new[id+2] = EX_new[sib_id+1];
          
          EY_new[id+0] = EY_new[id+1];
          EY_new[id+1] = EY_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id+2];
          EToB_new[id+1] = mesh.EToB[sib_id+0];
          EToB_new[id+2] = mesh.EToB[id+0];

          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 0;
          const dlong id_int_c = 1;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
          printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          }

          // Rule 1, IntFlag[e]==5
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==5)
          {    

           dlong marker=0; 
           EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+1]= EToNewV[e*4+1]; printf("Rule 1 Double Coarsening applied to e=%d \n",e);
          marker = 1;} 
          if (EToNewV[sib_e*4+1]!=-1) 
          {EToNewV[e*4+3]= EToNewV[sib_e*4+1];  printf("Rule 1 Double Coarsening applied from sib to e=%d \n",e);
          marker = 1;}
          if (EToNewV[e*4+2]==-1)
          { 
            EToNewV[e*4+2] = mesh.EToV[id+2] ;printf("Rule 1 First Coarsening applied to e=%d \n",e);
          if((EToRefLevel[e]>1 && RedFlag[e*(level+3)+EToRefLevel[e]-2]!=1)&& marker!=1) {EToNewV[e*4+2] = -1;}
          }  

          EToV_new[id+2] = v0_sib;

          EX_new[id+2] = EX_new[sib_id+0];
          EY_new[id+2] = EY_new[sib_id+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+2] = mesh.EToB[id+0];
                    EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

                          for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 4;
          const dlong id_int_c = 5;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse1_=%d,sibe_coarse1_=%d\n",e,sib_e);
          }

          // Rule 2, IntFlag[e]==3
          else if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==3)
          {    
            dlong marker=0; 
  EToNewV[e*4+0] = e ; 
          if (EToNewV[e*4+1]!=-1) 
          {EToNewV[e*4+1]= EToNewV[e*4+1]; printf("Rule 2 Double Coarsening applied to e=%d \n",e);
          marker=1;} 
          if (EToNewV[sib_e*4+1]!=-1) 
          {EToNewV[e*4+2]= EToNewV[sib_e*4+1];  printf("Rule 2 Double Coarsening applied from sib to e=%d \n",e);
          marker=1;}
          if (EToNewV[e*4+3]==-1)
          {EToNewV[e*4+3] = mesh.EToV[id+2] ;printf("Rule 2 First Coarsening applied to e=%d \n",e);
            if((EToRefLevel[e]>1 && RedFlag[e*(level+3)+EToRefLevel[e]-2]!=1) && marker!=1) {EToNewV[e*4+3] = -1;}}


          EToV_new[id+2] = v1_sib;

          EX_new[id+2] = EX_new[sib_id+1];
          EY_new[id+2] = EY_new[sib_id+1];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e] = EToRefLevel[sib_e]-1;
          Delete_Flag[sib_e] = 1;

          EToB_new[id+1] = mesh.EToB[id+0];
          EToB_new[sib_id+0] = 0;
          EToB_new[sib_id+1] = 0;
          EToB_new[sib_id+2] = 0;
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

          RefFlag[e]=0;
          //RefFlag[sib_e]=1;
          CombineFlag[e]=-1;
          CombineFlag[sib_e]=-1;
          del_vertex++;
          nn++;

           for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          const dlong id_p = e*mesh.Np;
          const dlong id_c = sib_e*mesh.Np;
          const dlong id_int_p = 2;
          const dlong id_int_c = 3;
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[2*n*mesh.Np+ip+id_int_p*mesh.Np*mesh.Np+mesh.Np];
                 qn2 += Ri*Qold[id_c+ip];  
            }        
          Q[id_p+n] = qn1+qn2;
          }                  
              CombineFlag[e]=0;
              CombineFlag[sib_e]=0;
                                    for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             // printf("Q[%d]=%f\n",id1+n,Q[id1+n]);
            }
          printf("e_coarse2_=%d,sibe_coarse2_=%d\n",e,sib_e);
          } 
        }
    }             
  }   

  // Loop for updating mesh information if coarsening done.
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        // Map ids of elements to new ones
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
          if (Delete_Flag[e]==1)
          {
            perm[e] = -1;
          }else{
            perm[e] = e_new;
            e_new++;
          }
        }

        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              
              printf("id=%d\n",e);
            } else{
                   // EToRefLevel[e_new] = EToRefLevel[e]; 
                    //IntFlag[e_new] = IntFlag[e];
                    if (perm[e]>= 0)
                   {
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = perm[e]*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id];
                    }

                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = perm[e]*mesh.Np+n; 
                    
                    Q_new[id_new] = Q[id];
                                    
                    }
                    
                      for (int i = 0; i < stride; ++i)
                    {
                      dlong oldVal = PToC[e*stride+i];
                      dlong newVal = -1;
                    
                      if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
                        newVal = perm[oldVal];
                    
                      PToC_new[perm[e]*stride+i] = newVal;
                    }
                  
                    for (int i = 0; i < (level+3); ++i)
                    {
                      IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
                      RedFlag_new[perm[e]*(level+3)+i]  = RedFlag[e*(level+3)+i];
                    }

                     EToRefLevel_new[perm[e]] = EToRefLevel[e];
                     RefFlag_new[perm[e]] = RefFlag[e];
                    dlong oldConf = ConfFlag[e];
                    dlong newConf = -1;
                              
                    if (oldConf >= 0 && oldConf < mesh.Nelements && perm[oldConf] >= 0)
                      newConf = perm[oldConf];
                              
                    ConfFlag_new[perm[e]] = newConf;   
                     
                     EToNewV_new[perm[e]*4+0] = perm[e];
                     EToNewV_new[perm[e]*4+1] = EToNewV[e*4+1];
                     EToNewV_new[perm[e]*4+2] = EToNewV[e*4+2];
                     EToNewV_new[perm[e]*4+3] = EToNewV[e*4+3];

                    e_new++;  

                    }
                    }

        }


        //deviceMemory<dlong> o_combineFlag = platform.malloc<dlong>(CombineFlag);
        //o_q.copyFrom(Q);
        //deviceMemory<dfloat> o_qOld = platform.malloc<dfloat>(Q);
        //o_IntFlag = platform.malloc<dlong>(IntFlag);  
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 
        //combineKernel(mesh.Nelements,Ncoarse,o_qOld ,o_q, o_combineFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_RM,level);        
        
        //o_q.copyTo(Q);

        for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
        {
          EToRefLevel_new[e] = 0; 
          RefFlag_new[e] = 0;
          ConfFlag_new[e] = -1;
          
          EToNewV_new[e*4+0] = -1;
          EToNewV_new[e*4+1] = -1;
          EToNewV_new[e*4+2] = -1;
          EToNewV_new[e*4+3] = -1;

          for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    Q_new[id] = 0.f;
                    }

          for (int n = 0; n < (level+3); ++n)
          {
          dlong id = e*(level+3)+n;
          IntFlag_new[id] = 0;
          RedFlag_new[id] = 0;
          }  
          
          for (int n = 0; n < stride; ++n)
          {
          dlong id = e*stride+n;
          PToC_new[id] = -1;
    
          }  
        }
        RefFlag= RefFlag_new;
        ConfFlag= ConfFlag_new;
        RedFlag= RedFlag_new;
        EToNewV = EToNewV_new;
        PToC = PToC_new;
        EToRefLevel = EToRefLevel_new;
        IntFlag = IntFlag_new;
        Q_new.copyTo(Q); 
        o_PToC = platform.malloc<dlong>(PToC);

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        //mesh.PmlSetup();
        //deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        //o_oldq = platform.malloc<dfloat>(Q);
        //o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
        o_q.copyFrom(Q);
      }         
}

void adaptivity_t::CoarseRed(deviceMemory<dfloat>& o_q,
                         memory<dfloat>& Q,
                         memory<dfloat>& Qold,
                         memory<dlong>& RefFlag,
                         memory<dlong>& RedFlag,
                         memory<dlong>& ConfFlag,
                         dlong Ncoarse,dlong level){

  //
  printf("Inside Red Coarsening\n");                          
  // Store old info & Allocate new arrays
  dlong const stride = 16*level;

  printf("Ncoarse inside coarsening loop=%d\n",Ncoarse);
  // Element to vertex & Element to boundary connectivity 

  memory<hlong>EToV_new(8*mesh.Nelements*mesh.Nverts);
  memory<int>EToB_new(8*mesh.Nelements*mesh.Nverts);

  // Vertex physical coordinates
  memory<dfloat>EX_new(8*mesh.Nelements*mesh.Nverts);
  memory<dfloat>EY_new(8*mesh.Nelements*mesh.Nverts);

  memory<dlong>RefFlag_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>RedFlag_new(8*mesh.Nelements*(level+3),-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>ConfFlag_new(8*mesh.Nelements,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>PToC_new(8*mesh.Nelements*stride,-1); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>EToRefLevel_new(8*mesh.Nelements,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dlong>IntFlag_new(8*mesh.Nelements*(level+3),0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  memory<dfloat>Q_new(mesh.Nelements*mesh.Np+mesh.totalHaloPairs*mesh.Np,0); // For bisection only!(2 children from 1 parent) for 4 levels of refinement max.
  
  // Permutation map
  dlong perm[8*mesh.Nelements]={};

  // A flag to ommit informations of deleted element
  hlong Delete_Flag[8*mesh.Nelements]={};

  // A flag to be used in combine kernel, initialized with zero values. 
  memory<dlong> CombineFlag(8*mesh.Nelements,0);

  // Copy old Element to Vertex Connectivity to the New One
  //#pragma omp parallel for
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
    EToRefLevel_new[e] = EToRefLevel[e];
    for (int i = 0; i < stride; ++i)
    {
      const dlong id = e*stride+i;
      PToC_new[id] = PToC[id];
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
  //#pragma omp parallel for
  for (int i = 0; i < ii; ++i)
  {

      int e = Coar[i];
      const hlong id = e*mesh.Nverts; 
       //printf("I am Here!!!!\n");  
      //printf("refflag=%d\n",RefFlag[e]); 
      if (EToRefLevel[e]<=0)
      {
      continue;  
      }
      dlong const L = EToRefLevel[e];
      if (L<=0) continue;
      dlong rep = PToC[e*(stride)+(L-1)*4+0];
      dlong sib = PToC[e*(stride)+(L-1)*4+1];
      // For Red ONLY!!!
      if (rep!=-1 && sib!=-1 && RefFlag[e]==-1 && RedFlag[e*(level+3)+(EToRefLevel[e])-1]== 1)
      {

        if (e==rep) // Means that child stored at parents location
      {  
        // Keep the Flag and check eligibility again

        // Sibling IDS
        hlong sib_e2 = PToC[e*stride+(EToRefLevel[e]-1)*4+1];
        hlong sib_e3 = PToC[e*stride+(EToRefLevel[e]-1)*4+2];
        hlong sib_e4 = PToC[e*stride+(EToRefLevel[e]-1)*4+3];

        hlong sib_id2 = sib_e2*mesh.Nfaces;
        hlong sib_id3 = sib_e3*mesh.Nfaces;
        hlong sib_id4 = sib_e4*mesh.Nfaces;

        
        // Neighbor IDs
hlong Neigh00 = mesh.EToE[sib_id2+0];
hlong Neigh01 = mesh.EToE[sib_id3+2];

hlong Neigh10 = mesh.EToE[sib_id3+0];
hlong Neigh11 = mesh.EToE[sib_id4+2];

hlong Neigh20 = mesh.EToE[sib_id4+0];
hlong Neigh21 = mesh.EToE[sib_id2+2];

// Guarded sibling reads
dlong Neigh00_sib = -1;
dlong Neigh01_sib = -1;
dlong Neigh10_sib = -1;
dlong Neigh11_sib = -1;
dlong Neigh20_sib = -1;
dlong Neigh21_sib = -1;

if (Neigh00 > -1 && EToRefLevel[Neigh00] > 0)
  Neigh00_sib = PToC[Neigh00*stride + (EToRefLevel[Neigh00]-1)*4 + 1];

if (Neigh01 > -1 && EToRefLevel[Neigh01] > 0)
  Neigh01_sib = PToC[Neigh01*stride + (EToRefLevel[Neigh01]-1)*4 + 1];

if (Neigh10 > -1 && EToRefLevel[Neigh10] > 0)
  Neigh10_sib = PToC[Neigh10*stride + (EToRefLevel[Neigh10]-1)*4 + 1];

if (Neigh11 > -1 && EToRefLevel[Neigh11] > 0)
  Neigh11_sib = PToC[Neigh11*stride + (EToRefLevel[Neigh11]-1)*4 + 1];

if (Neigh20 > -1 && EToRefLevel[Neigh20] > 0)
  Neigh20_sib = PToC[Neigh20*stride + (EToRefLevel[Neigh20]-1)*4 + 1];

if (Neigh21 > -1 && EToRefLevel[Neigh21] > 0)
  Neigh21_sib = PToC[Neigh21*stride + (EToRefLevel[Neigh21]-1)*4 + 1];

printf("In red coarsening check: e=%d,sib_e2=%d,sib_e3=%d,sib_e4=%d\n",e,sib_e2,sib_e3,sib_e4);
printf("In red coarsening check: Neigh00=%d,Neigh10=%d,Neigh20=%d\n",Neigh00,Neigh10,Neigh20);

// Pair validity checks
bool pair00 =
  (Neigh00 > -1 && Neigh01 > -1 &&
   EToRefLevel[Neigh00] > 0 && EToRefLevel[Neigh01] > 0 &&
   EToRefLevel[Neigh00] == EToRefLevel[Neigh01] &&
   (Neigh00_sib == Neigh01 || Neigh01_sib == Neigh00));

bool pair10 =
  (Neigh10 > -1 && Neigh11 > -1 &&
   EToRefLevel[Neigh10] > 0 && EToRefLevel[Neigh11] > 0 &&
   EToRefLevel[Neigh10] == EToRefLevel[Neigh11] &&
   (Neigh10_sib == Neigh11 || Neigh11_sib == Neigh10));

bool pair20 =
  (Neigh20 > -1 && Neigh21 > -1 &&
   EToRefLevel[Neigh20] > 0 && EToRefLevel[Neigh21] > 0 &&
   EToRefLevel[Neigh20] == EToRefLevel[Neigh21] &&
   (Neigh20_sib == Neigh21 || Neigh21_sib == Neigh20));
   bool side00_ok = false;
bool side10_ok = false;
bool side20_ok = false;

// side 00
if (Neigh00 == -1) {
  side00_ok = true;   // boundary side
}
else if (RedFlag[Neigh00*(level+3)+(EToRefLevel[Neigh00])-1] == 1) {
  side00_ok = true;   // keep your current red-neighbor rule
}
else if (pair00) {
  side00_ok = true;   // blue/green side is conforming by current sibling pair
}

// side 10
if (Neigh10 == -1) {
  side10_ok = true;
}
else if (RedFlag[Neigh10*(level+3)+(EToRefLevel[Neigh10])-1] == 1) {
  side10_ok = true;
}
else if (pair10) {
  side10_ok = true;
}

// side 20
if (Neigh20 == -1) {
  side20_ok = true;
}
else if (RedFlag[Neigh20*(level+3)+(EToRefLevel[Neigh20])-1] == 1) {
  side20_ok = true;
}
else if (pair20) {
  side20_ok = true;
}

// final decision
if (side00_ok && side10_ok && side20_ok) {
  RefFlag[e] = -1;
} else {
  RefFlag[e] = 0;
}
         
         if (Neigh00 == -1 || Neigh10 == -1 ||Neigh20 == -1 ){
         //RefFlag[e]=0;
         } 
         else if (RedFlag[Neigh00*(level+3)+(EToRefLevel[Neigh00])-1] == 1 
         && RedFlag[Neigh10*(level+3)+(EToRefLevel[Neigh10]-1)] == 1  
         && RedFlag[Neigh20*(level+3)+(EToRefLevel[Neigh20]-1)] == 1){RefFlag[e]=0;}
          
         if ((Neigh00!=-1 && Neigh01!=-1) && EToRefLevel[Neigh00]!=EToRefLevel[Neigh01] )
         {
         RefFlag[e]=0;
         } 
  
         if ((Neigh10!=-1 && Neigh11!=-1) && EToRefLevel[Neigh10]!=EToRefLevel[Neigh11])
         {
           RefFlag[e]=0;
         }

         if ((Neigh20!=-1 && Neigh21!=-1) && EToRefLevel[Neigh20]!=EToRefLevel[Neigh21])
         {
           RefFlag[e]=0;
         }
      }
      else {RefFlag[e]=0;} // Only use parent      
    } else {RefFlag[e]=0;}    
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
  printf("Flaged_Elements=%d\n",ii );
  // Coarsement Loop
  // Determine ids of new vertices and EToV

      const dlong Np = (mesh.N+1)*(mesh.N+2)/2;
      memory<dfloat> RM(4*Np*Np,0);
      o_RMRed.copyTo(RM);


  hlong nn = 0 ;
  for (int i = 0; i < ii; ++i)
  {
    int e = Coar[i];

    dlong id = e*mesh.Nverts; 
    dlong const L = EToRefLevel[e];
    if (L<=0) continue;
    dlong rep = PToC[e*(stride)+(L-1)*4+0];
    dlong sib = PToC[e*(stride)+(L-1)*4+1];
    if (rep!=-1 && sib!=-1 && RefFlag[e]==-1  && RedFlag[e*(level+3)+(EToRefLevel[e])-1]== 1  && L>0)
    {
      
      hlong sib_e2=0 ;
      hlong sib_e3=0 ;
      hlong sib_e4=0 ;
      hlong sib_id2=0;
      hlong sib_id3=0;
      hlong sib_id4=0;
      hlong v0=0 ;
      hlong v1=0 ;
      hlong v2=0 ;
      
      if (e==rep) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e2 = PToC[e*stride+(EToRefLevel[e]-1)*4+1];
         sib_e3 = PToC[e*stride+(EToRefLevel[e]-1)*4+2];
         sib_e4 = PToC[e*stride+(EToRefLevel[e]-1)*4+3];
         sib_id2 = sib_e2*mesh.Nverts;
         sib_id3 = sib_e3*mesh.Nverts;
         sib_id4 = sib_e4*mesh.Nverts;
         v0 = mesh.EToV[sib_id2+0]; 
         v1 = mesh.EToV[sib_id3+0]; 
         v2 = mesh.EToV[sib_id4+0];

      }      else // Means that child stored at added location
      { // Do nothing     
      }
      
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e2])
        {
          printf("Red Coarsened Element = %d \n",e);

          // Red Coarse
 
          EToV_new[id+0] = v0;
          EToV_new[id+1] = v1;  
          EToV_new[id+2] = v2;

          EX_new[id+0] = EX_new[sib_id2+0];
          EX_new[id+1] = EX_new[sib_id3+0];
          EX_new[id+2] = EX_new[sib_id4+0];
          
          EY_new[id+0] = EY_new[sib_id2+0];
          EY_new[id+1] = EY_new[sib_id3+0];
          EY_new[id+2] = EY_new[sib_id4+0];

          EToRefLevel[e] = EToRefLevel[e]-1;
          EToRefLevel[sib_e2] = EToRefLevel[sib_e2]-1;
          EToRefLevel[sib_e3] = EToRefLevel[sib_e3]-1;
          EToRefLevel[sib_e4] = EToRefLevel[sib_e4]-1;
          Delete_Flag[sib_e2] = 1;
          Delete_Flag[sib_e3] = 1;
          Delete_Flag[sib_e4] = 1;
          
          EToB_new[id+0] = mesh.EToB[sib_id2+0];
          EToB_new[id+1] = mesh.EToB[sib_id3+0];
          EToB_new[id+2] = mesh.EToB[sib_id4+0];

          EToB_new[sib_id2+0] = 0;
          EToB_new[sib_id2+1] = 0;
          EToB_new[sib_id2+2] = 0;

          EToB_new[sib_id3+0] = 0;
          EToB_new[sib_id3+1] = 0;
          EToB_new[sib_id3+2] = 0;

          EToB_new[sib_id4+0] = 0;
          EToB_new[sib_id4+1] = 0;
          EToB_new[sib_id4+2] = 0;      

          RefFlag[e]=0;
          RefFlag[sib_e2]=0;
          RefFlag[sib_e3]=0;
          RefFlag[sib_e4]=0;
          CombineFlag[e]=-1;
          CombineFlag[sib_e2]=-1;
          del_vertex++;
          nn = nn+3;

                for(int n=0;n<mesh.Np;++n){
        const dlong id_p = e*mesh.Np;
        Q[id_p+n] = 0.0;
      }                  

        for(int n=0;n<mesh.Np;++n){
          dfloat qn1=0.;
          dfloat qn2=0.;
          dfloat qn3=0.;
          dfloat qn4=0.;

          const dlong id_p = e*mesh.Np;
          const dlong id_c2 = sib_e2*mesh.Np;
          const dlong id_c3 = sib_e3*mesh.Np;
          const dlong id_c4 = sib_e4*mesh.Np;

            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[4*n*mesh.Np+ip+0*mesh.Np];
    
                 qn1 += Ri*Qold[id_p+ip];       
            }
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[4*n*mesh.Np+ip+1*mesh.Np];
                 qn2 += Ri*Qold[id_c2+ip];  
            } 
            for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[4*n*mesh.Np+ip+2*mesh.Np];
                 qn3 += Ri*Qold[id_c3+ip];  
            }  
           for(int ip=0;ip<mesh.Np;++ip){
            const dfloat Ri = RM[4*n*mesh.Np+ip+3*mesh.Np];
                 qn4 += Ri*Qold[id_c4+ip];  
            }     

          Q[id_p+n] = qn1+qn2+qn3+qn4;
          }                  
              CombineFlag[e]=0;
              RedFlag[e*(level+3)+(EToRefLevel[e])] = -1;
              RedFlag[sib_e2*(level+3)+(EToRefLevel[sib_e2])] = -1;
              RedFlag[sib_e3*(level+3)+(EToRefLevel[sib_e3])] = -1;
              RedFlag[sib_e4*(level+3)+(EToRefLevel[sib_e4])] = -1;
              CombineFlag[sib_e2]=0;
              CombineFlag[sib_e3]=0;
              CombineFlag[sib_e4]=0;

             if (mesh.EToE[sib_id2]!=-1)
              {
             
               
             if (RedFlag[(mesh.EToE[sib_id2]+0)*(level+3)+(EToRefLevel[mesh.EToE[sib_id2]])-1]==1)
             {
               ConfFlag[mesh.EToE[sib_id2]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id2]+0]=-1; RefFlag[mesh.EToE[sib_id2]+2]=-1; 
            printf("element %d marked element %d for green coarsening \n",e,mesh.EToE[sib_id2]);}

              }

             if (mesh.EToE[sib_id3]!=-1)
             {   
             if (RedFlag[(mesh.EToE[sib_id3]+0)*(level+3)+(EToRefLevel[mesh.EToE[sib_id3]])-1]==1)
             {
               ConfFlag[mesh.EToE[sib_id3]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id3]+0]=-1; RefFlag[mesh.EToE[sib_id3]+2]=-1;
            printf("element %d marked element %d for green coarsening \n",e,mesh.EToE[sib_id3]);}
             }

             if (mesh.EToE[sib_id4]!=-1)
             {
             if (RedFlag[(mesh.EToE[sib_id4]+0)*(level+3)+(EToRefLevel[mesh.EToE[sib_id4]])-1]==1)
             {
               ConfFlag[mesh.EToE[sib_id4]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id4]+0]=-1; RefFlag[mesh.EToE[sib_id4]+2]=-1;
            printf("element %d marked element %d for green coarsening \n",e,mesh.EToE[sib_id4]);}
             }


              printf("Neigh0=%d,Neigh1=%d,Neigh2=%d,\n",mesh.EToE[sib_id2]+0,mesh.EToE[sib_id3]+0,mesh.EToE[sib_id4]+0 );

                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
         
            }
        }
    }             
  }   

  // Loop for updating mesh information if coarsening done.
              
      if (Ncoarse!=0 && nn!=0)
      { 
        printf("inside coarsement\n"); 
        
        // Map ids of elements to new ones
        int e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
          if (Delete_Flag[e]==1)
          {
            perm[e] = -1;
          }else{
            perm[e] = e_new;
            e_new++;
          }
        }

        e_new = 0; 
        for (int e = 0; e < mesh.Nelements; ++e)
        {
                      if (Delete_Flag[e]==1)
            {
              // Do nothing
            } else{
                    if (perm[e]>= 0)
                   {
                    for (int n = 0; n < 3; ++n)
                    {
                    dlong id = e*mesh.Nverts+n;
                    dlong id_new = perm[e]*mesh.Nverts+n; 
          
                    mesh.EToB[id_new] = EToB_new[id];
                    mesh.EToV[id_new] = EToV_new[id];
                    mesh.EX[id_new]   = EX_new[id]; 
                    mesh.EY[id_new]   = EY_new[id];
                    }

                    for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    dlong id_new = perm[e]*mesh.Np+n; 
                    
                    Q_new[id_new] = Q[id];
                    
                    }
                    
                          for (int i = 0; i < stride; ++i)
                    {
                      dlong oldVal = PToC[e*stride+i];
                      dlong newVal = -1;
                    
                      if (oldVal >= 0 && oldVal < mesh.Nelements && perm[oldVal] >= 0)
                        newVal = perm[oldVal];
                    
                      PToC_new[perm[e]*stride+i] = newVal;
                    }
                  
                    for (int i = 0; i < (level+3); ++i)
                    {
                      IntFlag_new[perm[e]*(level+3)+i] = IntFlag[e*(level+3)+i];
                      RedFlag_new[perm[e]*(level+3)+i] = RedFlag[e*(level+3)+i];
                    }

                     EToRefLevel_new[perm[e]] = EToRefLevel[e];
                     RefFlag_new[perm[e]] = RefFlag[e];
                     dlong oldConf = ConfFlag[e];
                     dlong newConf = -1;

                     if (oldConf >= 0 && oldConf < mesh.Nelements && perm[oldConf] >= 0)
                       newConf = perm[oldConf];

                     ConfFlag_new[perm[e]] = newConf;   
                     
                    e_new++;  

                    }
                    }
 
            //        e_new++;
              //    }
        }


        //deviceMemory<dlong> o_combineFlag = platform.malloc<dlong>(CombineFlag);
        //o_q.copyFrom(Q);
        //deviceMemory<dfloat> o_qOld = platform.malloc<dfloat>(Q);
        //o_IntFlag = platform.malloc<dlong>(IntFlag);  
        //o_EToRefLevel = platform.malloc<dlong>(EToRefLevel); 
        //combineKernel(mesh.Nelements,Ncoarse,o_qOld ,o_q, o_combineFlag,o_IntFlag,o_EToRefLevel,o_PToC,o_RM,level);        
        
        //o_q.copyTo(Q);

        for (int e = mesh.Nelements-nn; e < mesh.Nelements; ++e)
        {
          EToRefLevel_new[e] = 0; 
          RefFlag_new[e] = 0;
          ConfFlag_new[e] = -1;
          


          for (int n = 0; n < mesh.Np; ++n)
                    {
                    dlong id = e*mesh.Np+n;
                    Q_new[id] = 0.f;
                    }

          for (int n = 0; n < (level+3); ++n)
          {
          dlong id = e*(level+3)+n;
          IntFlag_new[id] = 0;
          RedFlag_new[id] = -1;
          }   
          
          for (int n = 0; n < stride; ++n)
          {
          dlong id = e*stride+n;
          PToC_new[id] = -1;
          }   
        }
        RefFlag= RefFlag_new;
        ConfFlag = ConfFlag_new;
        RedFlag = RedFlag_new;
        PToC = PToC_new;
        EToRefLevel = EToRefLevel_new;
        IntFlag = IntFlag_new;
        Q_new.copyTo(Q); 
        o_PToC = platform.malloc<dlong>(PToC);

        mesh.Nelements = mesh.Nelements-nn;
        //mesh.Nnodes = mesh.Nnodes-del_vertex;
        mesh.o_EToB = platform.malloc<int>(mesh.EToB);
        printf("first coarsement done!!\n");
        printf("e_new=%d,Nelements=%d\n",e_new,mesh.Nelements);
        printf("Ncoarse inside coarsening loop=%lld\n",nn);
        printf("del_vertex_count=%lld\n",del_vertex);
        
        mesh = mesh.SetupUpdate(Ncoarse);
        //mesh.PmlSetup();
        //deviceMemory<dfloat> o_oldq = platform.reserve<dfloat>(mesh.Nelements*1*mesh.Np);
        //o_oldq = platform.malloc<dfloat>(Q);
        //o_q.copyFrom(o_oldq, mesh.Nelements*1*mesh.Np, 0, properties_t("async", true));
        o_q.copyFrom(Q);
      }         
}

}
