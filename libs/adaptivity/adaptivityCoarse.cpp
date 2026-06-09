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

void adaptivity_t::CoarseLocal(deviceMemory<dfloat>& o_q,
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
    if (RefFlag[e]==1 && (EToRefLevel[e]<level+3 && RedFlag[e]!=1) )
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
          ConfFlag_new[e] = 0;
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

void adaptivity_t::CoarseLocal1(deviceMemory<dfloat>& o_q,
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
    if (RefFlag[e]==-1 && (EToRefLevel[e]<level+3 && RedFlag[e]!=1) )
    {
    if (EToRefLevel[e]>0)
    {
      if (PToC[e*(level+3)+0]==e && PToC[e*(level+3)+EToRefLevel[e]]==-1)
      {
        printf("INSIDE NEW LOOP e=%d\n",PToC[2*e+0]*mesh.Np);
        printf("BEFORE PToC[e*(level+3)+0]=%d,PToC[e*(level+3)+EToRefLevel[e]]=%d, EToRefLevel[e]=%d\n", PToC[e*(level+3)+0],PToC[e*(level+3)+EToRefLevel[e]], EToRefLevel[e]);
        if(PToC[mesh.EToE[id+0]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+0]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+0]]=-1;}
        else if(PToC[mesh.EToE[id+1]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+1]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+1]]=-1;}
        else if(PToC[mesh.EToE[id+2]*(level+3)+EToRefLevel[e]]==e) {/*PToC[e*(level+3)+0]=mesh.EToE[id+2]; PToC[e*(level+3)+EToRefLevel[e]]=e; */RefFlag[mesh.EToE[id+2]]=-1;}
        printf("AFTER PToC[e*(level+3)+0]=%d,PToC[e*(level+3)+EToRefLevel[e]]=%d, EToRefLevel[e]=%d\n", PToC[e*(level+3)+0],PToC[e*(level+3)+EToRefLevel[e]], EToRefLevel[e]);
      }
      else if (PToC[e*(level+3)+0]==PToC[PToC[e*(level+3)+0]*(level+3)+0] && PToC[e*(level+3)+EToRefLevel[e]]==-1) {RefFlag[PToC[e*(level+3)+0]]=-1;}
    }
    }
  } 

    dlong ii = 0;
  for (int e = 0; e < mesh.Nelements; ++e)
  {
    if (RefFlag[e]==-1 && (EToRefLevel[e]<level+3 && RedFlag[e]!=1) )
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
          ConfFlag_new[e] = 0;
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
        ConfFlag= ConfFlag_new;
        RedFlag= RedFlag_new;
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
      if (EToRefLevel[e]>0)
      {
        printf("RED COARSENİNG BEFORE CHECK  e=%d,RedFlage=%d,parent=%d,child1=%d,child2=%d,child3=%d\n",e,RedFlag[e],PToC[e*(level+3)+0], PToC[e*(level+3)+1],PToC[e*(level+3)+2],PToC[e*(level+3)+3]);
      }
     
      // For Red ONLY!!!
      if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1 && RedFlag[e]== 1 && EToRefLevel[e]>0)
      {

        
        //printf("e_coarse=%d,e=%d,level[%d]=%d\n",PToC[4*e+0],e,PToC[4*e+0],EToRefLevel[e]);
        //printf("x1=%f,x2=%f,x3=%f\n",EX_new[id+0],EX_new[id+1],EX_new[id+2]);
        //printf("y1=%f,y2=%f,y3=%f\n",EY_new[id+0],EY_new[id+1],EY_new[id+2]);
        //hlong sib_e ;
        hlong e_parent;
        //hlong id;
        hlong sib_e2;
        hlong sib_e3;
        hlong sib_e4;
        //hlong v0 ;
        //hlong v1 ;
        //hlong v2 ;
        //hlong v0_sib;
        //hlong v1_sib;
        //hlong v2_sib;    
        //dlong swap = 0; 
        
        if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
        // Keep the Flag and check eligibility again

        // Sibling IDS
        hlong sib_e2 = PToC[e*(level+3)+1];
        hlong sib_e3 = PToC[e*(level+3)+2];
        hlong sib_e4 = PToC[e*(level+3)+3];

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

        printf("RED COARSENİNG CHECK  e=%d,Neigh00Lvl=%d,Neigh10Lvl=%d\n",e,EToRefLevel[Neigh00] ,EToRefLevel[Neigh01] );
        printf("RED COARSENİNG CHECK  e=%d,Neigh10Lvl=%d,Neigh11Lvl=%d\n",e,EToRefLevel[Neigh10] ,EToRefLevel[Neigh11] );
        printf("RED COARSENİNG CHECK  e=%d,Neigh20Lvl=%d,Neigh21Lvl=%d\n",e,EToRefLevel[Neigh20] ,EToRefLevel[Neigh21] );
        printf("RED COARSENİNG CHECK  e=%d,Neigh00=%d,Neigh10=%d\n",e,Neigh00 ,Neigh01 );
        printf("RED COARSENİNG CHECK  e=%d,Neigh10=%d,Neigh11=%d\n",e,Neigh10 ,Neigh11 );
        printf("RED COARSENİNG CHECK  e=%d,Neigh20=%d,Neigh21=%d\n",e,Neigh20 ,Neigh21 );
        printf("e=%d,x1=%f,x2=%f,x3=%f\n",e,mesh.EX[e*mesh.Nverts+0],mesh.EX[e*mesh.Nverts+1],mesh.EX[e*mesh.Nverts+2]);
        printf("e=%d,y1=%f,y2=%f,y3=%f\n",e,mesh.EY[e*mesh.Nverts+0],mesh.EY[e*mesh.Nverts+1],mesh.EY[e*mesh.Nverts+2]);
         if (Neigh00 == -1){
         //Do nothing
         }
         else if (RedFlag[Neigh00] != 1 && EToRefLevel[Neigh00]==EToRefLevel[Neigh01] && (PToC[Neigh00*(level+3)+EToRefLevel[Neigh00]]==Neigh01 || PToC[Neigh01*(level+3)+EToRefLevel[Neigh01]]==Neigh00))
         {
         if (Neigh10 == -1){
         //Do nothing
         } 
         else if (RedFlag[Neigh10] != 1 && EToRefLevel[Neigh10]==EToRefLevel[Neigh11] && (PToC[Neigh10*(level+3)+EToRefLevel[Neigh10]]==Neigh11 || PToC[Neigh11*(level+3)+EToRefLevel[Neigh11]]==Neigh10))
         {
         if (Neigh20 == -1){
         //Do nothing
         }  
         else if (RedFlag[Neigh20] != 1 && EToRefLevel[Neigh20]==EToRefLevel[Neigh21] && (PToC[Neigh20*(level+3)+EToRefLevel[Neigh20]]==Neigh21 || PToC[Neigh21*(level+3)+EToRefLevel[Neigh21]]==Neigh20))
         {
          RefFlag[e]=-1;
         }  else if(RedFlag[Neigh20] == 1 ){RefFlag[e] = -1;}  else {RefFlag[e]=0;}

         }  else if(RedFlag[Neigh10] == 1  ){RefFlag[e] = -1;}   else {RefFlag[e]=0;}
         }  else if(RedFlag[Neigh00] == 1){RefFlag[e] = -1;} else {RefFlag[e]=0;}    
         
         if (Neigh00 == -1 || Neigh10 == -1 ||Neigh20 == -1 ){
         //Do nothing
         } 
         else if (RedFlag[Neigh00] == 1 && RedFlag[Neigh10] == 1  && RedFlag[Neigh20] == 1){RefFlag[e]=0;}

         //if (Neigh00 == -1){
         if (EToRefLevel[Neigh00]!=EToRefLevel[Neigh01] )
         {
         RefFlag[e]=0;
         } 
         //}

         //if (Neigh10 == -1){
         if (EToRefLevel[Neigh10]!=EToRefLevel[Neigh11])
         {
           RefFlag[e]=0;
         }
         //}

         //if (Neigh20 == -1){
         if (EToRefLevel[Neigh20]!=EToRefLevel[Neigh21])
         {
           RefFlag[e]=0;
         }
         //}

  //       if (RedFlag[Neigh10] != 1 && EToRefLevel[Neigh10]==EToRefLevel[Neigh11] && (PToC[Neigh10*(level+3)+EToRefLevel[Neigh10]]==Neigh11 || PToC[Neigh11*(level+3)+EToRefLevel[Neigh11]]==Neigh10))
  //       {
  //        RefFlag[e]=-1;
  //       }  else if(RedFlag[Neigh10] != 1){RefFlag[e]=0;}          
//
  //       if (RedFlag[Neigh20] != 1 && EToRefLevel[Neigh20]==EToRefLevel[Neigh21] && (PToC[Neigh20*(level+3)+EToRefLevel[Neigh20]]==Neigh21 || PToC[Neigh21*(level+3)+EToRefLevel[Neigh21]]==Neigh20))
  //       {
  //        RefFlag[e]=-1;
  //       }  else if(RedFlag[Neigh20] != 1){RefFlag[e]=0;}           



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

    if (PToC[e*(level+3)+0]!=-1 && PToC[e*(level+3)+EToRefLevel[e]]!=-1 && RefFlag[e]==-1  && RedFlag[e]== 1  && EToRefLevel[e]>0)
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
      //printf("e=%d\n",PToC[2*e+0]*mesh.Np);
      
      if (e==PToC[e*(level+3)+0]) // Means that child stored at parents location
      {  
         
         // Extract Vertex Number of Elements to be Combined
         sib_e2 = PToC[e*(level+3)+1];
         sib_e3 = PToC[e*(level+3)+2];
         sib_e4 = PToC[e*(level+3)+3];
         sib_id2 = sib_e2*mesh.Nverts;
         sib_id3 = sib_e3*mesh.Nverts;
         sib_id4 = sib_e4*mesh.Nverts;
         v0 = mesh.EToV[sib_id2+0]; 
         v1 = mesh.EToV[sib_id3+0]; 
         v2 = mesh.EToV[sib_id4+0];

         //printf("PToC=%d,=%d,=%d\n", PToC[2*e+0],PToC[2*e+1],e);

      }      else // Means that child stored at added location
      { 
        printf("HEREEEEEEEEEEEEEEEEEEEEE\n");
         // Extract Vertex Number of Elements to be Combined
         //sib_e = PToC[e*(level+3)+0];
         //sib_id = sib_e*mesh.Nverts;
         //// swap id's 
         //hlong id_swap = id;
         //id = sib_id;
         //sib_id = id_swap;
         //// swap id's 
         //hlong e_swap = e;
         //e = sib_e;
         //sib_e = e_swap;
//
         //v0 = mesh.EToV[id+0]; 
         //v1 = mesh.EToV[id+1]; 
         //v2 = mesh.EToV[id+2];
         //v0_sib = mesh.EToV[sib_id+0]; 
         //v1_sib = mesh.EToV[sib_id+1]; 
         //v2_sib = mesh.EToV[sib_id+2];


     }

         //dlong id_del; // Id of the node that will be unnecessary
         //id_del = id+2;
      

         //dlong id_del2; // Id of the node that will be unnecessary (same node but for the sibling element)
         //id_del2 = sib_id+2;
         
        
        // Modify EToV with new vertex ids for bisection (3 different configurations)
        // & Calculate Physical Coordinates of new vertices
        // & Store boundary conditions of new faces
                if (EToRefLevel[e]==EToRefLevel[sib_e2])
        {

          // Red Coarse
         //if (IntFlag[e*(level+3)+EToRefLevel[e]-1]==1 )
          //{    
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
          //EToB_new[id+1] = mesh.EToB[sib_id+1];
          

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

          //const dlong id_int_p = 0;
          //const dlong id_int_c2 = 1;
          //const dlong id_int_c3 = 2;
          //const dlong id_int_c4 = 3;

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
              RedFlag[e] = 0;
              RedFlag[sib_e2] = 0;
              RedFlag[sib_e3] = 0;
              RedFlag[sib_e4] = 0;
              CombineFlag[sib_e2]=0;
              CombineFlag[sib_e3]=0;
              CombineFlag[sib_e4]=0;

             if (mesh.EToE[sib_id2]!=-1)
              {
             
               
             if (RedFlag[mesh.EToE[sib_id2]+0]==1)
             {
               ConfFlag[mesh.EToE[sib_id2]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id2]+0]=-1;}

              }

             if (mesh.EToE[sib_id3]!=-1)
             {   
             if (RedFlag[mesh.EToE[sib_id3]+0]==1)
             {
               ConfFlag[mesh.EToE[sib_id3]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id3]+0]=-1;}
             }

             if (mesh.EToE[sib_id4]!=-1)
             {
             if (RedFlag[mesh.EToE[sib_id4]+0]==1)
             {
               ConfFlag[mesh.EToE[sib_id4]+0] = e;
             } else {RefFlag[mesh.EToE[sib_id4]+0]=-1;}
             }
             //dlong  Neigh0;
             //dlong  Neigh1;
             //dlong  Neigh2;

              printf("Neigh0=%d,Neigh1=%d,Neigh2=%d,\n",mesh.EToE[sib_id2]+0,mesh.EToE[sib_id3]+0,mesh.EToE[sib_id4]+0 );
              //mesh.EToE[sib_id2]+0
              //mesh.EToE[sib_id3]+0
              //mesh.EToE[sib_id4]+0
              
             // ConfFlag[mesh.EToE[sib_id2]+0] = e;
             // ConfFlag[mesh.EToE[sib_id3]+0] = e;
             // ConfFlag[mesh.EToE[sib_id4]+0] = e;
             // ConfFlag[mesh.EToE[sib_id2]+1] = e;
             // ConfFlag[mesh.EToE[sib_id3]+1] = e;
             // ConfFlag[mesh.EToE[sib_id4]+1] = e;
             // ConfFlag[mesh.EToE[sib_id2]+2] = e;
             // ConfFlag[mesh.EToE[sib_id3]+2] = e;
             // ConfFlag[mesh.EToE[sib_id4]+2] = e;
          //printf("e_coarse0_=%d,sibe_coarse0_=%d\n",e,sib_e2);
                      for(int n=0;n<Np;++n){
              const dlong id_p = e*mesh.Np;          
              Qold[id_p+n] = Q[id_p+n];                
             //printf("Q[%d]=%f\n",id_p+n,Q[id_p+n]);
            }
          //}



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
          ConfFlag_new[e] = 0;
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
}
