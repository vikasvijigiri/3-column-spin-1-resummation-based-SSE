#include "SSEvariables.hpp"
#include "SSElattice.hpp"
#include "SSEupdates.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <cstdlib>
//#include <bits/stdc++.h>

//SSEvariables var;
static ran rann;
using namespace std;

void 
SSEupdates::diagonal_update(SSElattice* lattice){
        for (int p=0; p<Lc; ++p){
            if (op_s[p] == -1){
                int b = int(rann()*Nb);
                if (lattice[bsites[b][0]].S() != lattice[bsites[b][1]].S()){
                    if (prob_in > float(Lc-n1) or prob_in > rann()*(Lc-n1)){
                        op_s[p] = 2*b;
                        n1 += 1;
                    }    
                }        
            }
            else if (op_s[p] % 2 == 0){
                    if (prob_rm*(Lc-n1+1) > 1  or prob_rm*(Lc-n1+1) > rann()){
                        op_s[p]  = -1;
                        n1 -= 1;
                    }   
            }     
            else{
                int b = int(op_s[p]/2);
                lattice[bsites[b][0]].flip();
                lattice[bsites[b][1]].flip();
            }    
        }
    }

void 
SSEupdates::offdiagonal_update(SSElattice* lattice){
 
	std::fill_n (v_first, Ns, -1); std::fill_n (v_last, Ns, -1);
	 
        
	/*
	---------------- Linked list-storage-----------
	*/
	int v1,v2;
        for (int v0=0; v0<4*Lc; ++v0){
            if (v0 % 4 != 0) continue;
            int op = op_s[int(v0/4)];
            if (op != -1){
                int b = int(op/2);
                int i = bsites[b][0];
                int j = bsites[b][1];
                v1 = v_last[i];
                v2 = v_last[j];
               
                if (v1 != -1){
                    X[v1] = v0;
                    X[v0] = v1;
                }    
                else{
                    v_first[i] = v0;
                }
                if (v2 != -1){
                    X[v2] = v0+1;
                    X[v0+1] = v2;
                }    
                else{
                    v_first[j] = v0+1;
                }    
                v_last[i] = v0+2;
                v_last[j] = v0+3;
            }     
            else{
                for (int i=v0; i<v0+4; ++i) X[i] = 0;
            }
        }
        /*
        ---------------- PBC loops --------------------------//
        */
        for (int k=0; k<Ns; ++k){
            v1 = v_first[k];
            if (v1 != -1){
                v2 = v_last[k];
                X[v2] = v1;
                X[v1] = v2;
            }
        }   
        /*
        --------------- loops to be flipped ------------------//
        */
        for (int v0=0; v0<4*Lc; ++v0){
            if (v0 % 2 != 0) continue;
            if (X[v0] < 1) continue;
            v1 = v0;
            if (rann() < 0.5){
                while(1){
                    op_s[int(v1/4)] = (op_s[int(v1/4)] ^ 1);
                    X[v1] = -1;
                    v2 = (v1 ^ 1);
                    v1 = X[v2];
                    X[v2] = -1;
                    if (v1 == v0) break;
                }    
            }        
            else{
                while(1){
                    X[v1] = 0;
                    v2 = (v1 ^ 1);
                    v1 = X[v2];
                    X[v2] = 0;
                    if (v1 == v0) break;
 		}
 	    }	
        }
        /*
        ------- updating spins after operator flips ----------//
        */
        for (int i=0; i<Ns; ++i){
            if (v_first[i] /= -1){
                if (X[v_first[i]] == -1){
                    lattice[i].flip();
                }
                v_first[i] = -1; v_last[i] = -1;
            }        
            else{
                if (rann() < 0.5){
                    lattice[i].flip();
                }
            }  
        }    
    }   

//********************************************************************************************
//********************************************************************************************
void 
SSEupdates::mcstep (SSElattice * lattice, int i){

 if (i==0) {
 	//cout << "1 imp   " << i << endl;
	update(lattice);
	//cout << "2 imp   " << i << endl;
	linkoper(lattice);
	//cout << "3 imp   " << i << endl;
	updloop(lattice);
	//cout << "4 imp   " << i << endl;
 }
/*
 else {
    	allocate_vars();
    	dupdate();
    	linkoper();
    	gupdloop();
    	deallocate_vars();    
 }
*/
}



void 
SSEupdates::adjnl(){

	 int nl1;

	 lopers=lopers/nloops;
	 nl1=1+int(float(2*Lc)/lopers);
	 nl=(nl+nl1)/2;
	 lopers=0.;
	 nloops=0.;

}

void 
SSEupdates::checkl(){

	 int j,p,dl,l1;
	 int ** tmp;
	 bool * lstr;

	 dl=int(Lc/10)+2;
	 if (n1<Lc-int(dl/2)) return;
	 l1=Lc+dl;
	 lstr = new bool[l1];
	 tmp  = new int *[2];
	 for (int i=0; i<2; ++i)tmp[i]=new int [l1];	
	  
	 for (int i=0; i<l1; ++i){
	    lstr[i]=true;
	 }
	 for (int i=0; i<dl; ++i){
	    while (true) {
	       p=int(rann()*float(l1));
	       if (lstr[p]) {
		  lstr[p]=false;
		  break;
	       }
	    }
	 }
	 j=0;
	 for (int i=0; i<l1; ++i){
	    if (lstr[i]) {
	       tmp[0][i]=str[0][j];
	       tmp[1][i]=str[1][j];
	       j=j+1;	       	       
	    }
	    else {
	       tmp[0][i]=-1;
	       tmp[1][i]=-1;	       
	    }        
	 }
	 
	     //Free each sub-array
	 for(int i = 0; i < 2; ++i) {
		delete[] str[i];   
	 }
	 //Free the array of pointers
	 delete[] str;
	 Lc=l1;
	 str= new int *[2];
	 for (int i=0; i<2; ++i)str[i]=new int [Lc];
	 	 
	 for (int i=0; i<Lc; ++i){
	    str[0][i]=tmp[0][i];
	    str[1][i]=tmp[1][i];
	 }
	 delete[] lstr;
	 for(int i = 0; i < 2; ++i) {
		delete[] tmp[i];   
	 }
	 //Free the array of pointers
	 delete[] tmp;
	 //pvect1();

}

//********************************************************************************************
void 
SSEupdates::pvect0(){


 for (int j=0; j<jmax; ++j){
    for (int s2=0; s2<2; ++s2){
       for (int s1=0; s1<2; ++s1){
          wgt[s1][s2][j]=-0.25f*jz[j]*pow(-1,1+s1)*pow(-1,1+s2)-0.5f*hz*(pow(-1,1+s1)+pow(-1,1+s2))/z;
          if (wgt[s1][s2][j]>amax[j]) amax[j]=wgt[s1][s2][j];
       }
    }
 }

 for (int j=0; j<jmax; ++j){
    amax[j]+=1.0f;
    for (int s2=0; s2<2; ++s2){
       for (int s1=0; s1<2; ++s1){
          wgt[s1][s2][j]=amax[j]-wgt[s1][s2][j];
       }
    }
 }
 cout << wgt[0][0][0] << "  " << wgt[0][1][0] << "  " << wgt[1][0][0] << "  " << wgt[1][1][0]  <<  endl;
}


void 
SSEupdates::pvect2(){

 for (int j=0; j<jmax; ++j){
   for (int s2=0; s2<2; ++s2){
     for (int s1=0; s1<2; ++s1){
       awgt[s1][s2][j]=wgt[s1][s2][j];
       if (awgt[s1][s2][j]>1e-6) {
          dwgt[s1][s2][j]=1.0f/awgt[s1][s2][j];
          ////cout << " dwgt " << dwgt[s1][s2][j] << endl;
       }   
       else {
          dwgt[s1][s2][j]=1.e6;

       }
    }
    }
 }

}


void 
SSEupdates::pvect3(){
 for (int vx=1; vx<nvx; ++vx){
   for (int ic=0; ic<4; ++ic){
     for (int oc=0; oc<4; ++oc){
        vxprb[oc][ic][vx]=vxp[oc][ic][vx];
     }
   }
 }
}



void 
SSEupdates::initvrtx_dirloop(){


 int v0,v1,v2,v3,v4,v5,v6,v7,v8;
 int j, k, vx, ic, oc, vxn;
 int st[4],vxtyp[nvx]={0};
 double vxwgt[nvx]={0.},legwgt[4][4][nvx]={0.};

 int    jnn,jn[4],jo[4],jv[4];
 double jwn,fac,jw[4],mwgt[4][4];

 int jdog[5]={0,1,0,2,1},kdog[5]= {3,2,1,3,2}, icn, nj, nk, vxk;
 

 //Define legvx, vxoper, vxwgt, and vxleg
 for (int j=0; j<jmax; ++j){
    v0=8*j;
    v1=v0+1;
    v2=v0+2;
    v3=v0+3;
    v4=v0+4;
    v5=v0+5;
    v6=v0+6;
    v7=v0+7;
    v8=v0+8;

    legvx[1][0][0][1][j]=v1;
    legvx[0][1][1][0][j]=v2;
    legvx[1][0][1][0][j]=v3;
    legvx[0][1][0][1][j]=v4;
    legvx[1][1][1][1][j]=v5;
    legvx[0][0][0][0][j]=v6;
    legvx[1][1][0][0][j]=v7;
    legvx[0][0][1][1][j]=v8;

    for (int k=1; k<nvx; ++k){
       vxtyp[v0+k]=j;
    }

    vxoper[v1]=2;
    vxoper[v2]=2;
    vxoper[v3]=1;
    vxoper[v4]=1;
    vxoper[v5]=1;
    vxoper[v6]=1;
    vxoper[v7]=2;
    vxoper[v8]=2;

    vxwgt[v1]=jx[j];
    vxwgt[v2]=jx[j];
    vxwgt[v3]=wgt[1][0][j];
    vxwgt[v4]=wgt[0][1][j];
    vxwgt[v5]=wgt[1][1][j];
    vxwgt[v6]=wgt[0][0][j];
    vxwgt[v7]=jy[j];
    vxwgt[v8]=jy[j];

    for (int t2=0; t2<2; ++t2){
       for (int t1=0; t1<2; ++t1){
          for (int s2=0; s2<2; ++s2){
             for (int s1=0; s1<2; ++s1){
                vx=legvx[s1][s2][t1][t2][j];
                if (vx != 0) {
                   vxleg[0][vx]=s1;
                   vxleg[1][vx]=s2;
                   vxleg[2][vx]=t1;
                   vxleg[3][vx]=t2;
                }
             }
          }
       }
    }

 }

 //legwgt[4][4][nvx]={0.}; 	// legwgt(4,4,nvx)=0.d0 //
 //vxnew[4][4][nvx]={0};	// vxnew(4,4,nvx)=0     //
 //vxp[4][4][nvx]={0.}; 	// vxp(4,4,nvx)=0.d0    //

 //Create vxnew links
 for (int vx=1; vx<nvx; ++vx){
    j=vxtyp[vx];
    for (ic=0; ic<4; ++ic){
       for (oc=0; oc<4; ++oc){
          st[0]=vxleg[0][vx];
          st[1]=vxleg[1][vx];
          st[2]=vxleg[2][vx];
          st[3]=vxleg[3][vx];
          st[ic]=1-st[ic];
          st[oc]=1-st[oc];
          vxn=legvx[st[0]][st[1]][st[2]][st[3]][j];
          ////cout << "j value is " << vxn << endl;//pow(-1,1+st[0]) << " " << pow(-1,1+st[1]) << " " << pow(-1,1+st[2]) << " " << pow(-1,1+st[3]) << endl;          
          if (vxn != 0) {
             vxnew[oc][ic][vx]=vxn;
          }
       }
    }
 }

 //Choose a vertex and in channel
 for (int vx=1; vx<nvx; ++vx){
    for (int ic=0; ic<4; ++ic){
       //Now we have list of vertices in Directed Loop equations
       //Initialize vertex order
       j=0;
       for (oc=0; oc<4; ++oc){ 
          vxn=vxnew[oc][ic][vx];

          if (vxn != 0) {
             jn[j]=j;
             jo[j]=oc;
             jv[j]=vxn;
             jw[j]=vxwgt[vxn];
             j=j+1; 
                         
          }
       }
       if (j != 4) {        
        printf("Error in values of vxnew");
  
        // EXIT_FAILURE
        exit(1);}
       //Five moves to absolute order
       for (oc=0; oc<5; ++oc){ 
          j=jdog[oc];
          k=kdog[oc];
          if (jw[j]<jw[k]) {
             jnn=jn[j];
             jn[j]=jn[k];
             jn[k]=jnn;
             jwn=jw[j];
             jw[j]=jw[k];
             jw[k]=jwn;
          }
       }
       //{ implement single-bounce or bounce-free (see Eq. 25 of Sylju\r{a}sen)
       mwgt[0][0]=std::max(0.,jw[0]-jw[1]-jw[2]-jw[3]);
       mwgt[0][1]=std::min(jw[1],0.5*(jw[0]+jw[1]-jw[2]-jw[3]));
       mwgt[0][2]=std::min(jw[2],0.5*(jw[0]-jw[1]+jw[2]-jw[3]));
       mwgt[0][3]=jw[3];
       mwgt[1][2]=std::max(0.,0.5*(-jw[0]+jw[1]+jw[2]+jw[3]));
       mwgt[1][3]=0.f;
       mwgt[2][3]=0.f;
       ////cout << " " << jw[0] << " " << jw[1] << " " << jw[2] << " " << jw[3] << endl;
       for (int j=1; j<4; ++j){
          mwgt[j][j]=0.f;
          for (int k=0; k<j; ++k){
             mwgt[j][k]=mwgt[k][j];
          }
       }
       //Finally, assign weights to Directed Loop vertex elements
       for (int j=0; j<4; ++j){
          nj=jn[j];
          icn=jo[nj];
          vxn=jv[nj];
          for (int k=0; k<4; ++k){
             nk=jn[k];
             vxk=jv[nk];
             // Cycle through to find correct out channel
             // (the one that changes vxj to vxk)
             for (oc=0; oc<4; ++oc){ 
                if (vxnew[oc][icn][vxn]==vxk) break;
                if (oc==3) {        
		printf("Error in values of vxk");
	  
		// EXIT_FAILURE
		exit(1);} 
             }
             legwgt[oc][icn][vxn]=mwgt[j][k];
          }
       }
    }

 }

 //Convert matrix weights to probabilities
 for (int vx=1; vx<nvx; ++vx){
    fac=1.f/vxwgt[vx];
    for (int ic=0; ic<4; ++ic){
       for (int oc=0; oc<4; ++oc){ 
          vxp[oc][ic][vx]=fac*legwgt[oc][ic][vx];
          ////cout << "vxp " << vxp[oc][ic][vx] << endl;
       }
    }
 }

 //Create cumulative probabilities
 for (int vx=1; vx<nvx; ++vx){
    for (int ic=0; ic<4; ++ic){
       for (int oc=1; oc<4; ++oc){ 
          vxp[oc][ic][vx]+=vxp[oc-1][ic][vx];
       }
    }
 }

 //Normalize cumulative probabilities
 for (int vx=1; vx<nvx; ++vx){
    for (int ic=0; ic<4; ++ic){
       for (int oc=0; oc<4; ++oc){ 
          if (vxp[3][ic][vx]<1e-6) {
              vxp[oc][ic][vx]=0.f;
          }   
          else {
             vxp[oc][ic][vx]=vxp[oc][ic][vx]/vxp[3][ic][vx];
             if (vxp[oc][ic][vx]<1e-6) vxp[oc][ic][vx]=-1.f;
          }
       }
    }
 }


}





void 
SSEupdates::update(SSElattice * lattice){

 int b, o, ss1, ss2, typ;
 double p;

 plnk = new int * [4];
 slnk = new int * [4];
 bnd  = new int [Lc];    
 lpos = new int [Lc];
 lvtx = new int [Lc];    
 for (int i = 0; i < 4; i++) {plnk[i] = new int[Lc]; slnk[i] = new int[Lc];}
 
 
 
 for (int i=0; i<Lc; ++i){
    o=str[0][i];
    if (o==-1) {
       b = int(rann()*Nb);
       typ=btyp[b];
       ss1=(1+lattice[bsites[b][0]].S())/2;
       ss2=(1+lattice[bsites[b][1]].S())/2;
       p=awgt[ss1][ss2][typ]*(Beta*float(Nb)/float(Lc-n1));
       if (p>1.0f or rann()<p) {
          str[0][i]=1;
          str[1][i]=b;
          n1+=1;
       }
    }   
    else if (o==1) {   
       b=str[1][i];
       typ=btyp[b];
       ss1=(1+lattice[bsites[b][0]].S())/2;
       ss2=(1+lattice[bsites[b][1]].S())/2;
       p=dwgt[ss1][ss2][typ]*(float(Lc-n1+1)/(Beta*float(Nb)));
       if (p>1.0f or rann()<p) {
          str[0][i]=-1;
          str[1][i]=-1;
          n1-=1;
       }
    }   
    else {
       b=str[1][i];
       lattice[bsites[b][0]].flip();
       lattice[bsites[b][1]].flip();
    }
 }

}
    
    
void 
SSEupdates::linkoper(SSElattice * lattice){


 int  b, o, s1, s2, p1, p2, ss1, ss2, tt1, tt2, typ;
 int **lspn;
 
 lspn = new int * [2];
 for (int i = 0; i < 2; i++) lspn[i] = new int[Lc];
 
 for (int i=0; i<Ns; ++i){
    frst[i]=0;
    last[i]=0;
 }
 no=0;
 for (int i=0; i<Lc; ++i){
    o=str[0][i];
    if (o != -1) {
       b=str[1][i];
       bnd[no]=b;
       typ=btyp[b];
       s1=bsites[b][0];
       s2=bsites[b][1];
       ss1=lattice[s1].S();
       ss2=lattice[s2].S();
       if (o == 2) {
          lattice[s1].flip();
          lattice[s2].flip();
       }
       tt1=lattice[s1].S();
       tt2=lattice[s2].S();
       lpos[no]=i;
       lspn[0][no]=s1;
       lspn[1][no]=s2;
       lvtx[no]=legvx[(1+ss1)/2][(1+ss2)/2][(1+tt1)/2][(1+tt2)/2][typ];
       if (lvtx[no]==0){printf("oops spins not correctly converted into 0 or 1 or vice-versa!"); exit(1);}
       p1=last[s1];
       p2=last[s2];
       if (p1 != 0) {
          if (lspn[0][p1] == s1) {
             plnk[2][p1]=no;
             slnk[2][p1]=0;
             plnk[0][no]=p1;
             slnk[0][no]=2;
          }   
          else{
             plnk[3][p1]=no;
             slnk[3][p1]=0;
             plnk[0][no]=p1;
             slnk[0][no]=3;
          }
       }   
       else {
          frst[s1]=no;
          fspn[s1]=0;
       }
       if (p2 != 0) {
          if (lspn[0][p2] == s2) {
             plnk[2][p2]=no;
             slnk[2][p2]=1;
             plnk[1][no]=p2;
             slnk[1][no]=2;
          }
          else {
             plnk[3][p2]=no;
             slnk[3][p2]=1;
             plnk[1][no]=p2;
             slnk[1][no]=3;
          }
       }   
       else {
          frst[s2]=no;
          fspn[s2]=1;
       }
       last[s1]=no;
       last[s2]=no;
       no+=1;       
    }
 }

 for (int s1=0; s1<Ns; ++s1){
    p1=last[s1];
    if (p1 != 0) {
       p2=frst[s1];
       if (s1 == lspn[0][p1]) {
          plnk[2][p1]=p2;
          if (s1==lspn[0][p2]) {
             slnk[2][p1]=0;
          }
          else {
             slnk[2][p1]=1;
          }
       }   
       else {
          plnk[3][p1]=p2;
          if (s1==lspn[0][p2]) {
             slnk[3][p1]=0;
          }
          else{
             slnk[3][p1]=1;
          }
       }
       if (s1==lspn[0][p2]) {
          plnk[0][p2]=p1;
          if (s1==lspn[0][p1]) {
             slnk[0][p2]=2;
          }
          else {
             slnk[0][p2]=3;
          }
       }   
       else {
          plnk[1][p2]=p1;
          if (s1==lspn[0][p1]) {
             slnk[1][p2]=2;
          }
          else {
             slnk[1][p2]=3;
          }
       }
    }
 }
 for(int i = 0; i < 2; ++i)
 {
    delete[] lspn[i];
 }
}


    
void 
SSEupdates::updloop(SSElattice * lattice){

 int  ic, oc, nop, nop1, os, p, vx;
 double r;
 //cout << "hererereree   1" << endl; 
 nop=0; 
 for (int j=0; j<nl; ++j){
    ic=int(rann()*4);
    p =std::min(int(rann()*no)+1,no-1);
    nop1=0;
    while (1){
       nop1+=1;
       r=rann();
       vx=lvtx[p];
              //cout << "inside  1 " << p << "   " << nop1 << "   " << vx  << "    " << no << endl;
       for (oc=0; oc<4; ++oc){  
          if (r < vxprb[oc][ic][vx]) {
             lvtx[p]=vxnew[oc][ic][vx];
             break;
          }

          if (oc==3) { printf("some problem in vextex probs!");  exit(1);}
       }
       //cout << "inside  2 "<< oc << "     " << p << "        "<< lvtx[p] << " " << slnk[oc][p] << "    " << plnk[oc][p]<< endl;       
       ic=slnk[oc][p];
       os=vxleg[oc][lvtx[p]];
       p=plnk[oc][p];
       //cout << "inside  3 " << ic << "     " << p << "    " << os << endl;
       if (os == vxleg[ic][lvtx[p]]){ break;}
    }
    nop+=nop1;
           //cout << "loop    " << j<< endl;       
 }
 //cout << "hererereree   2" << endl;          
 for (int i=0; i<no; ++i){
    str[0][lpos[i]]=vxoper[lvtx[i]];
 }
 //cout << "hererereree   3" << endl; 
 for (int i=0; i<Ns; ++i){
    if (frst[i] != 0) {
       lattice[i].set_S(pow(-1,1+vxleg[fspn[i]][lvtx[frst[i]]]));
       ////cout << pow(-1,1+vxleg[fspn[i]][lvtx[frst[i]]]) << endl;
    }
 }
 lopers+=float(nop);
 nloops+=float(nl);
 //cout << "hererereree   4" << endl; 
 //*****************************************
 for(int i = 0; i < 4; ++i)
 {
	delete[] plnk[i];
	delete[] slnk[i]; 
		
 } 
 delete[] bnd;
 delete[] lpos;
 delete[] lvtx;
  //cout << "hererereree   5" << endl; 
}
    
     
               

void 
SSEupdates::adjust_cutoff_length(){
        int Lc1 = int(n1+n1/3);

        if (Lc1 > Lc){

        int* op_singcopy = new int [Lc]; 
        for (int i=0; i<Lc; ++i) op_singcopy[i] = op_s[i];
        delete[] op_s;
        op_s = new int [Lc1]; 
        for (int i=0; i<Lc; ++i) {op_s[i] = op_singcopy[i];}
        for (int i=Lc; i<Lc1; ++i) op_s[i]=-1;
        
        delete[] op_singcopy;

        Lc = Lc1;
        delete[] X; 
        X = new int [4*Lc]; 
        }
        }




void 
SSEupdates::initialize(){
     	Lc=88;
 	n1=0;
 	nl=5;
 	/*
        op_s = new int [Lc];
        std::fill_n (op_s, Lc, -1);
        v_first = new int [Ns];         	v_last  = new int [Ns];
        std::fill_n (v_first, Ns, -1); 	std::fill_n (v_last, Ns, -1);      
        X = new int [4*Lc];
        std::fill_n (X, 4*Lc, -1);
	*/
        str = new int * [2];
 	for (int i = 0; i < 2; i++) str[i] = new int[Lc];
        
        for (int i = 0; i < 2; i++){
        for (int j = 0; j < Lc; j++){
        str[i][j]=-1;
        }}
        
 	frst = new int [Ns]; 
 	last = new int [Ns]; 
 	fspn = new int [Ns];
        }


