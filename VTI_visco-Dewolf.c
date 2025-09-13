/************************************************************************/
/************************************************************************/
/************************************************************************/
/****** VTI-visco-De Wolf migration                    *******/

/************************************************************************/
/* First modification      Hachao Sun   $Date: 2023/05/27 22:14:15 $	*/
/* secord  modification    Hachao Sun   $Date: 2023/07/08 13:47:26 $	*/
/* third  modification     Hachao Sun   $Date: 2023/08/05 18:13:47 $	*/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/******** Credits: WTIT, Sun Huachao, shc0627@126.com,  *****************/
/************************************************************************/
/************************************************************************/

/************************************************************************/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h> 

#define PI 3.1415926
#define e 2.718281828
#include"complex.c"
#include"alloc.c"
#include"pfafft.c"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									",
" Required Parameters:							",
" nx,nz,nt	     number of samples			                ",
" fmax	sx   sz	    source samples	                                ",
" maxtrace  	    Maximum number of avenues	                        ",
" np          	    Number of threads					",
" ns ds             Number of sources 	Source interval			",
" dx,dz,dt	    horizontal sampling interval			",
" pxl,pxr        pad on left  ,pad on right				",  

NULL};


/************************ end self doc **********************************/

//Ricker 
float *ricker(float Freq,float dt,int *Npoint);

int main ()
{
   
        int nx;/* number of horizontal samples  */
        int nz;/* number of depth samples	*/
        int nt;/* number of time samples     */
        float dx;/* spatial sampling interval    */
        float dz;/* depth sampling interval 	*/
        float dt;/* time  sampling interval 	*/
        int sx;/* x source and geophone location       */
        int sz;/* z source and geophone location       */
        int snaptime,Fre;/* time of snapshot,single frequency output index 	*/
        float fmax;
        int pxl,pxr;/*pad on left  ,pad on right  */
        int ix,iz,it;
        int is,ns,ds,np;
        int iw,ik;	/* loop counters			*/
        int ntfft,nxfft;	/* fft size				*/
	int nw,nk;	        /* number of wave numbers		*/
	int ixshot;             /* dummy index		                */
        int maxtrace;           /**Total number of records ***/
	float *wl=NULL,*wtmp=NULL;	
	int ntw;
	float dw,dk;		/* wavenumber and frequency sampling interval */
	float fw,fk;		/* first wavenumber and frequency	*/
	float w,k,A,B,C,ia;		/* wavenumber and frequency		*/

	float w_out;
	float vbackground,vmin,vback,interface;/* background velocity       	*/
	float qbackground,qmin;
	double kz1,kz2,kz3,kz4;
	double phase1,phase2,phase3,phase4;

	float f1=0,f2=10,f3=70,f4=80;
	int nf1,nf2,nf3,nf4,truenw;  //Hamming  
	float w0=2.0*PI*30.0;
        float tmpp,tmppp;
        float over_relaxation_back_mo,over_relaxation_forward_mo;

        float **ve=NULL;        /*velocity model enlarge */
        float **re=NULL; 
        float *ip=NULL; 
        float **epsilu_enlage=NULL,**deta_enlage=NULL;
	complex *wlsp=NULL;

        /******Anisotropic parameter definition*********/
        float  epsilu_background,epsilu_min,deta_background,deta_min;
        double kz2_velocity,kz2_epsilu,kz2_deta;
        double velocity_disturbance,epsilu_disturbance,deta_disturbance;
	double A_k0,B_k0,phase_vti,De_aliasing;
        complex cshift2_velocity,cshift2_epsilu,cshift2_deta;

        complex velocity_disturbance_cshift1,epsilu_disturbance_cshift1,deta_disturbance_cshift1;
	complex B_k0_cshift1,phase_vti_cshift1;
        float ph1,ph2,pha,phb,a;
	double Refer_Fre,factor_qe,factor_back,forward_numerator_totle_pu,back_numerator_totle_pu;
        complex kz0_Numerator,kz0_cshift1,kz0_denominator,kz1_cshift1,kz0_totle,kz2_cshift1;
        float kz0_re,kz0_im,kz0_re1,kz0_im1;

	complex cshift1,cshift2,cshift3,cshift4,cshift5,cshift6;
        complex cshift3_dis,cshift1_dis,De_aliasing_angle;
	complex ask1,ask2,ask3,derivatives;
        complex tmp_kirchhoff;
        complex tmp_over;
        complex over_relaxation_forward;     
        complex over_relaxation_back;        
        complex over_relaxation_Kirchhoff;   

        complex Kirchhoff_numerator, Kirchhoff_denominator;
        complex over_relaxation_back_min,over_relaxation_forward_min;
        complex over_relaxation_forward1,over_relaxation_back1;
	complex forward_numerator_totle,back_numerator_totle,forward_denominator_totle,back_denominator_totle;
        complex forward_numerator,forward_denominator,back_numerator,back_denominator;  
                  

	/*Forward propagation,complex input,output*/ 
	complex **cp1=NULL;	
	complex **cq1=NULL;	  
	complex **cp2=NULL;	
	complex **cq2=NULL;	 
	complex **cp3=NULL;

        /*Back propagation,complex input,output*/
	complex **bp1=NULL;
	complex **bq1=NULL;	
	complex **bq2=NULL;
	complex **bp2=NULL;	
	complex **bp3=NULL;

        /*Back propagation,complex input,output*/
	complex **cvd1=NULL;	    complex **cvd_space1=NULL;
	complex **bvd1=NULL;	    complex **bvd_space1=NULL;

	complex **ced1=NULL;	    complex **ced_space1=NULL;
	complex **bed1=NULL;        complex **bed_space1=NULL;

	complex **cdd1=NULL;        complex **cdd_space1=NULL;
	complex **bdd1=NULL;	    complex **bdd_space1=NULL;

   
        /* get parameters */
        FILE *parfp=NULL;
        if((parfp=fopen("pars.txt","r"))==NULL){
            printf("ERROR: Cannot open parameters file:\n");
            getchar();
            exit(0);
         }
   
         fscanf(parfp,"nz=%d\n",&nz);
         fscanf(parfp,"nx=%d\n",&nx);
         fscanf(parfp,"dz=%f\n",&dz);
         fscanf(parfp,"dx=%f\n",&dx);
         fscanf(parfp,"sx=%d\n",&sx);
         fscanf(parfp,"sz=%d\n",&sz);
         fscanf(parfp,"nt=%d\n",&nt);
         fscanf(parfp,"dt=%f\n",&dt);
         fscanf(parfp,"fmax=%f\n",&fmax);
         fscanf(parfp,"pxl=%d\n",&pxl);
         fscanf(parfp,"pxr=%d\n",&pxr);
         fscanf(parfp,"maxtrace=%d\n",&maxtrace);
         fscanf(parfp,"ns=%d\n",&ns);
         fscanf(parfp,"ds=%d\n",&ds);
         fscanf(parfp,"np=%d\n",&np);

         /* printf parameters */
         printf("=======================\n");
         printf("    PARAMETERS:\n");
         printf("\tnz=%d\n",nz);
         printf("\tnx=%d\n",nx);
         printf("\tdz=%f\n",dz);
         printf("\tdx=%f\n",dx);
         printf("\tsx=%d\n",sx);
         printf("\tsz=%d\n",sz);
         printf("\tnt=%d\n",nt);
         printf("\tdt=%f\n",dt);
         printf("\tfmax=%f\n",fmax);
         printf("\tpxl=%d\n",pxl);
         printf("\tpxr=%d\n",pxr);
         printf("\tmaxtrace=%d\n",maxtrace);
         printf("\tns=%d\n",ns);
         printf("\tds=%d\n",ds);
         printf("\tnp=%d\n",np);
         printf("=======================\n");
         
         fclose(parfp);

         float **v=NULL;
         FILE *vfp=NULL;
         v=alloc2float(nx,nz);
        
         /* load velocity file */
	 if((vfp=fopen("vel.txt","r"))==NULL){
            printf("ERROR: Cannot open velocity file:\n");
            getchar();
            exit(0);
         }
         
         for(iz=0;iz<nz;iz++){
            for(ix=0;ix<nx;ix++){
                fscanf(vfp,"%f",&v[iz][ix]);
            }
         }
         fclose(vfp);

	 /* allocating space */
         float **q=NULL;
	 q=alloc2float(nx,nz);
	 /* load velocity file */
         FILE *qfp=NULL;
	 if((qfp=fopen("factor.txt","r"))==NULL){
            printf("ERROR: Cannot open velocity file:\n");
            getchar();
            exit(0);
         }
         for(iz=0;iz<nz;iz++){
            for(ix=0;ix<nx;ix++){
                fscanf(qfp,"%f",&q[iz][ix]);
            }
         }
         fclose(qfp);

         /* load record file */
         float **r=NULL;
         FILE *rfp=NULL;
         r=alloc2float(nt,maxtrace);
       
	 if((rfp=fopen("record.txt","r"))==NULL){
            printf("ERROR: Cannot open velocity file:\n");
            getchar();
            exit(0);
         }

         for(ix=0;ix<maxtrace;ix++){
             for(it=0;it<nt;it++){     
                  fscanf(rfp,"%f",&r[ix][it]);
             }
         }
         fclose(rfp); 

         float **epsilu=NULL;
         FILE *epsilufp=NULL;
         epsilu=alloc2float(nx,nz);
        
         /* load anisotropic parameter deta file */
	 if((epsilufp=fopen("deta.txt","r"))==NULL){
            printf("ERROR: Cannot open velocity file:\n");
            getchar();
            exit(0);
         }
         
         for(iz=0;iz<nz;iz++){
            for(ix=0;ix<nx;ix++){
                fscanf(epsilufp,"%f",&epsilu[iz][ix]);
            }
         }
         fclose(epsilufp);

         float **deta=NULL;
         FILE *detafp=NULL;
         deta=alloc2float(nx,nz);
        
         /* load anisotropic parameter epsilu file */
	 if((detafp=fopen("epsilu.txt","r"))==NULL){
            printf("ERROR: Cannot open velocity file:\n");
            getchar();
            exit(0);
         }
         
         for(iz=0;iz<nz;iz++){
            for(ix=0;ix<nx;ix++){
                fscanf(detafp,"%f",&deta[iz][ix]);
            }
         }
         fclose(detafp);

          time_t cur_time;
          time(&cur_time);
          printf("%s",ctime(&cur_time));

        /* enlarge the velocity model*/
         float **qe=NULL;        /*velocity model enlarge */
         qe=alloc2float(nx+pxl+pxr,nz); 

         for(iz=0;iz<nz;iz++){
             for(ix=pxl;ix<pxl+nx;ix++){
                qe[iz][ix]=q[iz][ix-pxl];
             }
             for(ix=0;ix<pxl;ix++){
                qe[iz][ix]=q[iz][0];
             }
             for(ix=pxl+nx;ix<pxl+nx+pxr;ix++){
                qe[iz][ix]=q[iz][nx-1];
             }              
          }

         /* enlarge the velocity model*/
         ve=alloc2float(nx+pxl+pxr,nz); 

         for(iz=0;iz<nz;iz++){
              for(ix=pxl;ix<pxl+nx;ix++){
                       ve[iz][ix]=v[iz][ix-pxl];
              }
              for(ix=0;ix<pxl;ix++){
                       ve[iz][ix]=v[iz][0];
              }
              for(ix=pxl+nx;ix<pxl+nx+pxr;ix++){
                       ve[iz][ix]=v[iz][nx-1];
              }              
         }

         /* enlarge the anisotropic parameter epsilu model*/
         epsilu_enlage=alloc2float(nx+pxl+pxr,nz); 

         for(iz=0;iz<nz;iz++){
              for(ix=pxl;ix<pxl+nx;ix++){
                       epsilu_enlage[iz][ix] = epsilu[iz][ix-pxl];
              }
              for(ix=0;ix<pxl;ix++){
                       epsilu_enlage[iz][ix] = epsilu[iz][0];
              }
              for(ix=pxl+nx;ix<pxl+nx+pxr;ix++){
                       epsilu_enlage[iz][ix] = epsilu[iz][nx-1];
              }              
         }

         /* enlarge the anisotropic parameter deta model*/
         deta_enlage=alloc2float(nx+pxl+pxr,nz); 

         for(iz=0;iz<nz;iz++){
              for(ix=pxl;ix<pxl+nx;ix++){
                       deta_enlage[iz][ix]=deta[iz][ix-pxl];
              }
              for(ix=0;ix<pxl;ix++){
                       deta_enlage[iz][ix]=deta[iz][0];
              }
              for(ix=pxl+nx;ix<pxl+nx+pxr;ix++){
                       deta_enlage[iz][ix]=deta[iz][nx-1];
              }              
         }

         float *window;          /*absorb window */
         window=alloc1float(nx+pxl+pxr);

         /*absorb window (hanning window)*/
         for(ix=0;ix<nx+pxl+pxr;ix++)
         {
              if(ix<=pxl){
                      window[ix]=0.5*(1+cos(PI*((pxl-ix))/pxl));
              }
              else if((ix>pxl)&&(ix<pxl+nx))
              {
                      window[ix]=1;
              }else{
                      window[ix]=0.5*(1+cos(PI*((ix-nx-pxl))/pxr));
              }
         }

         /* determine frequency sampling interval*/
         ntfft = npfar(nt);
         nw = ntfft/2+1;
         dw = 2.0*PI/(ntfft*dt);
         fw=0.0;
         printf("\tnw=%d\n",nw);
        
         /* determine wavenumber sampling (for complex to complex FFT) */
         nxfft = npfa(nx+pxl+pxr);
         nk = nxfft;
         dk = 2.0*PI/(nxfft*dx);
         fk = -PI/dx;
                     

         complex **databodyw=NULL;
         databodyw=alloc2complex(nw,maxtrace+pxl+pxr);
         for (ix=0;ix<maxtrace+pxl+pxr;++ix) {
               for(iw=0;iw<nw;++iw){
                      databodyw[ix][iw]=cmplx(0.0,0.0);
	       }
         } 

         /* enlarge the record*/
         re=alloc2float(ntfft,maxtrace); 
         for(ix=0;ix<maxtrace;ix++){
               for(it=0;it<ntfft;it++){    
                     if(it<nt){
                            re[ix][it]=r[ix][it];
                     }else{
                            re[ix][it]=0.0;
                     }
              }
        }

        /***********transform the databody from t domain to w domain (databodyw)*************/ 
        pfa2rc(1,1,ntfft,maxtrace,re[0], databodyw[0]);

        /* allocate space */
        wl=alloc1float(ntfft);
        wlsp=alloc1complex(nw);

	/* generate the Ricker wavelet */
	wtmp=ricker(fmax,dt,&ntw);

        /* zero out wl[] array */
         memset((void *) wl, 0, ntfft*sizeof(float));
		
        /* put in the wavelet in a centered fashion */ 		
        for(it=0;it<ntw;it++) {			
	       wl[(it-ntw/2+ntfft) % ntfft]=wtmp[it];	  		
        }	  	
	free1float(wtmp);

        /* fourier transform wl array */
	pfarc(-1,ntfft,wl,wlsp);
	free1float(wl);	
	                               		
	/* image */
        float **cresult_VTI_born=NULL;
	cresult_VTI_born = alloc2float(nz,nx+pxl+pxr);
        for(ix=0;ix<nx+pxl+pxr;++ix){
	     for(iz=0;iz<nz;iz++){
                  cresult_VTI_born[ix][iz]=0.0;
             }
        }

	/* image */
        float **image_VTI_born=NULL;
	image_VTI_born = alloc2float(nz,nx);
        for(ix=0;ix<nx;++ix){
	    for(iz=0;iz<nz;iz++){
                 image_VTI_born[ix][iz]=0.0;
            }
        }

        //OpenMP
	#pragma omp parallel for private(is,ix,iz,ik,iw,fw,nw,dw,fk,dk,ntfft,ixshot,nf1,nf2,nf3,nf4,truenw,forward_numerator_totle_pu,back_numerator_totle_pu,kz0_re1,kz0_im1,kz2_cshift1,De_aliasing_angle,De_aliasing,A_k0,B_k0,phase_vti,cshift2_velocity,cshift2_epsilu,cshift2_deta,epsilu_background,epsilu_min,deta_background,deta_min,kz2_velocity,kz2_epsilu,kz2_deta,velocity_disturbance,epsilu_disturbance,deta_disturbance,cshift3_dis,cshift1_dis,ask1,ask2,ask3,over_relaxation_forward_mo,over_relaxation_back_mo,over_relaxation_forward1,over_relaxation_back1,over_relaxation_back_min,over_relaxation_forward_min,forward_numerator_totle,back_numerator_totle,forward_denominator_totle,back_denominator_totle,interface,derivatives,kz1,kz2,kz3,A,B,C,ip,phase1,phase2,phase3,phase4,over_relaxation_forward,over_relaxation_back,over_relaxation_Kirchhoff,cshift1,cshift2,cshift3,cshift4,cshift5,cshift6,w,k,vbackground,vmin,tmpp,tmppp,tmp_kirchhoff,tmp_over,forward_numerator,forward_denominator,back_numerator,back_denominator,Kirchhoff_numerator, Kirchhoff_denominator,cq1,cq2,cp1,cp2,cp3,bq1,bq2,bp1,bp2,bp3,cvd_space1,ced_space1,cdd_space1,cvd1,ced1,cdd1,bvd_space1,bed_space1,bdd_space1,bvd1,bed1,bdd1,velocity_disturbance_cshift1,epsilu_disturbance_cshift1,deta_disturbance_cshift1,B_k0_cshift1,phase_vti_cshift1,ph1,ph2,pha,phb,a,Refer_Fre,factor_qe,factor_back,kz0_Numerator,kz0_cshift1,kz0_denominator,kz1_cshift1,kz0_totle,kz0_re,kz0_im,qbackground,qmin)  num_threads(np)       

        for(is=0;is<ns;is++){
              printf("\tis=%d\n",is+1);
              printf("******************************\n");

              ntfft = npfar(nt);
              nw = ntfft/2+1;
              dw = 2.0*PI/(ntfft*dt);
              fw=0.0;

              nxfft = npfa(nx+pxl+pxr);
              nk = nxfft;
              dk = 2.0*PI/(nxfft*dx);
              fk = -PI/dx;

	      /* allocate space Forward propagation */
	      cq1=alloc2complex(nxfft,nw);	
	      cp1 = alloc2complex(nx+pxl+pxr,nw);
              cp2 = alloc2complex(nx+pxl+pxr,nw);
              cp3 = alloc2complex(nx+pxl+pxr,nw);

	      /* allocate space Forward propagation about Anisotropic parameter*/
	      cvd_space1 = alloc2complex(nx+pxl+pxr,nw);
	      ced_space1 = alloc2complex(nx+pxl+pxr,nw);
	      cdd_space1 = alloc2complex(nx+pxl+pxr,nw);
	      bvd_space1 = alloc2complex(nx+pxl+pxr,nw);
	      bed_space1 = alloc2complex(nx+pxl+pxr,nw);
	      bdd_space1 = alloc2complex(nx+pxl+pxr,nw);

 	      /* allocate space back propagation */
              bq1=alloc2complex(nxfft,nw);
	      bp1 = alloc2complex(nx+pxl+pxr,nw);
	      bp2 = alloc2complex(nx+pxl+pxr,nw);
	      bp3 = alloc2complex(nx+pxl+pxr,nw);

 	      /* allocate space back propagation about Anisotropic parameter */
	      cvd1 = alloc2complex(nxfft,nw);
	      ced1 = alloc2complex(nxfft,nw);
	      cdd1 = alloc2complex(nxfft,nw);
	      bvd1 = alloc2complex(nxfft,nw);;	
	      bed1 = alloc2complex(nxfft,nw);
	      bdd1 = alloc2complex(nxfft,nw);


              /*******************************************/        
              for (ix=0; ix<nx+pxl+pxr; ix++) {
		    for (iw=0; iw<nw; iw++){     
			        cp1[iw][ix]=cmplx(0.0,0.0);
			        cp2[iw][ix]=cmplx(0.0,0.0);
			        cp3[iw][ix]=cmplx(0.0,0.0);
				bp1[iw][ix]=cmplx(0.0,0.0);
                                bp2[iw][ix]=cmplx(0.0,0.0);
                                bp3[iw][ix]=cmplx(0.0,0.0);

			        cvd_space1[iw][ix]=cmplx(0.0,0.0);
			        ced_space1[iw][ix]=cmplx(0.0,0.0);
			        cdd_space1[iw][ix]=cmplx(0.0,0.0);
			        bvd_space1[iw][ix]=cmplx(0.0,0.0);
			        bed_space1[iw][ix]=cmplx(0.0,0.0);
			        bdd_space1[iw][ix]=cmplx(0.0,0.0);

		    }
	      }
              /******************source***************/ 

	      /* compute the index of the frequency to be migrated */
	      fw=2.0*PI*f1;
	      nf1=fw/dw+0.5;
		 
	      fw=2.0*PI*f2;
	      nf2=fw/dw+0.5;

	      fw=2.0*PI*f3;
	      nf3=fw/dw+0.5;

	      fw=2.0*PI*f4;
	      nf4=fw/dw+0.5;  

	      /* the number of frequencies to migrated */
	      truenw=nf4-nf1+1;
	      fw=0.0+nf1*dw;

              nw=truenw;

                
              /* transpose the frequency domain data from     */
              /* data[ix][iw] to data[iw][ix] and apply a     */
              /* Hamming at the same time                     */
              for(ix=pxl+is*nx;ix<pxl+(is+1)*nx;++ix){
		     for (iw=0; iw<nw; iw++){
				float tmpp=0.0,tmppp=0.0;

				if(iw>=(nf1-nf1)&&iw<=(nf2-nf1)){
					tmpp=PI/(nf2-nf1);
					tmppp=tmpp*(iw-nf1)-PI;
					tmpp=0.54+0.46*cos(tmppp);
					bp1[iw][ix-is*nx]=crmul(databodyw[ix-pxl][iw+nf1],tmpp);
				} else {
					if(iw>=(nf3-nf1)&&iw<=(nf4-nf1)){
						tmpp=PI/(nf4-nf3);
						tmppp=tmpp*(iw-nf3);
						tmpp=0.54+0.46*cos(tmppp);
						bp1[iw][ix-is*nx]=crmul(databodyw[ix-pxl][iw+nf1],tmpp);
					} else {
						bp1[iw][ix-is*nx]=databodyw[ix-pxl][iw+nf1];}
		                }
				cp1[iw][ix-is*nx]=cmplx(0.0,0.0);
	               }
	        }

              /******************source***************/               
	      ixshot=is*ds+pxl+sx;		
	      for(iw=0;iw<nw;iw++){			                        
                        cp1[iw][ixshot]=wlsp[iw+nf1];
	      }              

              /********** loops over depth *********/
	      for(iz=sz;iz<nz;++iz){

	           /*****the Generalized kirchhoff imaging condition****/
                   for(ix=0;ix<nx+pxr+pxl;ix++){
	                for(iw=0;iw<nw;iw++){  
                              if(fabs(ix-ixshot)*dx<10*iz*dz){
			            tmp_kirchhoff=cmul(cp1[iw][ix],bp1[iw][ix]);
                              }else{
                                    tmp_kirchhoff=cmplx(0.0,0.0);   
                              }
                              cresult_VTI_born[ix][iz]=cresult_VTI_born[ix][iz]+tmp_kirchhoff.r/ntfft;
                        }
                  }


	          /*****get the average velocity****/
	          vmin=0.0;
                  //vmin=v[iz][0];
		  for(ix=0;ix<nx;ix++){
                        vmin=v[iz][ix]+vmin;
                        //if(vmin>v[iz][ix]) vmin=v[iz][ix];
		  }
                  vmin=vmin/nx;

	          /*****get the average Q****/
	          qmin=0.0;
		  for(ix=0;ix<nx;ix++){
			qmin=q[iz][ix]+qmin;
	          }
                  qmin=qmin/nx;

	          /*****get the average Anisotropic parameter epsilu****/
	          epsilu_min=0.0;
                  //epsilu_min=epsilu[iz][0];
		  for(ix=0;ix<nx;ix++){
			epsilu_min=epsilu[iz][ix]+epsilu_min;
                        //if(epsilu_min>epsilu[iz][ix]) epsilu_min=epsilu[iz][ix];
		  }
                  epsilu_min=epsilu_min/nx;

	          /*****get the average Anisotropic parameter deta****/
	          deta_min=0.0;
                  //deta_min=deta[iz][0];
		  for(ix=0;ix<nx;ix++){
			deta_min=deta[iz][ix]+deta_min;
                         //if(deta_min=deta[iz][ix]) deta_min=deta[iz][ix];
		  }
                  deta_min=deta_min/nx;

			
		  /* compute the shifted wavefield */
		  for (ik=0;ik<nx+pxl+pxr;++ik) {
			for (iw=0; iw<nw; ++iw) {
				cq1[iw][ik] = ik%2 ? cneg(cp1[iw][ik]) : cp1[iw][ik];
				bq1[iw][ik] = ik%2 ? cneg(bp1[iw][ik]) : bp1[iw][ik];
			}
		  }
 
		  /* zero out cq1[][] */
		  for (ik=nx+pxl+pxr; ik<nk; ++ik) { 
			for (iw=0; iw<nw; ++iw) {
				        cq1[iw][ik] = cmplx(0.0,0.0);
					bq1[iw][ik] = cmplx(0.0,0.0);
				        cvd1[iw][ik] = cmplx(0.0,0.0);
				        ced1[iw][ik] = cmplx(0.0,0.0);
				        cdd1[iw][ik] = cmplx(0.0,0.0);
				        bvd1[iw][ik] = cmplx(0.0,0.0);
				        bed1[iw][ik] = cmplx(0.0,0.0);
				        bdd1[iw][ik] = cmplx(0.0,0.0);
			}
		  }
 	

		  /* FFT to W-K domain */
		  pfa2cc(-1,1,nk,nw,cq1[0]);
		  pfa2cc(-1,1,nk,nw,bq1[0]);	

		  /**Background parameters**/
                  vbackground=vmin;
                  qbackground=qmin;
		  epsilu_background=epsilu_min;
		  deta_background=deta_min;

                  /*apply phase shift*/
		  for(ik=0,k=fk;ik<nk;++ik,k+=dk) {
			    for(iw=0,w=fw;iw<nw;++iw,w+=dw){
				 if(w==0.0) w=1.0e-10/dt; 
                                 kz1=(1.0-(((1.0+2.0*epsilu_background)*vbackground*vbackground*k*k/(w*w))/(1.0-2.0*(epsilu_background-deta_background)*vbackground*vbackground*k*k/(w*w))));
				 if(kz1>0.40){//angular range 
                                          Refer_Fre=(1.0-1.0/(qbackground*PI)*log(w/w0));
                                          kz0_cshift1=crmul(cmplx(1.0,1.0/(2.0*qbackground)),w*Refer_Fre/vbackground);
                                          kz1_cshift1=cmul(kz0_cshift1,kz0_cshift1);
                                          kz0_Numerator=csub(cmul(kz1_cshift1,kz1_cshift1),crmul(kz1_cshift1,(1.0+2.0*epsilu_background)*k*k));
                                          kz0_denominator=csub(kz1_cshift1,cmplx(2.0*(epsilu_background-deta_background)*k*k,0.0));

                                          kz2_cshift1=cdiv(kz0_Numerator,kz0_denominator);
                                         // kz2_cshift1=csqrt(cdiv(kz0_Numerator,kz0_denominator));

                                          //kz0_totle=cmul(kz0_cshift1,kz1_cshift1);
                                         // kz0_re=fabs(kz2_cshift1.r)*dz;
                                        //  kz0_im=fabs(kz2_cshift1.i)*dz;
                                         // cshift1=cmplx(cos(-kz0_re),sin(-kz0_re));
                                         // phase2=pow(e,kz0_im); 

                                         kz0_re=fabs(kz2_cshift1.r);
                                         kz0_im=fabs(kz2_cshift1.i);
                                          kz0_re1=sqrt((sqrt(kz0_re*kz0_re+kz0_im*kz0_im)+kz0_re)/2.0)*dz;
                                          kz0_im1=sqrt((sqrt(kz0_re*kz0_re+kz0_im*kz0_im)-kz0_re)/2.0)*dz;
                                          cshift1=cmplx(cos(-kz0_re1),sin(-kz0_re1));
                                          phase2=pow(e,kz0_im1); 
   
                             
                                          /******Background phase shift******/
                                          if(iz<=sz+10){
                                                 De_aliasing = w/(vbackground*sqrt(w*w/(vbackground*vbackground)+k*k*(1.6*1.6-1.0)));
                                                 if(De_aliasing>=0){
                                                        De_aliasing_angle=cmplx(De_aliasing,0.0);
                                                 }else{
                                                        De_aliasing_angle=cmplx(0.0,De_aliasing);
                                                 }
                                                 /******Background phase shift******/
                                                 cq1[iw][ik] = crmul(cmul(cq1[iw][ik],cshift1),phase2);
                                                 cq1[iw][ik] = cmul(cq1[iw][ik],De_aliasing_angle);
                                                 bq1[iw][ik] = crmul(cmul(bq1[iw][ik],cshift1),phase2);
                                                 bq1[iw][ik] = cmul(bq1[iw][ik],De_aliasing_angle);
                                           }else{
                                                 /******Background phase shift******/
                                                 cq1[iw][ik] = crmul(cmul(cq1[iw][ik],cshift1),phase2);
                                                 bq1[iw][ik] = crmul(cmul(bq1[iw][ik],cshift1),phase2);
                                           }


                                           B_k0_cshift1=cdiv(cmplx(k*k,0.0),cmul(kz0_cshift1,kz0_cshift1));
                                           /**Velocity perturbation in anisotropic parameters***/
                                           phase_vti_cshift1 = cmul(csqrt(kz2_cshift1),csub(cmplx(1.0,0.0),crmul(cmul(B_k0_cshift1,B_k0_cshift1),2.0*(epsilu_background-deta_background))));
                                           velocity_disturbance_cshift1=cdiv(cmul(kz0_cshift1,cadd(csub(cmplx(1.0,0.0),crmul(B_k0_cshift1,4.0*(epsilu_background-deta_background))),crmul(cmul(B_k0_cshift1,B_k0_cshift1),2.0*(1.0+2.0*epsilu_background)*(epsilu_background-deta_background)))),phase_vti_cshift1);
           
                                           cvd1[iw][ik]=cmul(cq1[iw][ik],velocity_disturbance_cshift1);
                                           bvd1[iw][ik]=cmul(bq1[iw][ik],velocity_disturbance_cshift1);

                                           /**epsilu perturbation in anisotropic parameters***/
                                           epsilu_disturbance_cshift1 = cdiv(crmul(cmul(kz1_cshift1,cmul(B_k0_cshift1,B_k0_cshift1)),-(1.0+2.0*deta_background)),phase_vti_cshift1);
                                           ced1[iw][ik]=cmul(cq1[iw][ik],epsilu_disturbance_cshift1);
                                           bed1[iw][ik]=cmul(bq1[iw][ik],epsilu_disturbance_cshift1);

                                           /**deta perturbation in anisotropic parameters***/
                                           deta_disturbance_cshift1 =cdiv(cmul(cmul(kz1_cshift1,B_k0_cshift1),csub(cmplx(1.0,0.0),crmul(B_k0_cshift1,(1.0+2.0*epsilu_background)))),phase_vti_cshift1);
                                           cdd1[iw][ik]=cmul(cq1[iw][ik],deta_disturbance_cshift1);
                                           bdd1[iw][ik]=cmul(bq1[iw][ik],deta_disturbance_cshift1);
				   } else {
                                                    /******Background phase shift******/
                                                    cq1[iw][ik] = cmplx(0.0,0.0);
                                                    bq1[iw][ik] = cmplx(0.0,0.0);

                                                    /**Velocity perturbation in anisotropic parameters***/
                                                    cvd1[iw][ik] = cmplx(0.0,0.0);
                                                    bvd1[iw][ik] = cmplx(0.0,0.0);

                                                    /**epsilu perturbation in anisotropic parameters***/
                                                    ced1[iw][ik] = cmplx(0.0,0.0);
                                                    bed1[iw][ik] = cmplx(0.0,0.0);

                                                    /**deta perturbation in anisotropic parameters***/
                                                    cdd1[iw][ik] = cmplx(0.0,0.0);
                                                    bdd1[iw][ik] = cmplx(0.0,0.0); 
				 }
 
		           }
		  }
				
		  /*IFFT*/
                  /******Background phase shift******/
		  pfa2cc(1,1,nk,nw,cq1[0]);
	          pfa2cc(1,1,nk,nw,bq1[0]);

                  /**Velocity perturbation***/
		  pfa2cc(1,1,nk,nw,cvd1[0]);
	          pfa2cc(1,1,nk,nw,bvd1[0]);

                  /**epsilu perturbation***/
		  pfa2cc(1,1,nk,nw,ced1[0]);
	          pfa2cc(1,1,nk,nw,bed1[0]);

                  /**deta perturbation***/
		  pfa2cc(1,1,nk,nw,cdd1[0]);
	          pfa2cc(1,1,nk,nw,bdd1[0]);

                  forward_numerator_totle = cmplx(0.0,0.0);
                  forward_denominator_totle = cmplx(0.0,0.0);
                  back_numerator_totle = cmplx(0.0,0.0);
                  back_denominator_totle = cmplx(0.0,0.0);

			
                  for(ix=0;ix<nx+pxl+pxr;++ix) {
		           for(iw=0,w=fw;iw<nw;w+=dw,++iw){
                                       if(w==0.0) w=1.0e-10/dt;
                                       /***********************************/
                                       factor_qe=(1.0-1.0/(qe[iz][ix]*PI)*log(w/w0));
                                       //factor_qe=pow((w/w0),1.0/(PI*qe[iz][ix]));
                                       factor_back=(1.0-1.0/(qbackground*PI)*log(w/w0));
                                       //factor_back=pow((w/w0),1.0/(PI*qbackground));

                                       kz2= -(1.0*factor_qe/(ve[iz][ix])-1.0*factor_back/(vbackground))*w*dz;
                                       cshift2=cmplx(cos(kz2),sin(kz2));
                                       phase3=pow(e,(1.0*factor_qe/(2.0*qe[iz][ix]*ve[iz][ix])-1.0*factor_back/(2.0*qbackground*vbackground))*w*dz);
                                        
                                       /**Velocity perturbation in anisotropic parameters***/
                                       kz2_velocity=-(1.0*factor_qe/(ve[iz][ix])-1.0*factor_back/(vbackground))*w*dz;
                                       //kz2_velocity=-w*(1.0/ve[iz][ix]-1.0/vbackground)*dz;
                                       //cshift3=cmplx(cos(-kz2_velocity)-1.0,sin(-kz2_velocity));
                                       phase4=pow(e,(1.0*factor_qe/(2.0*qe[iz][ix]*ve[iz][ix])-1.0*factor_back/(2.0*qbackground*vbackground))*w*dz);

                                       /**epsilu perturbation in anisotropic parameters***/
                                       kz2_epsilu=-(epsilu_enlage[iz][ix]-epsilu_background)*dz;
                                       //cshift2_epsilu=cmplx(cos(kz2_epsilu)-1.0,sin(kz2_epsilu));

                                       /**deta perturbation in anisotropic parameters***/
                                       kz2_deta=-(deta_enlage[iz][ix]-deta_background)*dz;
                                       //cshift2_deta=cmplx(cos(kz2_deta)-1.0,sin(kz2_deta));

                                       //Forward propagation Green
                                       /**Background field velocity correction**/
			               cq1[iw][ix]=crmul( cq1[iw][ix],1.0/nxfft);
				       cq1[iw][ix]=ix%2 ? cneg(cq1[iw][ix]) : cq1[iw][ix];
                                       cp2[iw][ix]=crmul(cmul(cq1[iw][ix],cshift2),phase3);

                                       /*Anisotropic velocity correction forward propagation**/
                                       cvd1[iw][ix]=crmul(cvd1[iw][ix],1.0/nxfft);
				       cvd1[iw][ix]=ix%2 ? cneg(cvd1[iw][ix]) : cvd1[iw][ix];
                                       cvd_space1[iw][ix]= csub(cmul(crmul(cmul(cvd1[iw][ix],cmplx(0.0,kz2_velocity)),phase4),cshift2),crmul(cmul(cp2[iw][ix],cmplx(0.0,kz2_velocity)),phase4));

                                       /*Anisotropic epsilu correction forward propagation**/
				       ced1[iw][ix]=crmul(ced1[iw][ix],1.0/nxfft);
				       ced1[iw][ix]=ix%2 ? cneg(ced1[iw][ix]) : ced1[iw][ix];
                                       ced_space1[iw][ix]=cmul(cmul(ced1[iw][ix],cmplx(0.0,kz2_epsilu)),cshift2);

                                       /*Anisotropic deta correction forward propagation**/
				       cdd1[iw][ix]=crmul(cdd1[iw][ix],1.0/nxfft);
				       cdd1[iw][ix]=ix%2 ? cneg(cdd1[iw][ix]) : cdd1[iw][ix];
                                       cdd_space1[iw][ix]=cmul(cmul(cdd1[iw][ix],cmplx(0.0,kz2_deta)),cshift2);

                                       /*forward velocity disturbance field*/ 
                                       cp1[iw][ix]=cadd(cadd(cadd(cp2[iw][ix],cvd_space1[iw][ix]),ced_space1[iw][ix]),cdd_space1[iw][ix]);
                                       cp1[iw][ix]=crmul(cp1[iw][ix],window[ix]);
                                       
                                       /**Background field velocity correction**/
			               bq1[iw][ix]=crmul( bq1[iw][ix],1.0/nxfft);
				       bq1[iw][ix]=ix%2 ? cneg(bq1[iw][ix]) : bq1[iw][ix];
                                       bp2[iw][ix]=crmul(cmul(bq1[iw][ix],cshift2),phase3);

                                       /*Anisotropic velocity correction forward propagation**/
                                       bvd1[iw][ix]=crmul(bvd1[iw][ix],1.0/nxfft);
				       bvd1[iw][ix]=ix%2 ? cneg(bvd1[iw][ix]) : bvd1[iw][ix];
                                       bvd_space1[iw][ix]= csub(cmul(crmul(cmul(bvd1[iw][ix],cmplx(0.0,kz2_velocity)),phase4),cshift2),crmul(cmul(bp2[iw][ix],cmplx(0.0,kz2_velocity)),phase4));

                                       /*Anisotropic epsilu correction forward propagation**/
				       bed1[iw][ix]=crmul(bed1[iw][ix],1.0/nxfft);
				       bed1[iw][ix]=ix%2 ? cneg(bed1[iw][ix]) : bed1[iw][ix];
                                       bed_space1[iw][ix]=cmul(cmul(bed1[iw][ix],cmplx(0.0,kz2_epsilu)),cshift2);

                                       /*Anisotropic deta correction forward propagation**/
				       bdd1[iw][ix]=crmul(bdd1[iw][ix],1.0/nxfft);
				       bdd1[iw][ix]=ix%2 ? cneg(bdd1[iw][ix]) : bdd1[iw][ix];
                                       bdd_space1[iw][ix]=cmul(cmul(bdd1[iw][ix],cmplx(0.0,kz2_deta)),cshift2);

                                       /*back velocity disturbance field*/ 
                                       bp1[iw][ix]=cadd(cadd(cadd(bp2[iw][ix],bvd_space1[iw][ix]),bed_space1[iw][ix]),bdd_space1[iw][ix]);
                                       bp1[iw][ix]=crmul(bp1[iw][ix],window[ix]);
                         }

                  } 


                 /*******************************************/            
        	 if(iz%50==0) printf("\tiz=%d\n",iz);	
	     }//iz off

             printf("******************************\n");
             free2complex(cp1);	                free2complex(cq1);
             free2complex(cp2);     

             free2complex(bp1);	                free2complex(bq1);
	     free2complex(bp2);  

             free2complex(cvd1);                free2complex(cvd_space1);     
             free2complex(bvd1);                free2complex(bvd_space1);    

             free2complex(ced1);                free2complex(ced_space1);     
             free2complex(bed1);                free2complex(bed_space1);    

             free2complex(cdd1);                free2complex(cdd_space1);  
             free2complex(bdd1);                free2complex(bdd_space1);        

        }//is off
       	  	
        /***********Storage of imaging results*********************/
         for(ix=pxl;ix<nx+pxl;ix++){
               for(iz=0;iz<nz;iz++){
                  image_VTI_born[ix-pxl][iz]=cresult_VTI_born[ix][iz];
               }
        }

              	            time_t cur_time_finish;
          time(&cur_time_finish);
          printf("%s",ctime(&cur_time_finish));

	/*****the Generalized kirchhoff imaging condition****/
        FILE *imagere_VTI_Born;	
        imagere_VTI_Born=fopen("image_re_VTI_visco_Dewolf.txt","w");

	for(iz=0;iz<nz;iz++){ 
         	for(ix=0;ix<nx;ix++){  	   	         
                        fprintf(imagere_VTI_Born,"%e ",image_VTI_born[ix][iz]);
               }            
               fprintf(imagere_VTI_Born,"\n");
        }
        fclose(imagere_VTI_Born);

       	
	/************free space********************/
        free2float(cresult_VTI_born);  	        free2float(image_VTI_born);
        free2float(epsilu_enlage);              free2float(deta_enlage); 
        free2float(epsilu);                     free2float(deta); 
        free2float(ve);                         free2float(v);
        free2float(re);                         free2float(r);  
        free2float(qe);                         free2float(q);  
        free1float(window);         
        free1complex(wlsp);  
        free2complex(databodyw);   
        return 0;
}

float * ricker(float Freq,float dt,int *Npoint)
{
	int i; /* they are the dummy counter*/
	float Bpar,t,u,*Amp;
	int Np1,N;
	
	if(Freq==0.0)Freq=30.0;
	if(dt==0.0)dt=0.004;
	Bpar=sqrt(6.0)/(PI*Freq);
	N=ceil(1.35*Bpar/dt);
	Np1=N;
	*Npoint=2*N+1;
	 
	Amp=alloc1float(*Npoint);
	
	Amp[Np1]=1.0;
  
	for(i=1;i<=N;i++) {
		t=dt*(float)i;
		u=2.0*sqrt(6.0)*t/Bpar;
		Amp[Np1+i]=Amp[Np1-i]=0.5*(2.0-u*u)*exp(-u*u/4.0);
	}

	return Amp;
}
