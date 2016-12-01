////////////////////////////////////////////////////////
////*********Transparent boundary condition*********////
////////////////////////////////////////////////////////
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<omp.h>
#include <string.h>
/******************************************************/
double randgen()
{
double r, R, t, b;
r=rand();
R=RAND_MAX;
t=r/R;
if (t>=0.2)
	{
	b=0;
	}
	else
	{
	b=1;
	}
return(b);
}
/*****************************************************/
main()
{
srand ((unsigned)time(NULL));
char buffer[1];
long int   a,jj,ii,gg,ff,mm,Nx,Ny,Nz,Ncw,nc;
double x0, y0, pi, lam, k0, w0, w, wc, dx, dy, dc, n0, dz, na, ng, rnd, Nsites, Px=0, Py=0, Ps=0, N=0, nn=1, R, theta;
/*****************************************************/
lam = 1;     /*wavelength*/
w0 = 3;    	/* Gaussian pulse width */
w = 100; 
dc = 0.5*lam;
R = 13000;
theta = 1;
Nsites = w/dc;  
Nsites = (int) Nsites; 
dx = 0.5*lam;   /*step size in X direction*/
nc = dc/dx;
nc = (int)nc;
dy = dx;
w = Nsites*nc*dx;	/*simulation window width*/
Nx = (int)(w/dx);      /*n number of samples in "x" direction*/
Ny = Nx;
wc = 4.1; 
Ncw = (int)(wc/dx);
double (*n)[Nx] = malloc(Nx*sizeof(double[Nx]));
double (*nS)[Nx] = malloc(Nx*sizeof(double[Nx]));
double (*nB)[Nx] = malloc(Nx*sizeof(double[Nx]));
double complex (*k1)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex (*k2)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex (*k3)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex (*k4)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex kx1, ky1;
double (*E0)[Nx] = malloc(Nx*sizeof(double[Nx]));
double complex (*E1)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex (*E2)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
/*******************************************************/
na = 1;        /*refractive index of Air */
ng = 1.5;        /*refractive index of core*/
/*define  the incident boundary condition */
x0 = (Nx/2)*dx;  /*center of Gaussian pulse*/
y0 = (Ny/2)*dx;
//printf("real center=%f\n",x0); 
#pragma omp parallel for
for(ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for   
    for (jj=0;jj<=Ny-1;jj++)
               E0[ii][jj]=exp(-(pow((ii+1)*dx-x0,2)+pow((jj+1)*dy-y0,2))/pow(w0,2));
/********************************************************/
for (gg=0;gg<Nsites;gg++)
	{
		for (ff=0;ff<Nsites;ff++)
			{
		rnd = randgen();
		for (ii=0;ii<=nc-1;ii++)
			     {   
				   for(jj=0;jj<=nc-1;jj++)   
				    	{
			        	        nS[ii+nc*gg][jj+nc*ff]=ng+rnd*(na-ng);
				    	}
			     }
			}
	}
/*********************************************************/
//for (ii=0;ii<=Nx-1;ii++){
//        for (jj=0;jj<=Ny-1;jj++){
//                if ((pow(ii-Nx/2,2)+pow(jj-Ny/2,2))>pow(Ncw,2)){
//                        nS[ii][jj]=na;}
//
//                }
//        }
/**********************************************************/
for (ii=0;ii<=Nx-1;ii++)
	{
		for(jj=0;jj<=Ny-1;jj++)
			{
				nB[ii][jj]=nS[ii][jj]*exp((ii*dx-x0)/R); 
			}
	}
/*********************************************************/
pi = 4*atanl(1);		/*pi*/
k0 = 2*pi/lam; 			/*wavevector*/
n0 = 0.5*(na+ng);	 	/*effective refractive index of media*/
dz = 0.01*k0*pow(dx,2)*n0;
printf("dz=%f\n",dz);	
Nz = (int)(5000/dz);        	/*number of samples in "z" direction */ 
double *SigmaS = malloc(Nz * sizeof (double));
double *Sigmax = malloc(Nz * sizeof (double));
double *Sigmay = malloc(Nz * sizeof (double));
double *power = malloc(Nz * sizeof (double));
printf("Nx = %ld, dc = %f, dx = %f, dz =%f,  Nsites = %f, w = %f\n", Nx, dc, dx, dz, Nsites, w);
/***********************************************************/
#pragma omp parallel for 
for (ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for 
    for(jj=0;jj<=Ny-1;jj++)
     E1[ii][jj]=E0[ii][jj];
/***********************************************************/
//#pragma omp parallel for 
//for (ii=0;ii<=Nx-1;ii++){
//#pragma omp parallel for 
//    for(jj=0;jj<=Ny-1;jj++){
//     k1[ii][jj]=0;
//     k2[ii][jj]=0;
//     k3[ii][jj]=0;
//     k4[ii][jj]=0;}}
/****************************************************************/
/********************Runge-Kutta method**************************/
/****************************************************************/
for (mm=1;mm<=Nz;mm++)
    {
//if (mm*dz>200 && mm*dz<(200+R*theta))
//	{	
		n=nB;
//	}
//	else 
//	{
//		n=nB; 	
//	}
Px=0;
Py=0; 
Ps=0; 
N=0; 
power[mm-1]=0;
SigmaS[mm-1]=0;
Sigmax[mm-1]=0;
Sigmay[mm-1]=0;
#pragma omp parallel for 
     for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
      for (jj=1;jj<=Ny-2;jj++)
                  k1[ii][jj] = 
(-I*nn/n0/k0/2)*((E1[ii+1][jj]-2*E1[ii][jj]+E1[ii-1][jj])/pow(dx,2)+(E1[ii][jj+1]-2*E1[ii][jj]+E1[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E1[ii][jj]);

#pragma omp parallel for 
    for(ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for 
        for(jj=0;jj<=Ny-1;jj++)
            E2[ii][jj]=E1[ii][jj]+0.5*dz*k1[ii][jj]; 
/***********************************************/
#pragma omp parallel for
       for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for
        for(jj=1;jj<=Ny-2;jj++)
        k2[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);
#pragma omp parallel for
        for (ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for
                  for(jj=0;jj<=Ny-1;jj++)      
                        E2[ii][jj]=E1[ii][jj]+0.5*dz*k2[ii][jj];
/***********************************************/
#pragma omp parallel for
       for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
        for (jj=1;jj<=Ny-2;jj++)
                k3[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);      
#pragma omp parallel for
                for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
             for (jj=1;jj<=Ny-2;jj++)
                                E2[ii][jj]=E1[ii][jj]+ dz*k3[ii][jj]; 
/***********************************************/
#pragma omp parallel for
for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for
    for (jj=1;jj<=Ny-2;jj++)
                k4[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);
#pragma omp parallel for
                for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for
            for (jj=1;jj<=Ny-2;jj++)                                
E2[ii][jj]=E1[ii][jj]+(k1[ii][jj]+2*k2[ii][jj]+2*k3[ii][jj]+k4[ii][jj])*dz/6;
/*************TBC Upper**************************/
//#pragma omp parallel for
for (jj=0;jj<=Ny-1;jj++)
        {   
        if (E2[2][jj]!=0){
        kx1=I/dx*clog(E2[1][jj]/E2[2][jj]);
        if (creal(kx1)<0)
            { kx1=0;}
                E2[0][jj]=E2[1][jj]*cexpl(-I*kx1*dx);
        }}                  
/*************TBC Lower*************************/
//#pragma omp parallel for
for (jj=0;jj<=Ny-1;jj++)
        {   
        if (E2[Nx-3][jj]!=0){
        kx1=I/dx*clog(E2[Nx-2][jj]/E2[Nx-3][jj]);
                if (creal(kx1)<0)
                        { kx1=0;}
                E2[Nx-1][jj]=E2[Nx-2][jj]*cexpl(-I*kx1*dx);
        }}
/*************TBC Left**************************/
//#pragma omp parallel for
for (ii=0;ii<=Nx-1;ii++)
        {
        if (E2[ii][2]!=0){   
        ky1=I/dy*clog(E2[ii][1]/E2[ii][2]);
                if (creal(ky1)<0)
                        { ky1=0;}
                E2[ii][0]=E2[ii][1]*cexpl(-I*ky1*dx);
        }}
/*************TBC Right**************************/
//#pragma omp parallel for
for (ii=0;ii<=Nx-1;ii++)
        {
        if(E2[ii][Ny-3]!=0){   
        ky1=I/dy*clog(E2[ii][Ny-2]/E2[ii][Ny-3]);
                if (creal(ky1)<0)
                        { ky1=0;}
                E2[ii][Ny-1]=E2[ii][Ny-2]*cexpl(-I*ky1*dx);
        }}
/***********************************************/
#pragma omp parallel for
        for (ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for
            for (jj=0;jj<=Ny-1;jj++)
                E1[ii][jj]=E2[ii][jj];
/***********************************************/
//#pragma omp parallel for
       for (ii=0;ii<=Nx-2;ii++ )
//#pragma omp parallel for
            for (jj=0;jj<=Ny-2;jj++)
			N = N +  pow(cabs(E2[ii][jj]),2) * pow(dx,2);
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
			Ps = Ps + pow(cabs(E2[ii][jj]),1) * dx * dx;
//#pragma omp parallel for
      //  for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
        //   for (jj=0;jj<=Ny-2;jj++ )
                        SigmaS[mm-1] = pow(N/pow(Ps,2), -0.5);
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
                        power[mm-1] = power[mm-1]  +  pow(cabs(E2[ii][jj]),2) * dx * dx;
/******************************************************/
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
			Px = Px  + ii * dx/N * pow(cabs(E2[ii][jj]),2) * dx * dx;
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
                        Sigmax[mm-1] = Sigmax[mm-1]  + pow((ii*dx-Px),2) 
* pow(cabs(E2[ii][jj]),2) /N* dx * dx;
/*******************************************************/
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
			Py = Py + jj * dy/N * pow(cabs(E2[ii][jj]),2) * dx * dy;
//#pragma omp parallel for
        for (ii=0;ii<=Nx-2;ii++)
//#pragma omp parallel for
           for (jj=0;jj<=Ny-2;jj++ )
                        Sigmay[mm-1] = Sigmay[mm-1] + pow((jj*dy-Py),2) 
* pow(cabs(E2[ii][jj]),2) /N* dx * dy;
/*******************************************************/
//printf("Step=%ld, Center=%f, Sigmax=%f, Sigmay=%f, SigmaS=%f, Propagation distance(mum)=%f, power=%f\n", mm, Ps, 2*sqrt(Sigmax[mm-1]), 2*sqrt(Sigmay[mm-1]), SigmaS[mm-1], mm*dz,power[mm-1]);
if (mm==1 || mm%500==0){
	sprintf(buffer, "%ld", mm);
	strcat(buffer,".dat");	 
        FILE *f;
        f = fopen(buffer,"w");
        for (ii=0;ii<=Nx-1;ii++)
            {
             fprintf(f,"\n");   
             for(jj=0;jj<=Ny-1;jj++)
            {
                fprintf(f,"%f", cabs(E2[ii][jj]));
                fprintf(f,"  ");
            }
            }   
                fclose(f);
	    }
FILE *fff;
fff = fopen("counter.dat", "w");
fprintf(fff,"%f", mm*dz);
fprintf(fff,"\n");
fclose(fff); 	
}
/***********Absolute value of the output********/       
//        FILE *f;
//        f = fopen("N5E.dat","w");
//        for (ii=0;ii<=Nx-1;ii++)
//            {
//             fprintf(f,"\n");   
//            for(jj=0;jj<=Ny-1;jj++)
//            {
//                fprintf(f,"%f", cabs(E2[ii][jj]));
//                fprintf(f,"  ");
//            }
//            }   
//                fclose(f);
/*************************************************/
      FILE *g1;
        g1 = fopen("SigmaS.dat","w");
        for (mm=1;mm<=Nz;mm++)
            {
             fprintf(g1,"%f  %f", mm*dz, SigmaS[mm-1]);
             fprintf(g1,"\n");   
            }   
                fclose(g1);
/**************************************************/
      FILE *g2;
        g2 = fopen("Sigmax.dat","w");
        for (mm=1;mm<=Nz;mm++)
            {
             fprintf(g2,"%f  %f", mm*dz, 2*sqrt(Sigmax[mm-1]));
             fprintf(g2,"\n");   
            }   
                fclose(g2);
/**************************************************/
      FILE *g3;
        g3 = fopen("Sigmay.dat","w");
        for (mm=1;mm<=Nz;mm++)
            {
             fprintf(g3,"%f  %f", mm*dz, 2*sqrt(Sigmay[mm-1]));
             fprintf(g3,"\n");   
            }   
                fclose(g3);
/**************************************************/
      FILE *H;
        H = fopen("power.dat","w");
        for (mm=1;mm<=Nz;mm++)
            {
             fprintf(H,"%f  %f", mm*dz, power[mm-1]);
             fprintf(H,"\n");   
            }   
                fclose(H);
}








