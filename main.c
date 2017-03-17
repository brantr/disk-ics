#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "disk_ics.h"


int main(int argc, char  **argv)
{
	DiskICs d;

	int i;
	int nz = 128;	   //grid size in z positive direction
	int nzout = 128;	   //grid size in z positive direction
	//int nz = 2560;
	//int nz = 12800;	   //grid size in z positive direction
	//int nz=128;
	double z_min = 0;  //bottom of z cell column in kpc
	double z_max = 10; //top of z cell column in kpc
	double dz;         // cell width
	double dzout;         // cell width

	double *z;			  // array of cell centers
	double *zout;			  // array of cell centers

	double *rho;		  // array of cell averaged densities
	double *rhog;		  // array of hydrostatic weight
	double *p; 		    // array of pressure
	double *gradp;		// array of pressure gradients
	double *gza;		  // array of accelerations
	double *drhoa;		  // array of accelerations
	double *rhoout;		  // array of cell averaged densities
	double *pout;		  // array of cell averaged densities
	double *gzaout;		  // array of cell averaged densities
	double *gradpout;		  // array of cell averaged densities

	double a,b;			//integration limits
	double SigmaCells;	//surface density in cells
	double rho_0;		//central density in Msun / kpc^3

	//double R = 13.9211647546;  //radius for comparison
	//double R = 6.85009694274;  //radius for comparison
	double R = 0.220970869121; //radius for comparison

	FILE *fp;	//file pointer
	char fname[200];	//filename

	//allow for command line arguments
	if(argc>=2)
		nz = atoi(argv[1]); //z-dir grid size is the second command line argument
	if(argc>=3)
		R = atof(argv[2]);	//radius is the third command line argument

	//compute cell width
	dz = (z_max - z_min)/((double) nz);
	dzout = (z_max - z_min)/((double) nzout);


	printf("1/2 Number of cells in (positive) z-direction = %d\n",nz);
	printf("Radius in kpc were z-column is located = %e\n",R);

	//set disk properties
	SetDiskProperties(&d);

	//print disk properties to screen
	//PrintDiskProperties(d);

	//print central surface density
	printf("Central Surface Density = %e [Msun / kpc^2]\n",SigmaDisk(0,d));

	//print surface density at R
	printf("Surface Density at R = %e is Sigma = %e\n",R,SigmaDisk(R,d));

	//allocate arrays
	z     = (double *) calloc(nz,sizeof(double));
	rho   = (double *) calloc(nz,sizeof(double));
	rhog  = (double *) calloc(nz,sizeof(double));
	p     = (double *) calloc(nz,sizeof(double));
	gradp = (double *) calloc(nz,sizeof(double));
	gza   = (double *) calloc(nz,sizeof(double));
	drhoa   = (double *) calloc(nz,sizeof(double));
	zout     = (double *) calloc(nzout,sizeof(double));
	rhoout     = (double *) calloc(nzout,sizeof(double));
	pout     = (double *) calloc(nzout,sizeof(double));
	gradpout     = (double *) calloc(nzout,sizeof(double));
	gzaout     = (double *) calloc(nzout,sizeof(double));

	//set z and (un-normalized) rho arrays
	//compute profile integration
	SigmaCells = 0.0;
	for(i=0;i<nz;i++)
	{
		//create z column of cell centers
		z[i] = 0.5*dz + ((double) i)*dz;
		//printf("i %d zl %e z %e zu %e\n",i,z[i]-0.5*dz,z[i],z[i]+0.5*dz);

		//compute surface density in each cell
		a = z[i] - 0.5*dz;
		b = z[i] + 0.5*dz;
		rho[i] = cell_surface_density(R, a, b, d)/dz;

		//compute acceleration in each cell
		gza[i] = gz(R, z[i], d);

		SigmaCells += rho[i]*dz;
	}
	printf("SigmaCells = %e\n",SigmaCells);

	//Compute central density
	rho_0 = 0.5 * SigmaDisk(R,d) / SigmaCells;

	//renormalize density
	//and compute rhog
	//and compute pressure
	SigmaCells = 0;
	for(i=0;i<nz;i++)
	{
		rho[i] *= rho_0;
		//if(rho[i]<1.0)
		//	rho[i] = 1.0;
		SigmaCells += rho[i]*dz;

		//rho * gz
		rhog[i] = rho[i] * gza[i];

		//p
		p[i] = pow(d.cs,2)*rho[i];

		//compute pressure gradient as cs^2 grad rho
		//a = rho_0*rho_disk_vert(R, z[i]-0.25*dz, d.B_h, d.R_h, d.B_d, d.R_d, d.Z_d);
		//b = rho_0*rho_disk_vert(R, z[i]+0.25*dz, d.B_h, d.R_h, d.B_d, d.R_d, d.Z_d);
		//gradp[i] = pow(d.cs,2)*(b-a)/(dz);

		if(i==0)
		{
			gradp[i] = pow(d.cs,2)*(rho[i+1]-rho[i])/(z[i+1]-z[i]);
		}else if(i==nz-1){
			gradp[i] = pow(d.cs,2)*(rho[i]-rho[i-1])/(z[i]-z[i-1]);
		}else{
			gradp[i] = pow(d.cs,2)*(rho[i+1]-rho[i-1])/(z[i+1]-z[i-1]);
		}

	}
	printf("Corrected SigmaCells %e 0.5*Sigma(%e) %e\n",SigmaCells,R,0.5*SigmaDisk(R,d));

	int iter = 0;
	int niter_max = 10000000;
	//output density profile with z-height
	sprintf(fname,"density_data/density.%06d.txt",iter);
	fp = fopen(fname,"w");
	for(i=0;i<nz;i++)
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",z[i],rho[i],rho[i]*gza[i],p[i],gradp[i]);
	fclose(fp);
	iter++;

	//OK, we know gradp and rhog are not equal, so try to correct
	//try moving density to be half-way to cs^2 drho/dz / g
	int flag = 1;
	double Lnorm;
	double g;
	double drho;
	double fac = 1.0e-6;
	while(flag)
	{
		//reset error
		Lnorm = 0;

		//begin with computing the difference between current rho and 
		//needed rho
		//we want -rho g = cs^2 drho/dz
		//so compute a = cs^2 drho/dz / -g
		//and record b = rho
		//then take 0.5*(b-a) and add it to rho
		SigmaCells = 0;
		for(i=0;i<nz-1;i++)
		{
			//try again adjusting density gradient
			//we want -rho g = cs^2 drho/dz
			//find drho = -rho g * dz / cs^2	
			a = rho[i]*gza[i] * dz / pow(d.cs,2);

			//find drho
			/*if(i==0)
			{
				b = rho[i+1]-rho[i];
			}else if(i==nz-1){
				b = rho[i] - rho[i-1];
			}else{
				b = 0.5*(rho[i+1]-rho[i-1]);
			}*/
			b = rho[i+1]-rho[i];

			//correct rho
			drho = (b-a);
			if(drho<0)
				drho *= 0.5;
			if(drho>0)
				drho *= 0.4;
			if(fabs(drho)>fac*rho[i+1])
			{
				drho = (drho/fabs(drho))*fac*rho[i+1];
			}
			drhoa[i] = drho;
			//printf("z %e rhogdz/cs^2 %e Drho %e rho[i] %e rho[i+1] %e drho %e\n",z[i],a,b,rho[i],rho[i+1],drho);


			/*
			//cs^2 drho/dz / (-g);
			//a = log(gradp[i]/(-1*gz(R, z[i], d)));
			//b = log(rho[i]);
			a = gradp[i]/(gza[i]);
			b = rho[i];		

			//correct rho
			//rho[i] *= exp(0.5*(b-a));
			drho = (a-b);
			if(drho<0)
				drho *= 0.5;
			if(drho>0)
				drho *= 0.4;
			if(fabs(drho)>0.5*rho[i])
			{
				drho = (drho/fabs(drho))*0.5*rho[i];
			}


			printf("iter %d z %e gradp %e gz %e gradp/gz %e rho %e (rho-gradp/gz)/gradp/gz % e change % e\n",iter,z[i],gradp[i],gza[i],a,b,(b-a)/a,drho);
			*/

			//correct rho
			//rho[i+1] -= drho;

			//add to surface density
			//SigmaCells += rho[i]*dz;

			//printf("z %e new rho %e\n",z[i],rho[i]);


		}
		for(i=1;i<nz;i++)
		{
			rho[i] -= drhoa[i-1];
			SigmaCells += rho[i]*dz;

		}

		//fix surface density and
		for(i=0;i<nz;i++)
		{
			rho[i] *= 0.5 * SigmaDisk(R,d) / SigmaCells;
			//printf("z %e rho %e\n",z[i],rho[i]);
		}
		

		//correct gradp
		SigmaCells = 0;
		for(i=0;i<nz-1;i++)
		{
			/*if(i==0)
			{
				gradp[i] = pow(d.cs,2)*(rho[i+1]-rho[i])/(z[i+1]-z[i]);
			}else if(i==nz-1){
				gradp[i] = pow(d.cs,2)*(rho[i]-rho[i-1])/(z[i]-z[i-1]);
			}else{
				gradp[i] = pow(d.cs,2)*(rho[i+1]-rho[i-1])/(z[i+1]-z[i-1]);
			}*/
			gradp[i] = pow(d.cs,2)*(rho[i+1]-rho[i])/(z[i+1]-z[i]);

			SigmaCells += rho[i]*dz;

			if(z[i]<5.0)
				Lnorm += fabs(rho[i]*gza[i]-gradp[i]);
		}
		//printf("iter %d SigmaDisk %e SigmaCells %e Lnorm %e\n",iter,0.5*SigmaDisk(R,d),SigmaCells,Lnorm);
		iter++;

		//output density profile with z-height
		//if( ((iter<100000)&&(!(iter%100))) || !(iter%(niter_max/100)))
		if( ((!(iter%4))&(iter<200)) || ((!(iter%40))&(iter<4000)) )
		{
			sprintf(fname,"density_data/density.%06d.txt",iter);
			fp = fopen(fname,"w");
			for(i=0;i<nz;i++)
				fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",z[i],rho[i],rho[i]*gza[i],p[i],gradp[i]);
			fclose(fp);
		}

		if(iter>niter_max)
		{
			flag = 0;
			sprintf(fname,"density_data/density.%06d.txt",iter);
			fp = fopen(fname,"w");
			for(i=0;i<nz;i++)
				fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",z[i],rho[i],rho[i]*gza[i],p[i],gradp[i]);
			fclose(fp);

			int iz;
			sprintf(fname,"test.out.txt");
			fp = fopen(fname,"w");


			for(iz=0;iz<nzout;iz++)
				zout[iz] = 0.5*dzout + ( (double) iz)*dzout;

			iz=0;
			for(i=0;i<nz;i++)
			{
				if(z[i]>zout[iz])
				{
					rhoout[iz] /= dzout;
					pout[iz]    = pow(d.cs,2)*rhoout[iz];
					gzaout[iz] = gz(R, zout[iz], d);
					if(iz<nzout-1)
						gradpout[iz] = pow(d.cs,2)*(rhoout[iz+1]-rhoout[iz])/(zout[iz+1]-zout[iz]);

					iz++;
					if(iz>=nzout)
						break;
				}
				rhoout[iz] += rho[i]*dz;
			}

			for(i=0;i<nzout-1;i++)
			{
				fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",zout[i],rhoout[i],rhoout[i]*gzaout[i],pout[i],gradpout[i]);

			}
			fclose(fp);

		}

		//printf("-------------- Lnorm = %e\n",Lnorm);
	}



	//free memory
	free(z);
	free(rho);
	free(rhog);
	free(p);
	free(gradp);
}