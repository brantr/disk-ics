#include <stdio.h>
#include <math.h>
#include "disk_ics.h"
#include "rk_int.h"

//surface density in solar masses / kpc^2
double SigmaDisk(double R, DiskICs d)
{
  return 0.25*d.M_d*exp(-R/d.R_g) / (2.*M_PI*d.R_g*d.R_g);
}

//# function for halo concentration, NFW profile
double nfw_func(double xx)
{
  return log(1+xx) - xx / (1.+xx);
}

//vertical acceleration in //kpc / kyr^2
double gz(double R, double z, DiskICs d)
{
  //x_pos = np.linspace(xmin+0.5*dx,xmax-0.5*dx,nx)

  //r_disk = 13.9211647546
  //#r_disk = 6.85009694274
  //#r_disk = 0.220970869121
  double r_halo = sqrt(z*z + R*R);
  double x = r_halo/d.R_h;
  double a_halo_z = 0;

  //# calculate acceleration due to NFW halo & Miyamoto-Nagai disk
  double a_halo = -1.*d.phi_0_h * nfw_func(x) / (x*r_halo); //kpc / kyr^2
  //double a_halo = -1.*d.phi_0_h*d.R_h * nfw_func(x) / (r_halo*r_halo); // phi_0_h not defined consistently across python scripts
  //double a_disk_z = -1.*d.G * d.M_d * z * (d.R_d + sqrt(z*z + d.Z_d*d.Z_d)) / ( pow(R*R + pow(d.R_d + sqrt(z*z + d.Z_d*d.Z_d),2),1.5) * sqrt(z*z + d.Z_d*d.Z_d) );

  //phi(R,z) miyamoto-nagai
  //-GM / sqrt( R^2 + (a + sqrt(b^2 + z^2)^2))



  //following checked against numerical derivative of phi_disk.
  double a = d.R_d;
  double b = d.Z_d;
  double a_disk_z = z * (1. + sqrt( pow(b/a,2) + pow(z/a,2) ))/( a*a * sqrt(pow(b/a,2) + pow(z/a,2)) * pow( pow(R/a,2) + pow(1 + sqrt(pow(b/a,2) + pow(z/a,2)),2) ,1.5));
  a_disk_z *= -1.*d.phi_0_d;


  //return a_disk_z;


  //# total acceleration is the sum of the halo + disk components


  double z1,z2;
  double phi1,phi2;
  if(z==0)
  {
    z1 = 0;
    z2 = 1.0e-3;
  }else{
    z1 = z-0.1*z;
    z2 = z+0.1*z;
  }
  //phi1 = d.phi_0_h*phi_halo(R,z1,d.R_h) + d.phi_0_d*phi_disk(R,z1,d.R_d,d.Z_d);
  //phi2 = d.phi_0_h*phi_halo(R,z2,d.R_h) + d.phi_0_d*phi_disk(R,z2,d.R_d,d.Z_d);
  phi1 = d.phi_0_d*phi_disk(R,z1,d.R_d,d.Z_d);
  phi2 = d.phi_0_d*phi_disk(R,z2,d.R_d,d.Z_d);
  //return a_disk_z;
  //return -dPhi/dz
  return -1.0*(phi2-phi1)/(z2-z1);
}

//# NFW halo potential with Phi_0_h divided out
double phi_halo(double R, double z, double R_h)
{
  return -1.*log(1. + sqrt(R*R + z*z)/R_h)/(sqrt(R*R + z*z)/R_h);
}

//# Miyamoto-Nagai disk with Phi_0_d divided out
double phi_disk(double R, double z, double R_d, double Z_d)
{
  return -1.0 / sqrt(pow(R/R_d,2) + pow(1. + sqrt(pow(z/R_d,2) + pow(Z_d/R_d,2)),2));
}

//un-normalized
//vertical density profile of the initial disk model
//before iterating for hydrostatic balance.
double rho_disk_vert(double R, double z, double B_h, double R_h, double B_d, double R_d, double Z_d)
{
  double A_h = phi_halo(R, 0, R_h);
  double A_d = phi_disk(R, 0, R_d, Z_d);
  //return exp(B_h*(A_h - phi_halo(z,R,R_h)) + B_d*(A_d - phi_disk(z,R,R_d,Z_d)));
  return exp(B_d*(A_d - phi_disk(z,R,R_d,Z_d)));

}

//un-normalized vertical disk density profile
//suitable for integration
double rho_disk_vert_integrand(double lnz, void *fp)
{
  double  z = exp(lnz);
  double *f = (double*) fp;
  double R   = f[0];
  double B_h = f[1];
  double R_h = f[2];
  double B_d = f[3];
  double R_d = f[4];
  double Z_d = f[5];

  return z*rho_disk_vert(R,z,B_h,R_h,B_d,R_d,Z_d);
}

//compute un-normalized cell surface density
double cell_surface_density(double R, double a, double b, DiskICs d)
{
  double lna, lnb;
  double fp[6];   //parameters
  double dxinit;
  double eps = 1.0e-6; //tolerance
  //prevent a=0
  if(a<1.0e-9)
  {
    lna = log(1.0e-9);
  }else{
    lna = log(a);
  }
  //assume b>0
  lnb = log(b);

  //initial step size
  dxinit = 0.01*(lnb-lna);

  //parameters
  fp[0] = R;
  fp[1] = d.B_h;
  fp[2] = d.R_h;
  fp[3] = d.B_d;
  fp[4] = d.R_d;
  fp[5] = d.Z_d;

  return integrate(rho_disk_vert_integrand,&fp[0],6,lna,lnb,dxinit,eps);
}

void SetDiskProperties(struct DiskICs *x)
{
  //some constants
  x->l_s = 3.086e21;// # length scale, centimeters in a kiloparsec
  x->m_s = 1.99e33;// # mass scale, g in a solar mass
  x->t_s = 3.154e10;// # time scale, seconds in a kyr
  x->d_s = x->m_s / pow(x->l_s,3);// # density scale, M_sun / kpc^3
  x->v_s = x->l_s / x->t_s;// # velocity scale, kpc / kyr
  x->p_s = x->d_s*pow(x->v_s,2); //# pressure scale, M_sun / kpc kyr^2
  x->G = 6.67259e-8;// # in cm^3 g^-1 s^-2
  x->mp = 1.67e-24;// # proton mass in grams
  x->G = x->G / pow(x->l_s,3) * x->m_s * pow(x->t_s,2) ;//# in kpc^3 / M_sun / kyr^2
  x->KB = 1.3806e-16;// # boltzmann constant in cm^2 g / s^2 K
  x->M_vir = 1e12;//
  x->M_d = 6.5e10;
  x->M_h = x->M_vir - x->M_d;
  x->R_vir = 261.;// # MW viral radius in kpc
  x->c_vir = 20;// 
  x->R_h = x->R_vir / x->c_vir;// # halo scale radius in kpc
  x->R_d = 3.5;// # stellar disk scale length in kpc
  x->Z_d = 3.5/5.0;// # disk scale height in kpc
  x->R_g = 2*x->R_d;// # gas disk scale length in kpc
  x->T = 1e4;// # gas temperature, 10^4 K
  x->cs = sqrt(x->KB*x->T/(0.6*x->mp))*x->t_s/x->l_s;// # isothermal sound speed
  x->v_to_kmps = x->l_s/x->t_s/100000.;//
  x->kmps_to_kpcpkyr = 1.0220122e-6;//
  x->gamma = 1.001;//
  //# define phi_0_h
  x->phi_0_h = x->G * x->M_h / (x->R_h * nfw_func(x->c_vir)); // (kpc/kyr)^2
  x->B_h = x->phi_0_h / pow(x->cs,2);  // unit free
  x->phi_0_d = x->G * x->M_d / x->R_d; // (kpc/kyr)^2
  x->B_d = x->phi_0_d / pow(x->cs,2);  // unit free
  return;
}

void PrintDiskProperties(struct DiskICs x)
{
  //some constants
  printf("x->l_s = %e;// # length scale, centimeters in a kiloparsec\n",x.l_s);
  printf("x->m_s = %e;// # mass scale, g in a solar mass\n",x.m_s);
  printf("x->t_s = %e;// # time scale, seconds in a kyr\n",x.t_s);
  printf("x->d_s = x->m_s / pow(x->l_s,3) = %e;// # density scale, M_sun / kpc^3\n",x.d_s);
  printf("x->v_s = x->l_s / x->t_s = %e;// # velocity scale, kpc / kyr\n",x.v_s);
  printf("x->p_s = x->d_s*pow(x->v_s,2) = %e; //# pressure scale, M_sun / kpc kyr^2\n",x.p_s);
  printf("x->G = %e;// # in cm^3 g^-1 s^-2\n",x.G * pow(x.l_s,3) / x.m_s / pow(x.t_s,2));
  printf("x->mp = %e;// # proton mass in grams\n",x.mp);
  printf("x->G = x->G / pow(x->l_s,3) * x->m_s * pow(x->t_s,2) = %e ;//# in kpc^3 / M_sun / kyr^2\n",x.G);
  printf("x->KB = %e;// # boltzmann constant in cm^2 g / s^2 K\n",x.KB);
  printf("x->M_vir = %e;//\n",x.M_vir);
  printf("x->M_d = %e;\n",x.M_d);
  printf("x->M_h = x->M_vir - x->M_d = %e;\n",x.M_h);
  printf("x->R_vir = %e;// # MW viral radius in kpc\n",x.R_vir);
  printf("x->c_vir = %e;// \n",x.c_vir);
  printf("x->R_h = x->R_vir / x->c_vir = %e;// # halo scale radius in kpc\n",x.R_h);
  printf("x->R_d = %e;// # stellar disk scale length in kpc\n",x.R_d);
  printf("x->Z_d = 3.5/5.0 = %e;// # disk scale height in kpc\n",x.Z_d);
  printf("x->R_g = 2*x->R_d = %e;// # gas disk scale length in kpc\n",x.R_g);
  printf("x->T = 1e4 = %e;// # gas temperature, 10^4 K\n",x.T);
  printf("x->cs = sqrt(x->KB*x->T/(0.6*x->mp))*x->t_s/x->l_s = %e;// # isothermal sound speed\n",x.cs);
  printf("x->v_to_kmps = x->l_s/x->t_s/100000. = %e;//\n",x.v_to_kmps);
  printf("x->kmps_to_kpcpkyr = 1.0220122e-6 = %e;//\n",x.kmps_to_kpcpkyr);
  printf("x->gamma = 1.001 = %e;//\n",x.gamma);
  return;
}