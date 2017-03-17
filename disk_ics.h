#ifndef DISK_ICS_H
#define DISK_ICS_H

struct DiskICs{

  //some constants
  double l_s;// = 3.086e21 # length scale, centimeters in a kiloparsec
  double m_s;// = 1.99e33 # mass scale, g in a solar mass
  double t_s;// = 3.154e10 # time scale, seconds in a kyr
  double d_s;// = m_s / l_s**3 # density scale, M_sun / kpc^3
  double v_s;// = l_s / t_s # velocity scale, kpc / kyr
  double p_s;// = d_s*v_s**2 # pressure scale, M_sun / kpc kyr^2
  double G;// = 6.67259e-8 # in cm^3 g^-1 s^-2
  double mp;// = 1.67e-24 # proton mass in grams
  //G = G / l_s**3 * m_s * t_s**2 # in kpc^3 / M_sun / kyr^2
  double KB;// = 1.3806e-16 # boltzmann constant in cm^2 g / s^2 K
  double M_vir;// = 1e12
  double M_d;// = 6.5e10
  double M_h;// = M_vir - M_d
  double R_vir;// = 261 # MW viral radius in kpc
  double c_vir;// = 20 
  double R_h;// = R_vir / c_vir # halo scale radius in kpc
  double R_d;// = 3.5 # stellar disk scale length in kpc
  double Z_d;// = 3.5/5.0 # disk scale height in kpc
  double R_g;// = 2*R_d # gas disk scale length in kpc
  double T;// = 1e4 # gas temperature, 10^4 K
  double cs;// = np.sqrt(KB*T/(0.6*mp))*t_s/l_s # isothermal sound speed
  double v_to_kmps;// = l_s/t_s/100000
  double kmps_to_kpcpkyr;// = 1.0220122e-6
  double gamma;// = 1.001
  double phi_0_h; //center of halo potential
  double phi_0_d; //center of disk potential
  double B_h;
  double B_d;
  
};

//set disk constant properties
void SetDiskProperties(struct DiskICs *d);

//print disk properties to screen
void PrintDiskProperties(struct DiskICs d);

//# function for halo concentration, NFW profile
double nfw_func(double xx);

//vertical acceleration in //kpc / kyr^2
double gz(double R, double z, DiskICs d);

//surface density in solar masses / kpc^2
double SigmaDisk(double R, DiskICs d);

//# NFW halo potential with Phi_0_h divided out
double phi_halo(double R, double z, double R_h);

//# Miyamoto-Nagai disk with Phi_0_d divided out
double phi_disk(double R, double z, double R_d, double Z_d);

//un-normalized
//vertical density profile of the initial disk model
//before iterating for hydrostatic balance.
double rho_disk_vert(double R, double z, double B_h, double R_h, double B_d, double R_d, double Z_d);

//un-normalized vertical disk density profile
//suitable for integration
double rho_disk_vert_integrand(double z, void *fp);

//compute un-normalized cell surface density
double cell_surface_density(double R, double a, double b, DiskICs d);

#endif //DISK_ICS_H