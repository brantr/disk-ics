/*! \file rk_int.h
 *  \brief Function declarations for numerical integration routines.
 */
#ifndef RK_INT_H
#define RK_INT_H

/*! \fn double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps)
 *  \brief Numerical integration routine using 5th-order Runge-Kutta.
 */
double integrate(double (*FUNC)(double,void*), void *fp ,int np,double a,double b,double dxinit, double eps);
/*! \fn void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp)
 *  \brief Runge-Kutta step in numerical integration routine
 */
void RUNGE5VAR(double *y,double dydx,double *x,double htry,double eps,double yscale,double *hnext,double (*DERIVS)(double,void*), void *fp);
/*! \fn void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp)
 *  \brief Function to advance the RK solution in the numerical integration.
 */
void RUNGE(double y,double dydx,double x,double h,double *yout,double *yerr,double (*DERIVS)(double,void*),void *fp);



/*! \fn double midpoint_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level)
 *  \brief Midpoint rule integration
 */
double midpoint_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);

/*! \fn double trapezoid_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);
 *  \brief Trapezoid rule integration
 */
double trapezoid_rule_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int level);

/*! \fn double romberg_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int m, int k);
 *  \brief Romberg integration
 */
double romberg_integration(double(*func)(double,void*), void *fp, int np, double a, double b, int m, int k);

#endif
