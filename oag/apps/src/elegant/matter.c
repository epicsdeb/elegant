/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: matter.c
 * contents: track_through_matter()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"

#define SQRT_3 (1.7320508075688772)
#define AMU (1.6605e-27)
#define SQR_PI (PI*PI)

void track_through_matter(
                          double **part, long np, MATTER *matter, double Po
                          )
{
  long ip;
  double L, Nrad, *coord, theta_rms, beta, P, gamma=0.0;
  double z1, z2, dx, dy, ds, t=0.0, dGammaFactor;
  double K1, K2=0.0, sigmaTotal, probScatter=0.0, dgamma;
  long nScatters=0;
  
  log_entry("track_through_matter");

  if (particleIsElectron==0)
    bomb("MATTER element doesn't work for particles other than electrons", NULL);
  
  if ((L=matter->length)==0)
    return;

  beta = Po/sqrt(sqr(Po)+1);
  if (matter->Xo==0)
    Nrad = theta_rms = 0;
  else {
    Nrad = matter->length/matter->Xo;
    theta_rms = 13.6/particleMassMV/Po/sqr(beta)*sqrt(Nrad)*(1+0.038*log(Nrad));
  }
  dGammaFactor = 1-exp(-Nrad);
  
  if (Nrad<1e-3) {
    if (matter->Z<1 || matter->A<1 || matter->rho<=0)
      bomb("MATTER element is too thin---provide Z, A, and rho for single-scattering calculation.", NULL);
    K1 = 4*sqr(matter->Z*particleRadius/(beta*Po));
    K2 = sqr(pow(matter->Z, 1./3.)/137.036/Po);
    sigmaTotal = K1*pow(PI, 3)/(sqr(K2)+K2*SQR_PI);
    probScatter = matter->rho/(AMU*matter->A)*matter->length*sigmaTotal;
    /* fprintf(stdout, "K1=%le, K2=%le, mean expected number of scatters is %le\n", K1, K2, probScatter); */
  }
  
  for (ip=0; ip<np; ip++) {
    coord = part[ip];
    if (Nrad) {
      if (!matter->elastic) {
        P = (1+coord[5])*Po;
        gamma = sqrt(sqr(P)+1);
        beta = P/gamma;
        t = coord[4]/beta;
      }
      if (Nrad>=1e-3) {
        /* multiple scattering */
        z1 = gauss_rn(0, random_2);
        z2 = gauss_rn(0, random_2);
        coord[0] += (dx=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[1]);
        coord[1] += z2*theta_rms;
        z1 = gauss_rn(0, random_2);
        z2 = gauss_rn(0, random_2);
        coord[2] += (dy=(z1/SQRT_3 + z2)*L*theta_rms/2 + L*coord[3]);
        coord[3] += z2*theta_rms;
        ds = sqrt(sqr(L)+sqr(dx)+sqr(dy));
      }
      else {
        long sections;
        double F, theta, phi, zs, dxp, dyp, L1, prob;
        
        sections = probScatter/matter->pLimit+1;
        L1 = L/sections;
        prob = probScatter/sections;
        ds = 0;
        while (sections-- > 0) {
          if (random_2(1)<prob) {
            nScatters ++;
            /* single-scattering computation */
            /* scatter occurs at location 0<=zs<=L */
            zs = L1*random_2(1);
            /* pick a value for CDF and get corresponding angle */
            F = random_2(1);
            theta = sqrt((1-F)*K2*SQR_PI/(K2+F*SQR_PI));
            phi = random_2(1)*PIx2;
            dxp = theta*sin(phi);
            dyp = theta*cos(phi);
            /* advance to location of scattering event */
            ds += zs*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            /* scatter */
            coord[1] += dxp;
            coord[3] += dyp;
            /* advance to end of slice */
            coord[0] += dxp*(L1-zs);
            coord[2] += dyp*(L1-zs);
            ds += (L1-zs)*sqrt(1+sqr(coord[1])+sqr(coord[3]));
          } else {
            ds += L1*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[0] += coord[1]*L1;
            coord[2] += coord[3]*L1;
          }            
        }
      }
      if (!matter->elastic) {
        dgamma = gamma*dGammaFactor;
        if (matter->energyStraggle) {
          double dgamma1;
          /* very simple-minded estimate: StDev(dE) = Mean(dE)/2 */
          while ((dgamma1 = dgamma*(1+0.5*gauss_rn(0, random_2)))<0)
            ;
          dgamma = dgamma1;
        }
        gamma -= dgamma;
        P = sqrt(sqr(gamma)-1);
        coord[5] = (P-Po)/Po;
        beta = P/gamma;
        coord[4] = t*beta+ds;
      }
      else
        coord[4] += ds;
    }
    else {
      coord[0] += L*coord[1];
      coord[2] += L*coord[3];
      coord[4] += L*sqrt(1+sqr(coord[1])+sqr(coord[3]));
    }
  }
/*
  if (Nrad<1e-3)
    fprintf(stdout, "%e scatters per particle\n",
            (1.0*nScatters)/np);
*/
  
  log_exit("track_through_matter");
}




