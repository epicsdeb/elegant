/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: track_through_rf_deflector()
 * purpose: track particles through an RF deflector
 * 
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

long find_nearby_array_entry(double *entry, long n, double key);
void set_up_rftm110(RFTM110 *rf_param, double **initial, long n_particles, double pc_central);
void set_up_rfdf(RFDF *rf_param, double **initial, long n_particles, double pc_central);
double linear_interpolation(double *y, double *t, long n, double t0, long i);
void set_up_mrfdf(MRFDF *rf_param, double **initial, long n_particles, double pc_central);

void track_through_rf_deflector(
                                double **final, 
                                RFDF *rf_param,
                                double **initial,
                                long n_particles,
                                double pc_central,
				double L_central,
				double zEnd,
				long pass
                                )
{
  double t_first;     /* time when first particle crosses cavity center */
  double t_part;      /* time at which a particle enters cavity */
  double Estrength;    /* |e.V.L/nSections|/(gap.m.c^2) */
  double x, xp, y, yp;
  double beta, px, py, pz, beta_z, pc, gamma, dpx;
  double omega, k, Ephase, voltFactor, t0;
  double cos_tilt, sin_tilt, dtLight, tLight;
  double length;
  long ip, is, n_kicks;

  omega = 2*PI*rf_param->frequency;
  k = omega/c_mks;

  if ((n_kicks = rf_param->n_kicks)<=0) {
    if ((n_kicks = (rf_param->length/(PIx2/k)*100))<1) {
      if (rf_param->length==0)
        n_kicks = 1;
      else
        n_kicks = 11;
    }
  }
  if (n_kicks%2==0)
    n_kicks += 1;
  
  if (rf_param->frequency==0) 
    bombElegant("RFDF cannot have frequency=0", NULL);
  if (rf_param->voltage==0 ||
      (rf_param->startPass>=0 && pass<rf_param->startPass) ||
      (rf_param->endPass>=0 && pass>rf_param->endPass)) {
    if (isSlave || !notSinglePart) 
      exactDrift(initial, n_particles, rf_param->length);
    return;
  }
    
  if (!rf_param->initialized || !rf_param->fiducial_seen)
      set_up_rfdf(rf_param, initial, n_particles, pc_central);

  gamma = sqrt(sqr(pc_central)+1);
  beta  = pc_central/gamma;
  if (pass==0) {
    rf_param->Ts = rf_param->t_first_particle;
    if (rf_param->alignWaveforms)
      rf_param->Ts = 0;
  }
  else
    rf_param->Ts += L_central/(beta*c_mks);
  t0 = rf_param->Ts;

  if (rf_param->n_Vpts) {
    long i_volt;
    if (rf_param->voltageIsPeriodic && t0>rf_param->V_tFinal) 
      t0 = fmod(t0-rf_param->V_tInitial, rf_param->V_tFinal-rf_param->V_tInitial)+rf_param->V_tInitial;

    /* find position within voltage waveform array */
    i_volt = find_nearby_array_entry(rf_param->t_Vf, rf_param->n_Vpts, t0);
    voltFactor = linear_interpolation(rf_param->Vfactor, rf_param->t_Vf, rf_param->n_Vpts, t0, i_volt);
  } else {
    voltFactor = 1;
  }

  voltFactor *= (1+rf_param->fse);

  voltFactor *= (gauss_rn_lim(1.0, rf_param->voltageNoise, 2, random_3) +
     (rf_param->voltageNoiseGroup
      ? rf_param->groupVoltageNoise*GetNoiseGroupValue(rf_param->voltageNoiseGroup)
      : 0));

  if (voltFactor==0) {
    exactDrift(initial, n_particles, rf_param->length);
    return;
  }
  
  t_first = rf_param->t_first_particle;
  length = rf_param->length/n_kicks;
  Ephase = (rf_param->phase
	    + gauss_rn_lim(0.0, rf_param->phaseNoise, 2, random_3)
	    + (rf_param->phaseNoiseGroup
	       ? rf_param->groupPhaseNoise*GetNoiseGroupValue(rf_param->phaseNoiseGroup)
	       : 0)
	    )*PI/180.0 
    + omega*(rf_param->time_offset - t_first);
#ifdef DEBUG
  fprintf(stderr, "t_first = %e s, Ephase = %e deg\n", t_first, Ephase*180/PI);
#endif
  if (!rf_param->standingWave)  {
    dtLight = length/c_mks/2;
  } else {
    dtLight = 0;
    Ephase -= omega*rf_param->length/2/c_mks;
  }
#ifdef DEBUG
  fprintf(stderr, "dtLight=%e, length=%e, omega*dtLight=%e, Ephase = %e deg\n", dtLight, length, omega*dtLight*180/PI, Ephase*180/PI);
#endif

  cos_tilt = cos(rf_param->tilt);
  sin_tilt = sin(rf_param->tilt);
  Estrength = voltFactor*(particleCharge*rf_param->voltage/n_kicks)/(particleMass*sqr(c_mks));

  if (isSlave || !notSinglePart) {
    if (rf_param->dx || rf_param->dy || rf_param->dz)
      offsetBeamCoordinates(initial, n_particles, rf_param->dx, rf_param->dy, rf_param->dz);
    if (rf_param->tilt)
      rotateBeamCoordinates(initial, n_particles, rf_param->tilt);
    for (ip=0; ip<n_particles; ip++) {
      x  = initial[ip][0];
      xp = initial[ip][1];
      y  = initial[ip][2];
      yp = initial[ip][3];
      pc = pc_central*(1+initial[ip][5]);
      pz = pc/sqrt(1+sqr(xp)+sqr(yp));
      px = xp*pz;
      py = yp*pz;
      beta = pc/sqrt(1+sqr(pc));
      t_part = initial[ip][4]/(c_mks*beta);
      tLight = 0;
#ifdef DEBUG
      fprintf(stderr, "start coord[%ld] = %e, %e, %e, %e, %e, %e\n",
	      ip, x, xp, y, yp, initial[ip][4], initial[ip][5]);
#endif
      for (is=0; is<=n_kicks; is++) {
	beta_z = pz/pc;
	if (is==0 || is==n_kicks) {
	  /* first half-drift and last half-drift */
	  t_part += (length*sqrt(1+sqr(xp)+sqr(yp))/(2*c_mks*beta_z));
          tLight = dtLight;
	  x += xp*length/2;
	  y += yp*length/2;
	  if (is==n_kicks)
	    break;
	} else {
	  t_part += (length*sqrt(1+sqr(xp)+sqr(yp))/(c_mks*beta_z));
          tLight += 2*dtLight; 
	  x += xp*length;
	  y += yp*length;
	}
#ifdef DEBUG
	fprintf(stdout, "ip=%ld  is=%ld  dphase=%f, phase=%f\n",
		ip, is, omega*(t_part-tLight)*180/PI, fmod((t_part-tLight)*omega+Ephase, PIx2)*180/PI);
#endif
	px += Estrength*cos((t_part-tLight)*omega + Ephase);
        if (rf_param->magneticDeflection)
          pz = sqrt(sqr(pc)-sqr(px)-sqr(py));
        pz += Estrength*k*x*sin((t_part-tLight)*omega + Ephase);
	xp = px/pz;
	yp = py/pz;
	pc = sqrt(sqr(px)+sqr(py)+sqr(pz));
      }
      beta = pc/sqrt(1+sqr(pc));
      final[ip][0] = x;
      final[ip][1] = xp;
      final[ip][2] = y;
      final[ip][3] = yp;
      final[ip][4] = t_part*c_mks*beta;
      final[ip][5] = (pc-pc_central)/pc_central;
      final[ip][6] = initial[ip][6];
#ifdef DEBUG
      fprintf(stderr, "stop  coord[%ld] = %e, %e, %e, %e, %e, %e\n",
	      ip, final[ip][0], final[ip][1], final[ip][2], final[ip][3],
	      final[ip][4], final[ip][5]);
#endif
    }
    if (rf_param->tilt)
      rotateBeamCoordinates(final, n_particles, -rf_param->tilt);
    if (rf_param->dx || rf_param->dy || rf_param->dz)
      offsetBeamCoordinates(final, n_particles, -rf_param->dx, -rf_param->dy, -rf_param->dz);
  }
  
}

/* routine: track_through_rftm110_deflector()
 * purpose: track particles through an RF deflector with more exact implementation
 *   for a tm110 mode.
 * 
 * Michael Borland, 2004
 */

void track_through_rftm110_deflector(
                                double **final, 
                                RFTM110 *rf_param,
                                double **initial,
                                long n_particles,
                                double pc_central,
				double L_central,
				double zEnd,
				long pass
                                )
{
  double t_first;     /* time when first particle crosses cavity center */
  double t_part;      /* time at which a particle enters cavity */
  double x, xp, y, yp, rho, k;
  double beta, px, py, pz, beta_x, beta_y, beta_z, pc;
  double omega, phase, phase0, Ez, cBx, cBy;
  double cos_phi, sin_phi, voltTimes2;
  double krho2, krho4, krho6, cos_2phi;
  long ip, i_volt;
  double gamma, t0;

  if (rf_param->frequency==0) 
    bombElegant("RFTM110 cannot have frequency=0", NULL);

  if (rf_param->voltage==0 ||
      (rf_param->startPass>=0 && pass<rf_param->startPass) ||
      (rf_param->endPass>=0 && pass>rf_param->endPass))
    return;

  if (!rf_param->initialized || !rf_param->fiducial_seen) 
    set_up_rftm110(rf_param, initial, n_particles, pc_central);

  gamma = sqrt(sqr(pc_central)+1);
  beta  = pc_central/gamma;
  if (pass==0) {
    rf_param->Ts = rf_param->t_first_particle;
    if (rf_param->alignWaveforms)
      rf_param->Ts = 0;
  }
  else
    rf_param->Ts += L_central/(beta*c_mks);
  t0 = rf_param->Ts;

  if (rf_param->n_Vpts) {
    if (rf_param->voltageIsPeriodic && t0>rf_param->V_tFinal) 
      t0 = fmod(t0-rf_param->V_tInitial, rf_param->V_tFinal-rf_param->V_tInitial)+rf_param->V_tInitial;

    /* find position within voltage waveform array */
    i_volt = find_nearby_array_entry(rf_param->t_Vf, rf_param->n_Vpts, t0);
    voltTimes2 = linear_interpolation(rf_param->Vfactor, rf_param->t_Vf, rf_param->n_Vpts, t0, i_volt);
    if (voltTimes2==0)
      return;
  } else {
    voltTimes2 = 1;
  }

  /* using 2*volt in expressions gives us theta=V/E */
  voltTimes2 *= 2*rf_param->voltage/(1e6*particleMassMV*particleRelSign)*
    (gauss_rn_lim(1.0, rf_param->voltageNoise, 2, random_3) +
     (rf_param->voltageNoiseGroup
      ? rf_param->groupVoltageNoise*GetNoiseGroupValue(rf_param->voltageNoiseGroup)
      : 0));

  omega = 2*PI*rf_param->frequency;
  k = omega/c_mks;
  t_first = rf_param->t_first_particle;
  phase0 = (rf_param->phase
	    + gauss_rn_lim(0.0, rf_param->phaseNoise, 2, random_3)
	    + (rf_param->phaseNoiseGroup
	       ? rf_param->groupPhaseNoise*GetNoiseGroupValue(rf_param->phaseNoiseGroup)
	       : 0)
	    )*PI/180.0 
    - omega*t_first;
  
  if(isSlave || !notSinglePart) {
    if (rf_param->tilt)
      rotateBeamCoordinates(initial, n_particles, rf_param->tilt);
    for (ip=0; ip<n_particles; ip++) {
      x  = initial[ip][0];
      y  = initial[ip][2];
      if ((rho = sqrt(x*x+y*y))>0) {
	cos_phi = x/rho;
	sin_phi = y/rho;
      } else {
	cos_phi = 1;
	sin_phi = 0;
      }
      xp = initial[ip][1];
      yp = initial[ip][3];

      pc = pc_central*(1+initial[ip][5]);
      pz = pc/sqrt(1+sqr(xp)+sqr(yp));
      px = xp*pz;
      py = yp*pz;
      beta = pc/sqrt(1+sqr(pc));
      t_part = initial[ip][4]/(c_mks*beta);
      /*     phase = omega*(t_part - t_first) + phase0;  */
      phase = omega*t_part + phase0;

      krho2 = sqr(k*rho);
      krho4 = sqr(krho2);
      krho6 = krho2*krho4;
      cos_2phi = 2*sqr(cos_phi)-1;
      Ez = -voltTimes2*k*rho*(192 - 24*krho2 + krho4)*cos_phi*sin(phase)/384.;
      cBx = voltTimes2*krho2*(384 - 32*krho2 + krho4)*
	cos_phi*cos(phase)*sin_phi/3072.;
      cBy = -voltTimes2*
	(-9216 + 2304*krho2 - 144*krho4 + 4*krho6 + 
	 1152*krho2*cos_2phi - 96*krho4*cos_2phi +
	 3*krho6*cos_2phi)*cos(phase)/18432.;
      beta_x = px/pc;
      beta_y = py/pc;
      beta_z = pz/pc;
      px += beta_z*cBy;
      py += - beta_z*cBx;
      pz += -(Ez + beta_x*cBy - beta_y*cBx) ;

      xp = px/pz;
      yp = py/pz;
      pc = sqrt(sqr(px)+sqr(py)+sqr(pz));
      beta = pc/sqrt(1+sqr(pc));
      final[ip][0] = x;
      final[ip][1] = xp;
      final[ip][2] = y;
      final[ip][3] = yp;
      final[ip][4] = t_part*c_mks*beta;
      final[ip][5] = (pc-pc_central)/pc_central;
      final[ip][6] = initial[ip][6];
    }
    if (rf_param->tilt)
      rotateBeamCoordinates(initial, n_particles, -rf_param->tilt);
  }
}


void set_up_rftm110(RFTM110 *rf_param, double **initial, long n_particles, double pc_central)
{
  long ip, i;
  double pc, beta;
  TABLE data;
  TRACKING_CONTEXT tContext;
#ifdef USE_KAHAN
  double error = 0.0; 
#endif

  if (!rf_param->fiducial_seen) {
    if (isSlave || !notSinglePart) {
      for (ip=rf_param->t_first_particle=0; ip<n_particles; ip++) {
        pc = pc_central*(1+initial[ip][5]);
        beta = pc/sqrt(1+sqr(pc));
#ifndef USE_KAHAN
        rf_param->t_first_particle += initial[ip][4]/beta/c_mks;
#else
        rf_param->t_first_particle = KahanPlus(rf_param->t_first_particle, initial[ip][4]/beta/c_mks, &error); 
#endif
      }
    }
#if USE_MPI
    if (USE_MPI && notSinglePart) {
      long n_total;
#ifndef USE_KAHAN 
      double tmp;
#endif
      if (isMaster) {
        n_particles = 0;
        rf_param->t_first_particle = 0.0;
      }
#ifndef USE_KAHAN 
      MPI_Allreduce(&(rf_param->t_first_particle), &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      rf_param->t_first_particle = tmp;
#else
      rf_param->t_first_particle = KahanParallel(rf_param->t_first_particle, error, MPI_COMM_WORLD); 
#endif
      
      MPI_Allreduce(&n_particles, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      n_particles = n_total; 
    }
#endif
    if (n_particles)
      rf_param->t_first_particle /= n_particles;
    rf_param->fiducial_seen = 1;
  }
  
  if (rf_param->initialized)
    return;

  rf_param->initialized = 1;

  getTrackingContext(&tContext);

  if (rf_param->voltageNoiseGroup) {
    DefineNoiseGroup(rf_param->voltageNoiseGroup);
    if (rf_param->phaseNoiseGroup==rf_param->voltageNoiseGroup) {
      printf("*** Warning: VOLTAGE_NOISE_GROUP and PHASE_NOISE_GROUP are identical for %s\n",
	     tContext.elementName);
      printf("This is probably a mistake!\n");
    }
  } else if (rf_param->groupVoltageNoise) {
    printf("Error: GROUP_VOLTAGE_NOISE is nonzero but VOLTAGE_NOISE_GROUP is zero for %s\n",
	   tContext.elementName);
    exitElegant(1);
  }

  if (rf_param->phaseNoiseGroup) {
    DefineNoiseGroup(rf_param->phaseNoiseGroup);
  } else if (rf_param->groupPhaseNoise) {
    printf("Error: GROUP_PHASE_NOISE is nonzero but PHASE_NOISE_GROUP is zero for %s\n",
	   tContext.elementName);
    exitElegant(1);
  }

  rf_param->initialized = 1;

  rf_param->Ts = 0;
  rf_param->n_Vpts = 0;

  if (rf_param->voltageWaveform) {
    if (!getTableFromSearchPath(&data, rf_param->voltageWaveform, 1, 0))
        bombElegant("unable to read voltage waveform for rftm110", NULL);

    if (data.n_data<=1)
        bombElegant("rftm110 voltage waveform contains less than 2 points", NULL);

    rf_param->t_Vf    = data.c1;
    rf_param->Vfactor = data.c2;
    rf_param->n_Vpts  = data.n_data;
    for (i=0; i<rf_param->n_Vpts-1; i++)
        if (rf_param->t_Vf[i]>rf_param->t_Vf[i+1])
            bombElegant("time values are not monotonically increasing in rftm110 voltage waveform", NULL);
    rf_param->V_tInitial = rf_param->t_Vf[0];
    rf_param->V_tFinal = rf_param->t_Vf[rf_param->n_Vpts-1];
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
  }
}

void set_up_rfdf(RFDF *rf_param, double **initial, long n_particles, double pc_central)
{
  long ip, i;
  double pc, beta;
  TABLE data;
  TRACKING_CONTEXT tContext;
#ifdef USE_KAHAN
  double error = 0.0; 
#endif

  if (!rf_param->fiducial_seen) {
    if (isSlave || !notSinglePart) {
      for (ip=rf_param->t_first_particle=0; ip<n_particles; ip++) {
        pc = pc_central*(1+initial[ip][5]);
        beta = pc/sqrt(1+sqr(pc));
#ifndef USE_KAHAN
        rf_param->t_first_particle += initial[ip][4]/beta/c_mks;
#else
        rf_param->t_first_particle = KahanPlus(rf_param->t_first_particle, initial[ip][4]/beta/c_mks, &error); 
#endif
      }
    }
#if USE_MPI
    if (USE_MPI && notSinglePart) {
      long n_total;
#ifndef USE_KAHAN 
      double tmp;
#endif
      if (isMaster) {
        n_particles = 0;
        rf_param->t_first_particle = 0.0;
      }
#ifndef USE_KAHAN 
      MPI_Allreduce(&(rf_param->t_first_particle), &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      rf_param->t_first_particle = tmp;
#else
      rf_param->t_first_particle = KahanParallel(rf_param->t_first_particle, error, MPI_COMM_WORLD); 
#endif
      
      MPI_Allreduce(&n_particles, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      n_particles = n_total; 
    }
#endif
    if (n_particles)
      rf_param->t_first_particle /= n_particles;
    rf_param->fiducial_seen = 1;
  }
  

  if (rf_param->initialized)
    return;

  rf_param->initialized = 1;

  getTrackingContext(&tContext);

  if (rf_param->voltageNoiseGroup) {
    DefineNoiseGroup(rf_param->voltageNoiseGroup);
    if (rf_param->phaseNoiseGroup==rf_param->voltageNoiseGroup) {
      printf("*** Warning: VOLTAGE_NOISE_GROUP and PHASE_NOISE_GROUP are identical for %s\n",
	     tContext.elementName);
      printf("This is probably a mistake!\n");
    }
  } else if (rf_param->groupVoltageNoise) {
    printf("Error: GROUP_VOLTAGE_NOISE is nonzero but VOLTAGE_NOISE_GROUP is zero for %s\n",
	   tContext.elementName);
    exitElegant(1);
  }

  if (rf_param->phaseNoiseGroup) {
    DefineNoiseGroup(rf_param->phaseNoiseGroup);
  } else if (rf_param->groupPhaseNoise) {
    printf("Error: GROUP_PHASE_NOISE is nonzero but PHASE_NOISE_GROUP is zero for %s\n",
	   tContext.elementName);
    exitElegant(1);
  }

  rf_param->initialized = 1;

  rf_param->Ts = 0;
  rf_param->n_Vpts = 0;

  if (rf_param->voltageWaveform) {
    if (!getTableFromSearchPath(&data, rf_param->voltageWaveform, 1, 0))
        bombElegant("unable to read voltage waveform for rftm110", NULL);

    if (data.n_data<=1)
        bombElegant("rftm110 voltage waveform contains less than 2 points", NULL);

    rf_param->t_Vf    = data.c1;
    rf_param->Vfactor = data.c2;
    rf_param->n_Vpts  = data.n_data;
    for (i=0; i<rf_param->n_Vpts-1; i++)
        if (rf_param->t_Vf[i]>rf_param->t_Vf[i+1])
            bombElegant("time values are not monotonically increasing in rftm110 voltage waveform", NULL);
    rf_param->V_tInitial = rf_param->t_Vf[0];
    rf_param->V_tFinal = rf_param->t_Vf[rf_param->n_Vpts-1];
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
  }
}

void track_through_multipole_deflector(
                                double **final, 
                                MRFDF *rf_param,
                                double **initial,
                                long n_particles,
                                double pc_central
                                )
{
  double t_first;     /* time when first particle reaches cavity */
  double x, xp, y, yp;
  double beta, px, py, pz, pc;
  double omega, k, t_part;
  double ptPhaseFactor, pzPhaseFactor;
  double dpx, dpy, dpz, phase;
  long ip, mode;

  if (rf_param->frequency==0) 
    bombElegant("MRFDF cannot have frequency=0", NULL);
  if (rf_param->factor==0)
    return;
    
  if (!rf_param->initialized || !rf_param->fiducial_seen)
      set_up_mrfdf(rf_param, initial, n_particles, pc_central);

  t_first = rf_param->t_first_particle;

  if (rf_param->tilt)
    rotateBeamCoordinates(initial, n_particles, rf_param->tilt);

  for (ip=0; ip<n_particles; ip++) {
    x  = initial[ip][0];
    xp = initial[ip][1];
    y  = initial[ip][2];
    yp = initial[ip][3];
    pc = pc_central*(1+initial[ip][5]);
    pz = pc/sqrt(1+sqr(xp)+sqr(yp));
    px = xp*pz;
    py = yp*pz;
    beta = pc/sqrt(1+sqr(pc));
    t_part = initial[ip][4]/(c_mks*beta);
    dpx = dpy = dpz = 0;
    
    for (mode=0; mode<5; mode++) {
      if (rf_param->frequency[mode]==0 || (rf_param->a[mode]==0 && rf_param->b[mode]==0))
        continue;
      
      omega = 2*PI*rf_param->frequency[mode];
      k = omega/c_mks;
      phase = rf_param->phase[mode]*PI/180.0  + omega*(t_part-t_first);
      ptPhaseFactor = cos(phase)/k;
      pzPhaseFactor = sin(phase);
      if (isSlave || !notSinglePart) {
        switch (mode) {
        case 0: /* dipole */
          dpx += rf_param->b[0]*ptPhaseFactor;
          dpy -= rf_param->a[0]*ptPhaseFactor;
          dpz += (rf_param->b[0]*x - rf_param->a[0]*y)*pzPhaseFactor;
          break;
        case 1: /* quadrupole */
          dpx += 2*(rf_param->b[1]*x-rf_param->a[1]*y)*ptPhaseFactor;
          dpy -= 2*(rf_param->b[1]*y+rf_param->a[1]*x)*ptPhaseFactor;
          dpz += (rf_param->b[1]*(x*x-y*y) - 2*rf_param->a[1]*x*y)*pzPhaseFactor;
          break;
        case 2: /* sextupole */
          dpx += (rf_param->b[2]*(3*x*x-3*y*y) + rf_param->a[2]*(-6*x*y))*ptPhaseFactor;
          dpy -= (rf_param->b[2]*6*x*y + rf_param->a[2]*(3*x*x-3*y*y))*ptPhaseFactor;
          dpz += (rf_param->b[2]*(ipow(x,3)-3*x*y*y) - rf_param->a[2]*(3*x*x*y-ipow(y,3)))*pzPhaseFactor;
          break;
        case 3: /* octupole  */
          dpx += (rf_param->b[3]*(4*ipow(x,3)-12*x*y*y) + rf_param->a[3]*(-12*x*x*y+4*ipow(y,3)))*ptPhaseFactor;
          dpy -= (rf_param->b[3]*(12*x*x*y - 4*ipow(y,3)) + rf_param->a[3]*(4*ipow(x,3)-12*x*y*y))*ptPhaseFactor;
          dpz += (rf_param->b[3]*(ipow(x,4)-6*ipow(x*y,2)+ipow(y,4)) - 4*rf_param->a[3]*(x*x-y*y)*x*y)*pzPhaseFactor;
          break;
        case 4: /* decapole */
          dpx += (rf_param->b[4]*(5*ipow(y,4)+5*ipow(x,4)-30*ipow(x*y,2)) + rf_param->a[4]*(-20*ipow(x,3)*y+20*x*ipow(y,3)))*ptPhaseFactor;
          dpy -= (rf_param->b[4]*(20*ipow(x,3)*y-20*x*ipow(y,3)) + rf_param->a[4]*(5*ipow(x,4)-30*ipow(x*y,2)+5*ipow(y,4)))*ptPhaseFactor;
          dpz += (rf_param->b[4]*(ipow(x,5)-10*ipow(x,3)*y*y+5*x*ipow(y,4)) - rf_param->a[4]*(ipow(y,5)+5*ipow(x,4)*y-10*y*ipow(x*y,2)))*pzPhaseFactor;
          break;
        }
      }
    }
    /* We assume the deflection is from magnetic fields, so p^2 doesn't change */
    px += dpx*rf_param->factor*particleCharge/(particleMass*sqr(c_mks));
    py += dpy*rf_param->factor*particleCharge/(particleMass*sqr(c_mks));
    pz = sqrt(sqr(pc)-sqr(px)-sqr(py));

    /* Now add the effect of longitudinal field */
    pz += dpz*rf_param->factor*particleCharge/(particleMass*sqr(c_mks));
    pc = sqrt(sqr(px)+sqr(py)+sqr(pz));

    xp = px/pz;
    yp = py/pz;
    
    beta = pc/sqrt(1+sqr(pc));
    final[ip][0] = x;
    final[ip][1] = xp;
    final[ip][2] = y;
    final[ip][3] = yp;
    final[ip][4] = t_part*c_mks*beta;
    final[ip][5] = (pc-pc_central)/pc_central;
    final[ip][6] = initial[ip][6];
    
  }

  if (rf_param->tilt)
    rotateBeamCoordinates(initial, n_particles, -rf_param->tilt);
}

void set_up_mrfdf(MRFDF *rf_param, double **initial, long n_particles, double pc_central)
{
  long ip;
  double pc, beta;
#ifdef USE_KAHAN
  double error = 0.0; 
#endif

  if (!rf_param->fiducial_seen) {
    if (isSlave || !notSinglePart) {
      for (ip=rf_param->t_first_particle=0; ip<n_particles; ip++) {
        pc = pc_central*(1+initial[ip][5]);
        beta = pc/sqrt(1+sqr(pc));
#ifndef USE_KAHAN
        rf_param->t_first_particle += initial[ip][4]/beta/c_mks;
#else
        rf_param->t_first_particle = KahanPlus(rf_param->t_first_particle, initial[ip][4]/beta/c_mks, &error); 
#endif
      }
    }
#if USE_MPI
    if (USE_MPI && notSinglePart) {
      long n_total;
#ifndef USE_KAHAN 
      double tmp;
#endif
      if (isMaster) {
        n_particles = 0;
        rf_param->t_first_particle = 0.0;
      }
#ifndef USE_KAHAN 
      MPI_Allreduce(&(rf_param->t_first_particle), &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      rf_param->t_first_particle = tmp;
#else
      rf_param->t_first_particle = KahanParallel(rf_param->t_first_particle, error, MPI_COMM_WORLD); 
#endif
      
      MPI_Allreduce(&n_particles, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      n_particles = n_total; 
    }
#endif
    if (n_particles)
      rf_param->t_first_particle /= n_particles;
    rf_param->fiducial_seen = 1;
  }
  

}
