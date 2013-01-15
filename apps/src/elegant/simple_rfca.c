/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: simple_rfca.c
 * contents: simple_rf_cavity()
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"

static char *fiducialModeChoice[4] = {
    "light", "tmean", "first", "pmaximum",
    };

long trackRfCavityWithWakes(double **part, long np, RFCA *rfca, double **accepted, 
                            double *P_central, double zEnd, long iPass, RUN *run,
                            CHARGE *charge, WAKE *wake, TRWAKE *trwake, LSCKICK *LSCKick,
                            long wakesAtEnd);


unsigned long parseFiducialMode(char *modeString)
{
    long code;

    if (!modeString)
        return FID_MODE_TMEAN;
    switch (code=match_string(modeString, fiducialModeChoice, 4, 0)) {
      case 0:
      case 1:
      case 2:
      case 3:
        return FID_MODE_LIGHT<<code;
      default:
        return 0;
        }
    }

double findFiducialTime(double **part, long np, double s0, double sOffset,
                        double p0, unsigned long mode)
{
  double tFid=0.0;
  
  if (mode&FID_MODE_LIGHT) 
    tFid =  (s0+sOffset)/c_mks;
  else if (mode&FID_MODE_FIRST) {
#if (!USE_MPI) 
    if (np)
     tFid = (part[0][4]+sOffset)/(c_mks*beta_from_delta(p0, np?part[0][5]:0.0)); 
    else
      bombElegant("0 particle for the FID_MODE_FIRST mode in findFiducialTime", NULL);
#else
    if (myid==1) {  /* If the first particle is lost, Pelegant and elegant will not have the same fiducial time */
      if (np)
	tFid = (part[0][4]+sOffset)/(c_mks*beta_from_delta(p0, np?part[0][5]:0.0)); 
      else
	bombElegant("0 particle for the FID_MODE_FIRST mode in findFiducialTime on processor 1", NULL);
    }	
    MPI_Bcast(&tFid, 1, MPI_DOUBLE, 1, MPI_COMM_WORLD);
    /*
    fprintf(stdout, "FID_MODE_FIRST mode is not supported in the current parallel version.\n");
    fprintf(stdout, "Please use serial version.\n");
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 9);
    */
#endif
  }
  else if (mode&FID_MODE_PMAX) {
    long ibest, i;
    double best;
#if (!USE_MPI)
    best = part[0][5];
#else  /* np could be 0 for some of the slave processors */ 
    if (notSinglePart) {
      if (!np || isMaster)
	best = -DBL_MAX;
      else
	best = part[0][5];
    }
    else
      best = part[0][5];     
#endif
    ibest = 0; 
    if (isSlave || !notSinglePart)
      for (i=1; i<np; i++)
        if (best<part[i][5]) {
          best = part[i][5];
          ibest = i;
        }
#if (!USE_MPI)
    tFid = (part[ibest][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ibest][5]));
#else
    if (notSinglePart) {
      if (USE_MPI) {
	double sBest;
	struct {
	  double val;
	  int rank;
	} in, out;      
	MPI_Comm_rank(MPI_COMM_WORLD, &(in.rank));
	in.val = best;
	/* find the global best value and its location, i.e., on which processor */
	MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	best = out.val;
	/* broadcast the part[ibest][4] corresponding to the global best value */
	if (in.rank==out.rank)
	  sBest = part[ibest][4];
	MPI_Bcast(&sBest, 1, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);
	tFid = (sBest+sOffset)/(c_mks*beta_from_delta(p0, best));
      }
    }
    else
      tFid = (part[ibest][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ibest][5]));
#endif
  }
  else if (mode&FID_MODE_TMEAN) {
    double tsum=0;
    long ip;
#ifdef USE_KAHAN
    double error = 0.0; 
#endif
    
    if (isSlave || !notSinglePart) {
      for (ip=tsum=0; ip<np; ip++) {
#ifndef USE_KAHAN     
	tsum += (part[ip][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ip][5]));
#else
        tsum = KahanPlus(tsum, (part[ip][4]+sOffset)/(c_mks*beta_from_delta(p0, part[ip][5])), &error); 
#endif	
      }
    }
#if (!USE_MPI)
    tFid = tsum/np;
#else
    if (notSinglePart) {
      if (USE_MPI) {
	double tmp; 
	long np_total;
      
	if (isMaster) {
	  tsum = 0.0;
	  np = 0;
	}
	MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); 
#ifndef USE_KAHAN
	MPI_Allreduce(&tsum, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	tmp = KahanParallel (tsum, error, MPI_COMM_WORLD);
#endif
	tFid = tmp/np_total; 
      }
    }
    else
      tFid = tsum/np; 
#endif
  }
  else
    bombElegant("invalid fiducial mode in findFiducialTime", NULL);
#ifdef DEBUG
  printf("Fiducial time (mode %lx): %21.15e\n", mode, tFid);
#endif
  return tFid;
}


long simple_rf_cavity(
    double **part, long np, RFCA *rfca, double **accepted, double *P_central, double zEnd
    )
{
  return trackRfCavityWithWakes(part, np, rfca, accepted, P_central, zEnd, 0,
                         NULL, NULL, NULL, NULL, NULL, 0);
}

long trackRfCavityWithWakes
  (
   double **part, long np, 
   RFCA *rfca, 
   double **accepted, 
   double *P_central, double zEnd,
   long iPass,
   RUN *run,
   CHARGE *charge,
   WAKE *wake,
   TRWAKE *trwake,
   LSCKICK *LSCKick,
   long wakesAtEnd
   )
{
    long ip, same_dgamma, nKicks, linearize, ik, matrixMethod;
    double timeOffset, dc4, x, xp;
    double P, gamma, gamma1, dgamma=0.0, dgammaMax=0.0, phase, length, dtLight, volt, To;
    double *coord, t, t0, omega, beta_i, tau, dt, tAve=0, dgammaAve=0;
    long useSRSModel = 0;
    static long been_warned = 0;
    double dgammaOverGammaAve = 0;
    long dgammaOverGammaNp = 0;
    long lockPhase = 0;
#ifdef USE_KAHAN
    double error = 0.0; 
#endif   
#if USE_MPI
    long np_total, np_tmp;
    long i;
    double error_sum=0.0, *sumArray, *errorArray;
#endif

    matrixMethod = 0;
    if (rfca->bodyFocusModel) {
      char *modelName[2] = { "none", "srs" };
      switch (match_string(rfca->bodyFocusModel, modelName, 2, 0)) {
      case 0:
        break;
      case 1:
        useSRSModel = 1;
        matrixMethod = 1;
        break;
      default:
        fprintf(stderr, "Error: bodyFocusModel=%s not understood for RFCA\n", rfca->bodyFocusModel);
        exitElegant(1);
        break;
      }
    }
    
    if (!been_warned) {        
        if (rfca->freq<1e3 && rfca->freq)  {
            fprintf(stdout, "\7\7\7warning: your RFCA frequency is less than 1kHz--this may be an error\n");
            fflush(stdout);
            been_warned = 1;
            }
        if (fabs(rfca->volt)<100 && rfca->volt) {
            fprintf(stdout, "\7\7\7warning: your RFCA voltage is less than 100V--this may be an error\n");
            fflush(stdout);
            been_warned = 1;
            }
        if (been_warned) {
            fprintf(stdout, "units of parameters for RFCA are as follows:\n");
            fflush(stdout);
            print_dictionary_entry(stdout, T_RFCA, 0, 0);
            }
        }
    if (isSlave) {
      if (!part)
        bombElegant("NULL particle data pointer (trackRfCavityWithWakes)", NULL);
    }
    if (isSlave || !notSinglePart) {
      for (ip=0; ip<np; ip++)
	if (!part[ip]) {
	  fprintf(stderr, "NULL pointer for particle %ld (trackRfCavityWithWakes)\n", ip);
	  fflush(stderr);
#if USE_MPI
          MPI_Abort(MPI_COMM_WORLD, 1);
#endif
	  abort();
	}
    }
    if (!rfca)
        bombElegant("NULL rfca pointer (trackRfCavityWithWakes)", NULL);

#if (!USE_MPI)
    if (np<=0) {
        log_exit("trackRfCavityWithWakes");
        return(np);
        }
#else
    if (notSinglePart) {
      if (isMaster)
	np_tmp = 0;
      else
	np_tmp = np;
      MPI_Allreduce(&np_tmp, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);   
      if (np_total<=0) {
	log_exit("trackRfCavityWithWakes");
	return(np_total);
      }
    }
    else if (np<=0) {
        log_exit("trackRfCavityWithWakes");
        return(np);
    }
#endif 
    if (rfca->change_t && rfca->Q)
        bombElegant("incompatible RF cavity parameters: change_t!=0 && Q!=0", NULL);

    length = rfca->length;

    if (rfca->volt==0 && !wake && !trwake && !LSCKick) {
      if (rfca->length) {
	if (isSlave || !notSinglePart) {
	  for (ip=0; ip<np; ip++) {
	    coord = part[ip];
	    coord[0] += coord[1]*length;
	    coord[2] += coord[3]*length;
	    coord[4] += length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
	  }
	}
      }
      log_exit("trackRfCavityWithWakes");
      return(np);
    }

    omega = PIx2*rfca->freq;
    volt  = rfca->volt/(1e6*particleMassMV*particleRelSign);
    if (omega)
        tau = rfca->Q/omega;
    else
        tau = 0;

    nKicks = length?rfca->nKicks:1;
    if (nKicks<=0) {
      matrixMethod = 1;
      nKicks = 1;
    }
    length /= nKicks;
    volt /= nKicks;
    dtLight = length/c_mks;
    
    if (rfca->phase_reference==0) {
      rfca->phase_reference = unused_phase_reference();
#if defined(DEBUG)
      fprintf(stdout, "RFCA assigned to phase reference %ld\n", rfca->phase_reference);
#endif
    }
    

    switch (get_phase_reference(&phase, rfca->phase_reference)) {
        case REF_PHASE_RETURNED:
            break;
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
            if (!rfca->fiducial_seen) {
                unsigned long mode;
                if (!(mode = parseFiducialMode(rfca->fiducial)))
                    bombElegant("invalid fiducial mode for RFCA element", NULL);
                if (rfca->tReference!=-1)
                  t0 = rfca->tReference;
                else 
                  t0 = findFiducialTime(part, np, zEnd-rfca->length, length/2, *P_central, mode);
                rfca->phase_fiducial = -omega*t0;
                rfca->fiducial_seen = 1;
                }
            set_phase_reference(rfca->phase_reference, phase=rfca->phase_fiducial);
#if defined(DEBUG)
            fprintf(stdout, "RFCA fiducial phase is %e\n", phase);
#endif
            break;
        default:
            bombElegant("unknown return value from get_phase_reference()", NULL);
            break;
        }

    if (omega) {
        t0 = -rfca->phase_fiducial/omega;
        To = PIx2/omega;
        }
    else
        t0 = To = 0;
    phase += rfca->phase*PI/180.0;

    same_dgamma = 0;
    if (omega==0 && tau==0) {
        dgamma = volt*sin(phase);
        dgammaMax = volt;
        same_dgamma = 1;
        }

    timeOffset = 0;
    if (isSlave || !notSinglePart) {
      if (omega && rfca->change_t) {
	coord = part[0];
	P     = *P_central*(1+coord[5]);
	beta_i = P/(gamma=sqrt(sqr(P)+1));
	t     = coord[4]/(c_mks*beta_i);
	if (omega!=0 && t>(0.9*To) && rfca->change_t)
	  timeOffset = ((long)(t/To+0.5))*To;
      }
    }

    linearize = rfca->linearize;
    lockPhase = rfca->lockPhase;

    if (linearize || lockPhase) {
      tAve = 0;
      if (nKicks!=1)
        bombElegant("Must use n_kicks=1 for linearized rf cavity", NULL);

      if (isSlave || !notSinglePart) {
	for (ip=0; ip<np; ip++) {
	  coord = part[ip];
	  P     = *P_central*(1+coord[5]);
	  beta_i = P/(gamma=sqrt(sqr(P)+1));
	  t     = (coord[4]+length/2)/(c_mks*beta_i)-timeOffset;
#ifndef USE_KAHAN	  
	  tAve += t;
#else
          tAve = KahanPlus(tAve, t, &error); 
#endif
	}
      }
#if (!USE_MPI)
      tAve /= np;
#else
      if (notSinglePart) {
	if (USE_MPI) {
	  double tAve_total = 0.0;
#ifndef USE_KAHAN
	  MPI_Allreduce(&tAve, &tAve_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  tAve = tAve_total/np_total;
#else
	  sumArray = malloc(sizeof(double) * n_processors);
	  errorArray = malloc(sizeof(double) * n_processors);

	  MPI_Allgather(&tAve,1,MPI_DOUBLE,sumArray,1,MPI_DOUBLE,MPI_COMM_WORLD);
	  /* collect errors from all processors */
	  MPI_Allgather(&error,1,MPI_DOUBLE,errorArray,1,MPI_DOUBLE,MPI_COMM_WORLD);
	  for (i=1; i<n_processors; i++) {
	    tAve_total = KahanPlus(tAve_total, sumArray[i], &error_sum);
	    tAve_total = KahanPlus(tAve_total, errorArray[i], &error_sum);
	  }
	  tAve = tAve_total/np_total;
	  
	  free(sumArray);
	  free(errorArray);
#endif
	}
      }
      else
        tAve /= np;
#endif
      if (lockPhase) {
	phase = PI/180*rfca->phase;
	dgammaAve = volt*sin(phase);
      }
      else 
	dgammaAve = volt*sin(omega*(tAve-timeOffset)+phase);
    }
    if (isSlave || !notSinglePart) {
      for (ip=0; ip<np; ip++) {
	coord = part[ip];
	coord[0] -= rfca->dx;
	coord[2] -= rfca->dy;
      }
    }

    if (!matrixMethod) {
      double *inverseF;
      inverseF = calloc(sizeof(*inverseF), np);
      
      for (ik=0; ik<nKicks; ik++) {
        dgammaOverGammaAve = dgammaOverGammaNp = 0;
	if (isSlave || !notSinglePart) {
	  for (ip=0; ip<np; ip++) {
	    coord = part[ip];
	    if (coord[5]==-1)
	      continue;
	    if (length)
	      /* compute distance traveled to center of this section */
	      dc4 = length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
	    else 
	      dc4 = 0;
          
	    /* compute energy kick */
	    P     = *P_central*(1+coord[5]);
	    beta_i = P/(gamma=sqrt(sqr(P)+1));
	    t     = (coord[4]+dc4)/(c_mks*beta_i) - timeOffset;
	    if (ik==0 && timeOffset && rfca->change_t) 
	      coord[4] = t*c_mks*beta_i-dc4;
	    if ((dt = t-t0)<0)
	      dt = 0;
	    if  (!same_dgamma) {
	      if (!linearize)
		dgamma = volt*sin(omega*(t-(lockPhase?tAve:0)-ik*dtLight)+phase)*(tau?sqrt(1-exp(-dt/tau)):1);
	      else
		dgamma = dgammaAve +  volt*omega*(t-tAve)*cos(omega*(tAve-timeOffset)+phase);
	    }
	    if (gamma) {
	      dgammaOverGammaNp ++;
	      dgammaOverGammaAve += dgamma/gamma;
	    }
          
	    if (length) {
	      if (rfca->end1Focus && ik==0) {
		/* apply focus kick */
                inverseF[ip] = dgamma/(2*gamma*length);
		coord[1] -= coord[0]*inverseF[ip];
		coord[3] -= coord[2]*inverseF[ip];
	      } 
	      /* apply initial drift */
	      coord[0] += coord[1]*length/2;
	      coord[2] += coord[3]*length/2;
	      coord[4] += length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
	    } 

	    /* apply energy kick */
	    add_to_particle_energy(coord, t, *P_central, dgamma);
	    if ((gamma1 = gamma+dgamma)<=1)
	      coord[5] = -1;
	    else 
	      /* compute inverse focal length for exit kick */
	      inverseF[ip] = -dgamma/(2*gamma1*length);
	  }
	}
        if (!wakesAtEnd) {
          /* do wakes */
          if (wake) 
            track_through_wake(part, np, wake, P_central, run, iPass, charge);
          if (trwake)
            track_through_trwake(part, np, trwake, *P_central, run, iPass, charge);
          if (LSCKick) {
#if !USE_MPI
            if (dgammaOverGammaNp)
              dgammaOverGammaAve /= dgammaOverGammaNp;           
#else
	    if (notSinglePart) {
              double t1 = dgammaOverGammaAve;
              long t2 = dgammaOverGammaNp;
	      MPI_Allreduce (&t1, &dgammaOverGammaAve, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
              MPI_Allreduce (&t2, &dgammaOverGammaNp, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	      if (dgammaOverGammaNp)
		dgammaOverGammaAve /= dgammaOverGammaNp; 
	    } else if (dgammaOverGammaNp)
	      dgammaOverGammaAve /= dgammaOverGammaNp;
    
#endif
            addLSCKick(part, np, LSCKick, *P_central, charge, length, dgammaOverGammaAve);
          }
        }
        if (length) {
	  if(isSlave || !notSinglePart) {
	    /* apply final drift and focus kick if needed */
	    for (ip=0; ip<np; ip++) {
	      coord = part[ip];
	      coord[0] += coord[1]*length/2;
	      coord[2] += coord[3]*length/2;
	      coord[4] += length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
	      if (rfca->end2Focus && (ik==nKicks-1)) {
		coord[1] -= coord[0]*inverseF[ip];
		coord[3] -= coord[2]*inverseF[ip];
	      }
	    }
	  }
        }
        
        if (wakesAtEnd) {
          /* do wakes */
          if (wake) 
            track_through_wake(part, np, wake, P_central, run, iPass, charge);
          if (trwake)
            track_through_trwake(part, np, trwake, *P_central, run, iPass, charge);
          if (LSCKick) {
            if (dgammaOverGammaNp)
#if !USE_MPI
              dgammaOverGammaAve /= dgammaOverGammaNp;           
#else
	    if (notSinglePart) {
              double t1 = dgammaOverGammaAve;
              long t2 = dgammaOverGammaNp;
	      MPI_Allreduce (&t1, &dgammaOverGammaAve, 1, MPI_DOUBLE, MPI_SUM, workers);
              MPI_Allreduce (&t2, &dgammaOverGammaNp, 1, MPI_LONG, MPI_SUM, workers);  
              dgammaOverGammaAve /= dgammaOverGammaNp; 
	    } else
	      dgammaOverGammaAve /= dgammaOverGammaNp;
#endif
            addLSCKick(part, np, LSCKick, *P_central, charge, length, dgammaOverGammaAve);
          }
        }
      }
      free(inverseF);
    } else {
      double sin_phase=0.0, cos_phase, inverseF;
      double R11=1, R21=0, R22, R12, dP, ds1;

      for (ik=0; ik<nKicks; ik++) {
	dgammaOverGammaAve = dgammaOverGammaNp = 0;
	if (isSlave || !notSinglePart) {
          for (ip=0; ip<np; ip++) {
            coord = part[ip];
        
            /* use matrix to propagate particles */
            /* compute energy change using phase of arrival at center of cavity */
            P     = *P_central*(1+coord[5]);
            beta_i = P/(gamma=sqrt(sqr(P)+1));
            ds1 = length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            t     = (coord[4]+ds1)/(c_mks*beta_i)-timeOffset;
            if (timeOffset && rfca->change_t) 
              coord[4] = t*c_mks*beta_i-ds1;
            if ((dt = t-t0)<0)
              dt = 0;
            if  (!same_dgamma) {
              if (!linearize) {
                sin_phase = sin(omega*(t-(lockPhase?tAve:0)-ik*dtLight)+phase);
                cos_phase = cos(omega*(t-(lockPhase?tAve:0)-ik*dtLight)+phase);
                dgamma = (dgammaMax=volt*(tau?sqrt(1-exp(-dt/tau)):1))*sin_phase;
              } else {
                cos_phase = cos(omega*tAve+phase);
                sin_phase = omega*(t-tAve)*cos_phase;
                dgamma = (dgammaMax=volt*(tau?sqrt(1-exp(-dt/tau)):1))*sin_phase +
                  dgammaAve;
              }
            }
            
            if (rfca->end1Focus && length && ik==0) {
              /* apply end focus kick */
              inverseF = dgamma/(2*gamma*length);
              coord[1] -= coord[0]*inverseF;
              coord[3] -= coord[2]*inverseF;
            } 
            
            dP = sqrt(sqr(gamma+dgamma)-1) - P;
	    if (gamma) {
	      dgammaOverGammaNp ++;
	      dgammaOverGammaAve += dgamma/gamma;
	    }
            
            if (useSRSModel) {
              /* note that Rosenzweig and Serafini use gamma in places
               * where they should probably use momentum, but I'll keep
               * their expressions for now.
               */
              double alpha, sin_alpha, gammaf;
              gammaf = gamma+dgamma;
              if (fabs(sin_phase)>1e-6)
                alpha = log(gammaf/gamma)/(2*SQRT2*sin_phase);
              else
                alpha = dgammaMax/gamma/(2*SQRT2);
              R11 = cos(alpha);
              R22 = R11*gamma/gammaf;
              R12 = 2*SQRT2*gamma*length/dgammaMax*(sin_alpha=sin(alpha));
              R21 = -sin_alpha*dgammaMax/(length*gammaf*2*SQRT2);
            } else {
              /* my original treatment used momentum for all 
               * computations, which I still think is correct
               */
              R22 = 1/(1+dP/P);
              if (fabs(dP/P)>1e-14)
                R12 = length*(P/dP*log(1+dP/P));
              else
                R12 = length;
            }
            
            coord[4] += ds1;
            x = coord[0];
            xp = coord[1];
            coord[0] = x*R11 + xp*R12;
            coord[1] = x*R21 + xp*R22;
            x = coord[2];
            xp = coord[3];
            coord[2] = x*R11 + xp*R12;
            coord[3] = x*R21 + xp*R22;
            coord[4] += length/2*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[5] = (P+dP-(*P_central))/(*P_central);
            
            if ((gamma += dgamma)<=1)
              coord[5] = -1;
            if (rfca->end2Focus && length && ik==(nKicks-1)) {
              inverseF = -dgamma/(2*gamma*length);
              coord[1] -= coord[0]*inverseF;
              coord[3] -= coord[2]*inverseF;
            }
            /* adjust s for the new particle velocity */
            coord[4] = (P+dP)/gamma*coord[4]/beta_i;
          }
        }
        
        /* do wakes */
        if (wake) 
          track_through_wake(part, np, wake, P_central, run, iPass, charge);
        if (trwake)
          track_through_trwake(part, np, trwake, *P_central, run, iPass, charge);
        if (LSCKick) {
          if (dgammaOverGammaNp)
#if !USE_MPI
            dgammaOverGammaAve /= dgammaOverGammaNp;           
#else
          if (USE_MPI) {
            double t1 = dgammaOverGammaAve;
            long t2 = dgammaOverGammaNp;
            MPI_Allreduce (&t1, &dgammaOverGammaAve, 1, MPI_DOUBLE, MPI_SUM, workers);
            MPI_Allreduce (&t2, &dgammaOverGammaNp, 1, MPI_LONG, MPI_SUM, workers);  
            dgammaOverGammaAve /= dgammaOverGammaNp; 
          }
#endif
          addLSCKick(part, np, LSCKick, *P_central, charge, length, dgammaOverGammaAve);
        }
      }     
    }
    
    if (isSlave || !notSinglePart) {
      for (ip=0; ip<np; ip++) {
	coord = part[ip];
	coord[0] += rfca->dx;
	coord[2] += rfca->dy;
      }    
      np = removeInvalidParticles(part, np, accepted, zEnd, *P_central);
    }
    if (rfca->change_p0)
      do_match_energy(part, np, P_central, 0);

    return(np);
}    

void add_to_particle_energy(double *coord, double timeOfFlight, double Po, double dgamma)
{
  double gamma, gamma1, PRatio, P, P1, Pz1, Pz;

  P = Po*(1+coord[5]);                    /* old momentum */
  gamma1 = (gamma=sqrt(P*P+1)) + dgamma;  /* new gamma */
  if (gamma1<=1)
    gamma1 = 1+1e-7;
  P1 = sqrt(gamma1*gamma1-1);             /* new momentum */
  coord[5] = (P1-Po)/Po;                  

  /* adjust s for the new particle velocity */
  coord[4] = timeOfFlight*c_mks*P1/gamma1;

  /* adjust slopes so that Px and Py are conserved */
  Pz = P/sqrt(1+sqr(coord[1])+sqr(coord[3]));
  Pz1 = sqrt(Pz*Pz + gamma1*gamma1 - gamma*gamma);
  PRatio = Pz/Pz1;
  coord[1] *= PRatio;
  coord[3] *= PRatio;

}

long track_through_rfcw
  (double **part, long np, RFCW *rfcw, double **accepted, double *P_central, double zEnd,
   RUN *run, long i_pass, CHARGE *charge
   )
{
  static long warned = 0;
  if (rfcw->cellLength<=0) 
    bombElegant("invalid cell length for RFCW", NULL);
  if (rfcw->length==0 && !warned) {
    fprintf(stdout, "** Warning: length of RFCW element is zero. Wakefields will scale to 0!\n");
    warned = 1;
  }
  /* set up the RFCA, TRWAKE, and WAKE structures */
  rfcw->rfca.length = rfcw->length;
  rfcw->rfca.volt = rfcw->volt;
  rfcw->rfca.phase  = rfcw->phase ;
  rfcw->rfca.freq = rfcw->freq;
  rfcw->rfca.Q = rfcw->Q;
  rfcw->rfca.change_p0 = rfcw->change_p0;
  rfcw->rfca.change_t = rfcw->change_t;
  rfcw->rfca.tReference = -1;
  rfcw->rfca.end1Focus = rfcw->end1Focus;
  rfcw->rfca.end2Focus = rfcw->end2Focus;
  if (rfcw->bodyFocusModel) 
    SDDS_CopyString(&rfcw->rfca.bodyFocusModel, rfcw->bodyFocusModel);
  else
    rfcw->rfca.bodyFocusModel = NULL;
  rfcw->rfca.nKicks = rfcw->length?rfcw->nKicks:1;
  rfcw->rfca.dx = rfcw->dx;
  rfcw->rfca.dy = rfcw->dy;
  rfcw->rfca.linearize = rfcw->linearize;
  if (!rfcw->initialized) {
    rfcw->rfca.phase_reference = rfcw->phase_reference;
    if (rfcw->fiducial)
      SDDS_CopyString(&rfcw->rfca.fiducial, rfcw->fiducial);
    else 
      rfcw->rfca.fiducial = NULL;
  }

  rfcw->trwake.charge = 0;
  rfcw->trwake.xfactor = rfcw->trwake.yfactor = 1;
  rfcw->trwake.n_bins = rfcw->n_bins;
  rfcw->trwake.interpolate = rfcw->interpolate;
  rfcw->trwake.smoothing = rfcw->smoothing;
  rfcw->trwake.SGHalfWidth = rfcw->SGHalfWidth;
  rfcw->trwake.SGOrder = rfcw->SGOrder;
  /* misalignment is taken care of by code before and after wake call */
  rfcw->trwake.dx = 0;
  rfcw->trwake.dy = 0;
  rfcw->trwake.xDriveExponent = rfcw->trwake.yDriveExponent = 1;
  rfcw->trwake.xProbeExponent = rfcw->trwake.yProbeExponent = 0;
  if (!rfcw->initialized && rfcw->includeTrWake) {
    rfcw->trwake.initialized = 0;
    if (rfcw->wakeFile) {
      if (rfcw->trWakeFile || rfcw->zWakeFile)
        SDDS_Bomb("You can't give wakeFile along with trWakeFile or zWakeFile for RFCW element");
      SDDS_CopyString(&rfcw->trWakeFile, rfcw->wakeFile);
      SDDS_CopyString(&rfcw->zWakeFile, rfcw->wakeFile);
    }
    
    if (rfcw->WxColumn || rfcw->WyColumn) {
      if (!rfcw->trWakeFile)
        SDDS_Bomb("no input file for transverse wake for RFCW element");
      SDDS_CopyString(&rfcw->trwake.inputFile, rfcw->trWakeFile);
      if (!rfcw->tColumn)
        SDDS_Bomb("no tColumn value for wake for RFCW element");
      SDDS_CopyString(&rfcw->trwake.tColumn, rfcw->tColumn);
      if (rfcw->WxColumn) 
        SDDS_CopyString(&rfcw->trwake.WxColumn, rfcw->WxColumn);
      if (rfcw->WyColumn)
        SDDS_CopyString(&rfcw->trwake.WyColumn, rfcw->WyColumn);
    } else
      rfcw->WxColumn = rfcw->WyColumn = NULL;
  }
  
  rfcw->wake.charge = 0;
  rfcw->wake.n_bins = rfcw->n_bins;
  rfcw->wake.interpolate = rfcw->interpolate;
  rfcw->wake.smoothing = rfcw->smoothing;
  rfcw->wake.SGHalfWidth = rfcw->SGHalfWidth;
  rfcw->wake.SGOrder = rfcw->SGOrder;
  rfcw->wake.change_p0 = rfcw->change_p0;
  if (!rfcw->initialized && rfcw->includeZWake) {
    if (rfcw->WzColumn) {
      if (!rfcw->zWakeFile)
        SDDS_Bomb("no input file for z wake for RFCW element");
      SDDS_CopyString(&rfcw->wake.inputFile, rfcw->zWakeFile);
      if (!rfcw->tColumn)
        SDDS_Bomb("no tColumn value for wake for RFCW element");
      SDDS_CopyString(&rfcw->wake.tColumn, rfcw->tColumn);
      SDDS_CopyString(&rfcw->wake.WColumn, rfcw->WzColumn);
      rfcw->wake.initialized = 0;
    } else
      rfcw->wake.WColumn = NULL;
  }

  rfcw->LSCKick.bins = rfcw->LSCBins;
  rfcw->LSCKick.interpolate = rfcw->LSCInterpolate;
  rfcw->LSCKick.lowFrequencyCutoff0 = rfcw->LSCLowFrequencyCutoff0;
  rfcw->LSCKick.lowFrequencyCutoff1 = rfcw->LSCLowFrequencyCutoff1;
  rfcw->LSCKick.highFrequencyCutoff0 = rfcw->LSCHighFrequencyCutoff0;
  rfcw->LSCKick.highFrequencyCutoff1 = rfcw->LSCHighFrequencyCutoff1;
  rfcw->LSCKick.radiusFactor = rfcw->LSCRadiusFactor;
  
  rfcw->initialized = 1;

  if (rfcw->WzColumn && rfcw->includeZWake)
    rfcw->wake.factor = rfcw->length/rfcw->cellLength/(rfcw->rfca.nKicks?rfcw->rfca.nKicks:1);
  if ((rfcw->WxColumn || rfcw->WyColumn) && rfcw->includeTrWake)
    rfcw->trwake.factor = rfcw->length/rfcw->cellLength/(rfcw->rfca.nKicks?rfcw->rfca.nKicks:1);
  np = trackRfCavityWithWakes(part, np, &rfcw->rfca, accepted, P_central, zEnd,
                              i_pass, run, charge, 
                              (rfcw->WzColumn && rfcw->includeZWake) ? &rfcw->wake : NULL,
                              ((rfcw->WxColumn || rfcw->WyColumn) && rfcw->includeTrWake) ? &rfcw->trwake : NULL,
                              rfcw->doLSC ? &rfcw->LSCKick : NULL,
                              rfcw->wakesAtEnd);
  return np;
}



/* See H. Wiedemann, Particle Accelerator Physics I, 8.2.2 */
double rfAcceptance_Fq(double q) 
{
  return 2*(sqrt(q*q-1)-acos(1/q));
}

double solveForOverVoltage(double F, double q0)
{
  return zeroNewton(&rfAcceptance_Fq, F, q0, 1e-6, 1000, 1e-12);
}

