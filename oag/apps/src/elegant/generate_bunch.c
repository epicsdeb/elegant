/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: generate_bunch()
 * purpose: generate a 6-dimensional bunch in a variety of distributions
 * The distributions are calculated first in normalized coordinates, then transformed:
 *     u = x/sqrt(beta)    u' = x'/sqrt(beta) + alpha/beta*u
 *
 *     psi(x, x') = N*exp(-(x'^2*beta + 2*alpha*x'*x + x^2*gamma)/(2*eps)) = 
 *                  N*exp(-(u'^2*beta^2 + u^2)/(2*eps))
 *     hence, sigma(u) = sqrt(eps), sigma(u') = sigma(u)/beta
 * 
 * The limits for the gaussian distribution are usually put on the sigma values for u' and u,
 * but can also by put on the invariant.
 *
 * Momentum effects are included via the equations (for x and y):
 *        x -> x + dp*eta      x' -> x' + dp*eta'
 *
 * The longitudinal distribution can be rotated to give s-dp correlation.  The "coupling angle" is used as in:
 *       s -> s + dp*tan(angle)
 *
 * Michael Borland, 1989, 1991
 */
#include "mdb.h"
#include "track.h"

void enforce_beta_alpha_emit(double **coord, long n_part, long offset, double beta, double alpha, double emit);
void zero_centroid(double **particle, long n_particles, long coord);
long dynap_distribution(double **particle, long n_particles, double sx, double sy,
            long nx, long ny);
#if SDDS_MPI_IO
long dynap_distribution_p(double **particle, long n_particles, double sx, double sy,
            long nx, long ny);
#endif
void gaussian_distribution(double **particle, long n_particles, 
    long offset, double s1, double s2, long symmetrize, long *haltonID, double limit,
    double limit_invar, double beta, long halo);
void hard_edge_distribution(double **particle, long n_particles, 
    long offset, double max1, double max2, long symmetrize, long *haltonID, double cutoff);
void uniform_distribution(double **particle, long n_particles, long offset, 
    double max1, double max2, long symmetrize, long *haltonID, double cutoff);
void shell_distribution(double **particle, long n_particles, long offset, 
    double s1, double s2, long symmetrize, double cutoff);
void line_distribution(double **particle, long n_particles, long offset, 
    double max1, double max2);
void transform_from_normalized_coordinates(double **part, long n_part, long offset, double beta, double alpha);
void enforce_sigma_values(double **coord, long n_part, long offset, double s1d, double s2d);
void couple_coordinates(double **coord, long n_part, long offset, double angle);
void gaussian_4d_distribution(double **particle, long n_particles, long offset, 
    double s1, double s2, double s3, double s4, double cutoff, long halo);
void uniform_4d_distribution(double **particle, long n_particles, long offset, 
    double max1, double max2, double max3, double max4, double cutoff);

long generate_bunch(
    double **particle,
    long n_particles,
    TRANSVERSE *x_plane, 
    TRANSVERSE *y_plane,
    LONGITUDINAL *longit,
    long *enforce_rms_params,
    long limit_invar,
    long symmetrize,
    long *haltonID,
    long *doRandomizeOrder,
    long limit_in_4d,
    double Po
    )
{
    long i_particle, first_call=1, i, j;
    double s1, s2, s3, s4, delta_p;
    double s56, beta, emit, alpha=0.0;
    double *randomizedData = NULL;
#if !SDDS_MPI_IO
    /* for Pelegant regression test */
    static long initial_saved=0;
    static double **first_particle_address;  
    static long total_n_particles = 0;
    
    if (!initial_saved) { 
      first_particle_address = particle;  /* Save the beginning address of particle array */
      remaining_sequence_No = orig_sequence_No;
      initial_saved = 1;
    }
    else
      remaining_sequence_No--;
    total_n_particles += n_particles;
#else
    long sum=0, tmp, my_offset, *offset = tmalloc(n_processors*sizeof(*offset)), total_particles=0;
    if (isSlave) {
      MPI_Allgather (&n_particles, 1, MPI_LONG, offset, 1, MPI_LONG, workers);
      tmp = offset[0];
      for (i=1; i<n_processors; i++) {
	sum += tmp;
	tmp = offset[i];
	offset[i] = sum; 
      }
      offset[0] = 0; 
      my_offset = offset[myid-1];
      particleID += my_offset;
      total_particles =  sum;
    }
#endif
    if (x_plane->beam_type==DYNAP_BEAM) {
      if (first_call) {
#if SDDS_MPI_IO
        n_particles = dynap_distribution_p(particle, n_particles,
                                         sqrt(x_plane->emit*x_plane->beta),
                                         sqrt(y_plane->emit*y_plane->beta),
                                         (long)(x_plane->cutoff+0.5),
                                         (long)(y_plane->cutoff+0.5));
#else
	n_particles = total_n_particles;
	particle = first_particle_address;
	remaining_sequence_No=1; /*This is a special case of generating particles in parallel */
        n_particles = dynap_distribution(particle, n_particles, 
                                         sqrt(x_plane->emit*x_plane->beta),        
                                         sqrt(y_plane->emit*y_plane->beta),
                                         (long)(x_plane->cutoff+0.5),
                                         (long)(y_plane->cutoff+0.5));
#endif
        set_beam_centroids(particle, 0, n_particles, x_plane->cent_posi, x_plane->cent_slope);
        set_beam_centroids(particle, 2, n_particles, y_plane->cent_posi, y_plane->cent_slope);
        set_beam_centroids(particle, 4, n_particles, longit->cent_s, longit->cent_dp);
	for (i_particle=0; i_particle<n_particles; i_particle++) 
	  particle[i_particle][6] = particleID++;
      }
#if SDDS_MPI_IO
      /* prepare for the next bunch */
      particleID += (total_particles - n_particles - my_offset + 1);
#endif
      first_call = 0;
    }
    else {
      if (!limit_in_4d) {
        /* make mono-energetic distribution in x plane */
        s1 = sqrt(x_plane->emit);
        s2 = s1/x_plane->beta;
	switch (x_plane->beam_type) {
        case GAUSSIAN_BEAM:
          gaussian_distribution(particle, n_particles, 0, s1, s2, 
                                symmetrize && (enforce_rms_params[0] || n_particles%4==0), 
                                haltonID, x_plane->cutoff,
                                (limit_invar?sqr(x_plane->cutoff)*x_plane->emit:0.0), x_plane->beta, 
                                0);
          break;
        case GAUSSIAN_HALO_BEAM:
          gaussian_distribution(particle, n_particles, 0, s1, s2, 
                                symmetrize && (enforce_rms_params[0] || n_particles%4==0), 
                                haltonID, x_plane->cutoff,
                                (limit_invar?sqr(x_plane->cutoff)*x_plane->emit:0.0), x_plane->beta,
                                1);
          break;
        case HARD_EDGE_BEAM:
          hard_edge_distribution(particle, n_particles, 0, s1, s2, 
                                 symmetrize && (enforce_rms_params[0] || n_particles%4==0), 
                                 haltonID, x_plane->cutoff);
          break;
        case UNIFORM_ELLIPSE:
          uniform_distribution(particle, n_particles, 0, s1, s2, 
                               symmetrize && (enforce_rms_params[0] || n_particles%4==0), 
                               haltonID, x_plane->cutoff);
          break;
        case SHELL_BEAM:
          shell_distribution(particle, n_particles, 0, s1, s2, 
                             symmetrize && (enforce_rms_params[0] || n_particles%4==0), x_plane->cutoff);
          break;
        case LINE_BEAM:
          line_distribution(particle, n_particles, 0, s1, s2);
          break;
        default:
          printf("fatal internal error: unknown x beam-type %ld\n",
                 x_plane->beam_type);
          abort();
          break;
        }
#if !SDDS_MPI_IO
	if(remaining_sequence_No<=1) {
	  zero_centroid(first_particle_address, total_n_particles, 0);
	  zero_centroid(first_particle_address, total_n_particles, 1);
	  if (enforce_rms_params[0])
	    enforce_sigma_values(first_particle_address, total_n_particles, 0, s1, s2);
	  transform_from_normalized_coordinates(first_particle_address, total_n_particles, 0, x_plane->beta, x_plane->alpha);
	}
#else
        zero_centroid(particle, n_particles, 0);
        zero_centroid(particle, n_particles, 1);
        if (enforce_rms_params[0])
          enforce_sigma_values(particle, n_particles, 0, s1, s2);
        transform_from_normalized_coordinates(particle, n_particles, 0, x_plane->beta, x_plane->alpha);
#endif
        
        /* make mono-energetic distribution in y plane */
        s1 = sqrt(y_plane->emit);
        s2 = s1/y_plane->beta;
        switch (y_plane->beam_type) {
        case GAUSSIAN_BEAM:
          gaussian_distribution(particle, n_particles, 2, s1, s2, 
                                symmetrize && (enforce_rms_params[1] || n_particles%4==0), 
                                haltonID,
                                y_plane->cutoff, (limit_invar?sqr(y_plane->cutoff)*y_plane->emit:0.0),
                                y_plane->beta, 0);
          break;
        case GAUSSIAN_HALO_BEAM:
          gaussian_distribution(particle, n_particles, 2, s1, s2, 
                                symmetrize && (enforce_rms_params[1] || n_particles%4==0), 
                                haltonID,
                                y_plane->cutoff, (limit_invar?sqr(y_plane->cutoff)*y_plane->emit:0.0),
                                y_plane->beta, 1);
          break;
        case HARD_EDGE_BEAM:
          hard_edge_distribution(particle, n_particles, 2, s1, s2, 
                                 symmetrize && (enforce_rms_params[1] || n_particles%4==0), 
                                 haltonID, y_plane->cutoff);
          break;
        case UNIFORM_ELLIPSE:
          uniform_distribution(particle, n_particles, 2, s1, s2, 
                               symmetrize && (enforce_rms_params[1] || n_particles%4==0), 
                               haltonID, y_plane->cutoff);
          break;
        case SHELL_BEAM:
          shell_distribution(particle, n_particles, 2, s1, s2, 
                             symmetrize && (enforce_rms_params[1] || n_particles%4==0), y_plane->cutoff);
          break;
        case LINE_BEAM:
          line_distribution(particle, n_particles, 2, s1, s2);
          break;
        default:
          printf("fatal internal error: unknown y beam-type %ld\n",
                 y_plane->beam_type);
          abort();
          break;
        }
#if !SDDS_MPI_IO
        if(remaining_sequence_No<=1) {
          zero_centroid(first_particle_address, total_n_particles, 2);
          zero_centroid(first_particle_address, total_n_particles, 3);
          if (enforce_rms_params[1])
            enforce_sigma_values(first_particle_address, total_n_particles, 2, s1, s2);
          transform_from_normalized_coordinates(first_particle_address, total_n_particles, 2, y_plane->beta, y_plane->alpha);
        }
#else
        zero_centroid(particle, n_particles, 2);
        zero_centroid(particle, n_particles, 3);
        if (enforce_rms_params[1])
          enforce_sigma_values(particle, n_particles, 2, s1, s2);
        transform_from_normalized_coordinates(particle, n_particles, 2, y_plane->beta, y_plane->alpha);
#endif
      }
      else {
        if (y_plane->beam_type!=x_plane->beam_type || y_plane->cutoff!=x_plane->cutoff)
          bomb("distribution types and cutoffs for x and y planes must be the same for limit_in_4d", NULL);
        s1 = sqrt(x_plane->emit);
        s2 = s1/x_plane->beta;
        s3 = sqrt(y_plane->emit);
        s4 = s3/y_plane->beta;
        switch (x_plane->beam_type) {
        case GAUSSIAN_BEAM:
          gaussian_4d_distribution(particle, n_particles, 0, s1, s2, s3, s4, x_plane->cutoff, 0);
          break;
        case GAUSSIAN_HALO_BEAM:
          gaussian_4d_distribution(particle, n_particles, 0, s1, s2, s3, s4, x_plane->cutoff, 1);
          break;
        case UNIFORM_ELLIPSE:
          uniform_4d_distribution(particle, n_particles, 0, s1, s2, s3, s4, x_plane->cutoff);
          break;
        default:
          bomb("limit_in_4d is available only for gaussian and uniform beam distributions", NULL);
          break;
        }
#if !SDDS_MPI_IO
        if(remaining_sequence_No<=1) {
          zero_centroid(first_particle_address, total_n_particles, 0);
          zero_centroid(first_particle_address, total_n_particles, 1);
          zero_centroid(first_particle_address, total_n_particles, 2);
          zero_centroid(first_particle_address, total_n_particles, 3);
          if (enforce_rms_params[0])
            enforce_sigma_values(first_particle_address, total_n_particles, 0, s1, s2);
          transform_from_normalized_coordinates(first_particle_address, total_n_particles, 0, x_plane->beta, x_plane->alpha);
          if (enforce_rms_params[1])
            enforce_sigma_values(first_particle_address, total_n_particles, 2, s1, s2);
          transform_from_normalized_coordinates(first_particle_address, total_n_particles, 2, x_plane->beta, x_plane->alpha);
        }
#else
        zero_centroid(particle, n_particles, 0);
        zero_centroid(particle, n_particles, 1);
        zero_centroid(particle, n_particles, 2);
        zero_centroid(particle, n_particles, 3);
        if (enforce_rms_params[0])
          enforce_sigma_values(particle, n_particles, 0, s1, s2);
        if (enforce_rms_params[1])
          enforce_sigma_values(particle, n_particles, 2, s3, s4);
        transform_from_normalized_coordinates(particle, n_particles, 0, x_plane->beta, x_plane->alpha);
        transform_from_normalized_coordinates(particle, n_particles, 2, y_plane->beta, y_plane->alpha);
#endif
      }
    }
    
    /* make energy-distance distribution */
    /* three possibilities for inputs:
     * 1. sigma_s, sigma_dp, and dp_s_coupling
     * 2. sigma_s, sigma_dp, and alpha_z
     * 3. emit_z, beta_z, and alpha_z
     */
    if (longit->emit) {
      /* #3 */
      beta = longit->beta;
      alpha = longit->alpha;
      emit = longit->emit;
      s1 = sqrt(emit);
      s2 = s1/beta;
    } else {
      if (!longit->alpha)
        s56 = longit->sigma_dp*longit->sigma_s*longit->dp_s_coupling;
      else
        s56 = -longit->sigma_dp*longit->sigma_s*longit->alpha/sqrt(1+sqr(longit->alpha));
      emit = sqrt(sqr(longit->sigma_dp*longit->sigma_s) - sqr(s56));
#ifdef IEEE_MATH
      if (isnan(emit) || isinf(emit))
        emit = 0;
#endif
      if (emit) {
        beta  = sqr(longit->sigma_s)/emit;
        if (longit->alpha)
          alpha = longit->alpha;
        else
          alpha = -s56/emit;
        s1 = sqrt(emit);
        s2 = s1/beta;
      }
      else {
        s1 = longit->sigma_s;
        s2 = longit->sigma_dp;
        beta = 0;
      }
    }
    switch (longit->beam_type) {
    case GAUSSIAN_BEAM:
      gaussian_distribution(particle, n_particles, 4, s1, s2, 
                            symmetrize && (enforce_rms_params[2] || n_particles%4==0),
                            haltonID, longit->cutoff, 
                            (limit_invar?sqr(longit->cutoff)*emit:0.0), beta, 0);
      break;
    case GAUSSIAN_HALO_BEAM:
      gaussian_distribution(particle, n_particles, 4, s1, s2, 
                            symmetrize && (enforce_rms_params[2] || n_particles%4==0),
                            haltonID, longit->cutoff, 
                            (limit_invar?sqr(longit->cutoff)*emit:0.0), beta, 1);
      break;
    case HARD_EDGE_BEAM:
      hard_edge_distribution(particle, n_particles, 4, s1, s2, 
                             symmetrize && (enforce_rms_params[2] || n_particles%4==0), 
                             haltonID, longit->cutoff);
      break;
    case UNIFORM_ELLIPSE:
      uniform_distribution(particle, n_particles, 4, s1, s2, 
                           symmetrize && (enforce_rms_params[2] || n_particles%4==0),
                           haltonID, longit->cutoff);
      break;
    case SHELL_BEAM:
      shell_distribution(particle, n_particles, 4, s1, s2, 
                         symmetrize && (enforce_rms_params[2] || n_particles%4==0), longit->cutoff);
      break;
    case LINE_BEAM:
      line_distribution(particle, n_particles, 4, s1, s2);
      break;
    default:
      printf("fatal internal error: unknown longitudinal beam-type %ld\n",
             longit->beam_type);
      abort();
      break;
    }
#if !SDDS_MPI_IO
    if(remaining_sequence_No<=1) {
      zero_centroid(first_particle_address, total_n_particles, 4);
      zero_centroid(first_particle_address, total_n_particles, 5);
      if (enforce_rms_params[2])
	enforce_sigma_values(first_particle_address, total_n_particles, 4, s1, s2);
      if (emit)
	transform_from_normalized_coordinates(first_particle_address, total_n_particles, 4, beta, alpha);
    }
#else
    zero_centroid(particle, n_particles, 4);
    zero_centroid(particle, n_particles, 5);
    if (enforce_rms_params[2])
      enforce_sigma_values(particle, n_particles, 4, s1, s2);
    if (emit)
      transform_from_normalized_coordinates(particle, n_particles, 4, beta, alpha);
#endif
    if (x_plane->beam_type==DYNAP_BEAM)
      return(n_particles);

    for (i=0; i<3; i++) {
      switch (doRandomizeOrder[i]) {
      case 1:
        /* randomize two coordinates independently */
        if (!(randomizedData=malloc(sizeof(*randomizedData)*n_particles)))
          bomb("memory allocation failure (generate_bunch)", NULL);
        for (j=0; j<n_particles; j++)
          randomizedData[j] = particle[j][2*i];
        randomizeOrder((char*)randomizedData, sizeof(*randomizedData), n_particles, 0, 
                       random_4);
        for (j=0; j<n_particles; j++)
          particle[j][2*i] = randomizedData[j];

        for (j=0; j<n_particles; j++)
          randomizedData[j] = particle[j][2*i+1];
        randomizeOrder((char*)randomizedData, sizeof(*randomizedData), n_particles, 0, 
                       random_4);
        for (j=0; j<n_particles; j++)
          particle[j][2*i+1] = randomizedData[j];
        free(randomizedData);
        break;
      case 2:
        /* randomize two coordinates together */
        if (!(randomizedData=malloc(sizeof(*randomizedData)*n_particles*2)))
          bomb("memory allocation failure (generate_bunch)", NULL);
        for (j=0; j<n_particles; j++) {
          randomizedData[2*j] = particle[j][2*i];
          randomizedData[2*j+1] = particle[j][2*i+1];
        }
        randomizeOrder((char*)randomizedData, 2*sizeof(*randomizedData), n_particles, 0, 
                       random_4);
        for (j=0; j<n_particles; j++) {
          particle[j][2*i] = randomizedData[2*j];
          particle[j][2*i+1] = randomizedData[2*j+1];
        }
        free(randomizedData);
        break;
      default:
        break;  
      }
    }
    
    /* incorporate dispersion and centroid shifts into (x, x', y, y') */
    /* also add particle ID */
#if !SDDS_MPI_IO
    /* This should be involved after both centroid and sigma have been forced */
  if(remaining_sequence_No<=1) {
    for (i_particle=0; i_particle<total_n_particles; i_particle++) {
      if (longit->chirp)
        first_particle_address[i_particle][5] += longit->chirp*first_particle_address[i_particle][4];
      first_particle_address[i_particle][4] += longit->cent_s;
      first_particle_address[i_particle][5] += longit->cent_dp;
      delta_p = first_particle_address[i_particle][5];
      first_particle_address[i_particle][0] += delta_p*x_plane->eta + x_plane->cent_posi;
      first_particle_address[i_particle][1] += delta_p*x_plane->etap + x_plane->cent_slope;
      first_particle_address[i_particle][2] += delta_p*y_plane->eta + y_plane->cent_posi;
      first_particle_address[i_particle][3] += delta_p*y_plane->etap + y_plane->cent_slope;
      first_particle_address[i_particle][6] = particleID++;
    }
    initial_saved = 0; /* Prepare for next bunch */
    total_n_particles = 0;
  }
#else
    for (i_particle=0; i_particle<n_particles; i_particle++) {
      if (longit->chirp)
        particle[i_particle][5] += longit->chirp*particle[i_particle][4];
      particle[i_particle][4] += longit->cent_s;
      particle[i_particle][5] += longit->cent_dp;
      delta_p = particle[i_particle][5];
      particle[i_particle][0] += delta_p*x_plane->eta + x_plane->cent_posi;
      particle[i_particle][1] += delta_p*x_plane->etap + x_plane->cent_slope;
      particle[i_particle][2] += delta_p*y_plane->eta + y_plane->cent_posi;
      particle[i_particle][3] += delta_p*y_plane->etap + y_plane->cent_slope;
      particle[i_particle][6] = particleID++;
    }
    /* prepare for the next bunch */
    particleID += (total_particles - n_particles - my_offset);
#endif

#ifdef IEEE_MATH
    /* check for invalid data */
    for (i_particle=0; i_particle<n_particles; i_particle++) {
      long i;
      for (i=0; i<6; i++)
        if (isnan(particle[i_particle][i]) || isinf(particle[i_particle][i]))
          bomb("invalid particle data generated in generate_bunch()!", NULL);
    }
#endif

    first_call = 0;
#if SDDS_MPI_IO
    tfree(offset);
#endif
    return(n_particles);
  }

void gaussian_distribution(
                           double **particle, 
                           long n_particles, 
                           long offset, 
                           double s1, 
                           double s2,
                           long symmetrize,
                           long *haltonID,
                           double limit,
                           double limit_invar,
                           double beta,
                           long halo
                           )
{
  double x1, x2, limit1, limit2;
  long i_particle, flag;

  limit1 = s1*limit;
  limit2 = s2*limit;
  if (haltonID[offset] && haltonID[offset+1]) {
    double s12[2], buffer[2];
    long dim;
#if SDDS_MPI_IO
    /* To generate n particles on m processors, each processor will 
       generate myid*n/m particles, but only the last n/m will be used */
    long i, start_particle, *particle_array = tmalloc(n_processors*sizeof(*particle_array));

    MPI_Allgather (&n_particles, 1, MPI_LONG, particle_array, 1, MPI_LONG, MPI_COMM_WORLD);
    for (i=1; i<n_processors; i++) {
      particle_array[i] += particle_array[i-1] ; 
    }
    start_particle = particle_array[myid]-n_particles; /* The first particle for a processor */
    n_particles = particle_array[myid]; 
#endif
    s12[0] = s1;
    s12[1] = s2;
    for (i_particle=0; i_particle<n_particles; i_particle++) {
      do {
        for (dim=0; dim<2; dim++) {
          buffer[dim] = nextHaltonSequencePoint(haltonID[offset+dim]);
          if (!convertSequenceToGaussianDistribution(&buffer[dim], 1, 0)) {
            dim--;
            continue;
          }
          buffer[dim] *= s12[dim];
        }
        if (limit_invar)
          flag = (sqr(buffer[0])+sqr(beta*buffer[1]))<=limit_invar;
        else
          flag = fabs(buffer[0])<=limit1 && fabs(buffer[1])<=limit2;
        if ((!halo && flag) || (halo && !flag))
          break;
      } while (1);
#if SDDS_MPI_IO
      if (i_particle>=start_particle) {
	for (dim=0; dim<2; dim++)
	  particle[i_particle-start_particle][offset+dim] = buffer[dim];
      }
#else
      for (dim=0; dim<2; dim++)
        particle[i_particle][offset+dim] = buffer[dim];
#endif
    }
#if SDDS_MPI_IO
    tfree(particle_array);
#endif
  }
  else {
    for (i_particle=0; i_particle<n_particles; i_particle++) {
      do {
        particle[i_particle][0+offset] = x1 = gauss_rn_lim(0.0, s1, 0.0, random_4);
        particle[i_particle][1+offset] = x2 = gauss_rn_lim(0.0, s2, 0.0, random_4);
        if (limit_invar)
          flag = (sqr(x1)+sqr(beta*x2))<=limit_invar;
        else
          flag = fabs(x1)<=limit1 && fabs(x2)<=limit2;
        if ((!halo && flag) || (halo && !flag))
          break;
      } while (1);
      if (symmetrize) {
        if (++i_particle>=n_particles)
          break;
        particle[i_particle][0+offset] = x1;
        particle[i_particle][1+offset] = -x2;
        if (++i_particle>=n_particles)
          break;
        particle[i_particle][0+offset] = -x1;
        particle[i_particle][1+offset] = x2;
        if (++i_particle>=n_particles)
          break;
        particle[i_particle][0+offset] = -x1;
        particle[i_particle][1+offset] = -x2;
      }
    }
  }
}

void uniform_distribution(
                          double **particle, 
                          long n_particles, 
                          long offset, 
                          double max1, 
                          double max2, 
                          long symmetrize,
                          long *haltonID,
                          double cutoff
                          )
{
  double x1, x2;
  long i_particle;
  double range1, range2; 
  double rnd1, rnd2=0.0;
#if SDDS_MPI_IO
  /* To generate n particles on m processors, each processor will 
     generate myid*n/m particles, but only the last n/m will be used */
  long i, start_particle, *particle_array = tmalloc(n_processors*sizeof(*particle_array));
  
  MPI_Allgather (&n_particles, 1, MPI_LONG, particle_array, 1, MPI_LONG, MPI_COMM_WORLD);
  for (i=1; i<n_processors; i++) {
    particle_array[i] += particle_array[i-1] ; 
  }
  start_particle = particle_array[myid]-n_particles; /* The first particle for a processor */
  n_particles = particle_array[myid]; 
#endif

  log_entry("uniform_distribution");

  range1 = 2*max1*cutoff;
  range2 = 2*max2*cutoff;

  for (i_particle=0; i_particle<n_particles; i_particle++) {
    do {
      if (haltonID[offset] && haltonID[offset+1]) {
        if ((rnd1=nextHaltonSequencePoint(haltonID[offset]))<0 ||
            (rnd2=nextHaltonSequencePoint(haltonID[offset+1]))<0)
          bomb("problem determining Halton sequence", NULL);
        rnd1 -= 0.5;
        rnd2 -= 0.5;
      } else {
        rnd1 = random_4(1)-.5;
        rnd2 = random_4(1)-.5;
      }
    } while (rnd2*rnd2+rnd1*rnd1 > .25) ;
#if !SDDS_MPI_IO
    x1 = range1*rnd1;
    x2 = range2*rnd2;
    particle[i_particle][0+offset] = x1;
    particle[i_particle][1+offset] = x2;
    if (symmetrize) {
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = x1;
      particle[i_particle][1+offset] = -x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = -x1;
      particle[i_particle][1+offset] = x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = -x1;
      particle[i_particle][1+offset] = -x2;
    }
#else
    if (i_particle>=start_particle) {
      x1 = range1*rnd1;
      x2 = range2*rnd2;
      particle[i_particle-start_particle][0+offset] = x1;
      particle[i_particle-start_particle][1+offset] = x2;
      if (symmetrize) {
	if (++i_particle>=n_particles)
	  break;
	particle[i_particle-start_particle][0+offset] = x1;
	particle[i_particle-start_particle][1+offset] = -x2;
	if (++i_particle>=n_particles)
        break;
	particle[i_particle-start_particle][0+offset] = -x1;
	particle[i_particle-start_particle][1+offset] = x2;
	if (++i_particle>=n_particles)
	  break;
	particle[i_particle-start_particle][0+offset] = -x1;
	particle[i_particle-start_particle][1+offset] = -x2;
      }
    }
#endif
  }
#if SDDS_MPI_IO
  tfree(particle_array);
#endif
  log_exit("uniform_distribution");
}

void shell_distribution(
                        double **particle, 
                        long n_particles, 
                        long offset, 
                        double s1, 
                        double s2, 
                        long symmetrize,
                        double cutoff
                        )
{
  double x1, x2;
  long i_particle;
  double angle, dangle;

  log_entry("shell_distribution");

  if (n_particles<=1) 
    dangle = 0;
  else
    dangle = PIx2/(n_particles-1);
  angle = -dangle;

  s1 *= cutoff;
  s2 *= cutoff;
  for (i_particle=0; i_particle<n_particles; i_particle++) {
    x1 = s1*cos(angle += dangle);
    x2 = s2*sin(angle);
    particle[i_particle][0+offset] = x1;
    particle[i_particle][1+offset] = x2;
  }
  log_exit("shell_distribution");
}

void hard_edge_distribution(
                            double **particle, 
                            long n_particles, 
                            long offset, 
                            double max1, 
                            double max2, 
                            long symmetrize,
                            long *haltonID,
                            double cutoff
                            )
{
  double x1, x2;
  long i_particle;
  double range1; 
  double range2; 
  double rnd1, rnd2;
#if SDDS_MPI_IO
    /* To generate n particles on m processors, each processor will 
       generate myid*n/m particles, but only the last n/m will be used */
    long i, start_particle, *particle_array = tmalloc(n_processors*sizeof(*particle_array));

    MPI_Allgather (&n_particles, 1, MPI_LONG, particle_array, 1, MPI_LONG, MPI_COMM_WORLD);
    for (i=1; i<n_processors; i++) {
      particle_array[i] += particle_array[i-1] ; 
    }
    start_particle = particle_array[myid]-n_particles; /* The first particle for a processor */
    n_particles = particle_array[myid]; 
#endif

  log_entry("hard_edge_distribution");

  range1 = 2*max1*cutoff;
  range2 = 2*max2*cutoff;

  for (i_particle=0; i_particle<n_particles; i_particle++) {
    if (haltonID[offset] && haltonID[offset+1]) {
      rnd1 = nextHaltonSequencePoint(haltonID[offset]);
      rnd2 = nextHaltonSequencePoint(haltonID[offset+1]);
      rnd1 -= 0.5;
      rnd2 -= 0.5;
    } else {
      rnd1 = random_4(1)-.5;
      rnd2 = random_4(1)-.5;
    }
#if SDDS_MPI_IO
    if (i_particle<start_particle)
      continue;
    x1 = range1*rnd1;
    x2 = range2*rnd2;
    particle[i_particle-start_particle][0+offset] = x1;
    particle[i_particle-start_particle][1+offset] = x2;
    if (symmetrize) {
      if (++i_particle>=n_particles)
        break;
      particle[i_particle-start_particle][0+offset] = x1;
      particle[i_particle-start_particle][1+offset] = -x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle-start_particle][0+offset] = -x1;
      particle[i_particle-start_particle][1+offset] = x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle-start_particle][0+offset] = -x1;
      particle[i_particle-start_particle][1+offset] = -x2;
    }
#else
    x1 = range1*rnd1;
    x2 = range2*rnd2;
    particle[i_particle][0+offset] = x1;
    particle[i_particle][1+offset] = x2;
    if (symmetrize) {
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = x1;
      particle[i_particle][1+offset] = -x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = -x1;
      particle[i_particle][1+offset] = x2;
      if (++i_particle>=n_particles)
        break;
      particle[i_particle][0+offset] = -x1;
      particle[i_particle][1+offset] = -x2;
    }
#endif
  }
  log_exit("hard_edge_distribution");

#if SDDS_MPI_IO
  tfree(particle_array);
#endif
}


void line_distribution(
    double **particle, 
    long n_particles, 
    long offset, 
    double max1, 
    double max2
    )
{
    long ip;
    double x1, x2, dx1, dx2;

    x1 = -max1;
    x2 = -max2;
    if (n_particles>1) {
        dx1 = 2*max1/(n_particles-1);
        dx2 = 2*max2/(n_particles-1);
        }
    else
        x1 = x2 = dx1 = dx2 = 0;
    for (ip=0; ip<n_particles; ip++) {
        particle[ip][offset] = x1;
        particle[ip][offset+1] = x2;
        x1 += dx1;
        x2 += dx2;
        }
    }


void enforce_sigma_values(double **coord, long n_part, long offset, double s1d, double s2d)
{
    double s1a, s2a;
    double f1, f2;
    long i;
#if SDDS_MPI_IO
    double s1a_total, s2a_total; 
    long total_particles;
#endif

    log_entry("enforce_sigma_values");

#if !SDDS_MPI_IO
    if (n_part<2) {
        log_exit("enforce_sigma_values");
        return;
        }
#endif
    s1a = s2a = 0;
    for (i=0; i<n_part; i++) {
        s1a += sqr(coord[i][offset+0]);
        s2a += sqr(coord[i][offset+1]);
        }
#if !SDDS_MPI_IO
    s1a = sqrt(s1a/n_part);
    s2a = sqrt(s2a/n_part);
#else
    if (isSlave) {
      MPI_Allreduce (&s1a, &s1a_total, 1, MPI_DOUBLE, MPI_SUM, workers);
      MPI_Allreduce (&s2a, &s2a_total, 1, MPI_DOUBLE, MPI_SUM, workers);
      MPI_Allreduce (&n_part, &total_particles, 1, MPI_LONG, MPI_SUM, workers);
      s1a = sqrt(s1a_total/total_particles);
      s2a = sqrt(s2a_total/total_particles);
    }
#endif

    if (s1a!=0)
        f1 = s1d/s1a;
    else
        f1 = 0;
    if (s2a!=0)
        f2 = s2d/s2a;
    else
        f2 = 0;
    for (i=0; i<n_part; i++) {
        coord[i][offset+0] *= f1;
        coord[i][offset+1] *= f2;
        }
    log_exit("enforce_sigma_values");
    }

void zero_centroid(double **particle, long n_particles, long coord)
{
    long i;
    double sum;

    log_entry("zero_centroid");

    if (!n_particles) {
        log_exit("zero_centroid");
        return;
        }
    for (i=sum=0; i<n_particles; i++)
        sum += particle[i][coord];
#if !SDDS_MPI_IO
    sum /= n_particles;
#else
    if (isSlave) {
      double total_sum; 
      long total_particles;
      MPI_Allreduce (&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, workers);
      MPI_Allreduce (&n_particles, &total_particles, 1, MPI_LONG, MPI_SUM, workers);
      sum = total_sum/total_particles;
    }
#endif
    for (i=0; i<n_particles; i++)
        particle[i][coord] -= sum;
    log_exit("zero_centroid");
    }

long dynap_distribution(double **particle, long n_particles, double sx, double sy,
        long nx, long ny)
{
    long ix, iy, ip;
    double x, y, dx=0.0, dy=0.0;
    double *coord;

    log_entry("dynap_distribution");

    if ((nx<1 || ny<1) || (nx==1 && ny==1 && n_particles!=1)) {
        if ((nx = sqrt(n_particles*1.0))<=1) 
            bomb("too few particles requested for dynamic-aperture beam type", NULL);
        }
    if ((ny = n_particles/nx)<1)
        bomb("invalid particle parameters for dynamic-aperture beam type--increase number of particles", NULL);

    if (nx>1)
        dx = 2*sx/(nx-1);
    if (ny>1)
        dy = sy/(ny-1);
    x = -sx;
    for (ix=ip=0; ix<nx; ix++) {
        y = 0;
        for (iy=0; iy<ny; iy++, ip++) {
            coord = particle[ip];
            coord[0] = x;
            coord[2] = y;
            coord[1] = coord[3] = coord[4] = coord[5] = 0;
            y += dy;
            }
        x += dx;
        }
    log_exit("dynap_distribution");
    return(nx*ny);
    }

#if SDDS_MPI_IO
long dynap_distribution_p(double **particle, long n_particles, double sx, double sy,
        long nx, long ny)
{
    long ix, iy, ip;
    double x, y, dx=0.0, dy=0.0;
    double *coord;
    long n_particles_total, remainder, start_particle, end_particle;

    log_entry("dynap_distribution");
    MPI_Allreduce(&n_particles, &n_particles_total, 1,  MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if ((nx<1 || ny<1) || (nx==1 && ny==1 && n_particles_total!=1)) {
        if ((nx = sqrt(n_particles_total*1.0))<=1) 
            bomb("too few particles requested for dynamic-aperture beam type", NULL);
        }
    if ((ny = n_particles_total/nx)<1)
        bomb("invalid particle parameters for dynamic-aperture beam type--increase number of particles", NULL);
    if (nx>1)
        dx = 2*sx/(nx-1);
    if (ny>1)
        dy = sy/(ny-1);
    x = -sx;
/* distribute particles to different processors according to the ip */
    if (isSlave) {
       n_particles_total = nx*ny;	
       n_particles = n_particles_total/(n_processors-1);
       remainder =  n_particles_total%(n_processors-1);	
       if (myid<=remainder) {
	 n_particles++;
	 start_particle = n_particles*(myid-1);
	 end_particle = n_particles*myid-1;
       }
       else {
	start_particle = n_particles*(myid-1) + remainder;
 	end_particle = n_particles*myid + remainder-1;
      }
    for (ix=ip=0; ix<nx; ix++) {
        y = 0;
        for (iy=0; iy<ny; iy++, ip++) {
	  if(ip>=start_particle && ip<=end_particle) {
            coord = particle[ip-start_particle];
            coord[0] = x;
            coord[2] = y;
            coord[1] = coord[3] = coord[4] = coord[5] = 0;
	  }
            y += dy;
            }
        x += dx;
        }
  }
    log_exit("dynap_distribution");

    if (isMaster)
       return 0;
    else	
        return(n_particles);
}
#endif

void transform_from_normalized_coordinates(double **part, long n_part, long offset, double beta, double alpha)
{
    long i;
    double *coord, u, up;
    double c1, c2;

    log_entry("transform_from_normalized_coordinates");

    c1 = sqrt(beta);
    c2 = alpha/c1;

    for (i=0; i<n_part; i++) {
        coord = part[i];
        u  = coord[0+offset];
        up = coord[1+offset];
        coord[0+offset] = c1*u;
        coord[1+offset] = c1*up - c2*u;
        }
    log_exit("transform_from_normalized_coordinates");
    }

void couple_coordinates(double **particle, long n_particles, long offset, double angle)
{
    long i;
    double slope;

    log_entry("couple_coordinates");

    slope = tan(angle);
    for (i=0; i<n_particles; i++)
        particle[i][offset] += particle[i][offset+1]*slope;
    log_exit("couple_coordinates");
    }

void set_beam_centroids(
    double **particle, 
    long offset, 
    long n_particles, 
    double cent_posi, 
    double cent_slope
    )
{
    long i;
    for (i=0; i<n_particles; i++) {
        particle[i][offset]   += cent_posi;
        particle[i][offset+1] += cent_slope;
        }
    }

void gaussian_4d_distribution(
                              double **particle, 
                              long n_particles, 
                              long offset, 
                              double s1, 
                              double s2,
                              double s3, 
                              double s4,
                              double cutoff,
                              long halo
    )
{
    double x1, x2, x3, x4, cutoff2;
    long i_particle, flag;

    log_entry("gaussian_4d_distribution");

    cutoff2 = sqr(cutoff);
    for (i_particle=0; i_particle<n_particles; i_particle++) {
        do {
            x1 = gauss_rn_lim(0.0, 1.0, 0.0, random_4);
            x2 = gauss_rn_lim(0.0, 1.0, 0.0, random_4);
            x3 = gauss_rn_lim(0.0, 1.0, 0.0, random_4);
            x4 = gauss_rn_lim(0.0, 1.0, 0.0, random_4);
            flag = sqr(x1)+sqr(x2)+sqr(x3)+sqr(x4)<=cutoff2;
            if ((!halo && flag) || (halo && !flag))
              break;
          } while (1);
        particle[i_particle][0+offset] = x1*s1;
        particle[i_particle][1+offset] = x2*s2;
        particle[i_particle][2+offset] = x3*s3;
        particle[i_particle][3+offset] = x4*s4;
        }
    log_exit("gaussian_4d_distribution");
    }

void uniform_4d_distribution(
    double **particle, 
    long n_particles, 
    long offset, 
    double max1, 
    double max2, 
    double max3, 
    double max4, 
    double cutoff
    )
{
    long i_particle;
    double rnd1, rnd2, rnd3, rnd4;

    log_entry("uniform_4d_distribution");

    max1 *= 2*cutoff;
    max2 *= 2*cutoff;
    max3 *= 2*cutoff;
    max4 *= 2*cutoff;

    for (i_particle=0; i_particle<n_particles; i_particle++) {
        do {
            rnd1 = random_4(1)-.5;
            rnd2 = random_4(1)-.5;
            rnd3 = random_4(1)-.5;
            rnd4 = random_4(1)-.5;
            } while (rnd1*rnd1+rnd2*rnd2+rnd3*rnd3+rnd4*rnd4 > .25) ;
        particle[i_particle][0+offset] = max1*rnd1;
        particle[i_particle][1+offset] = max2*rnd2;
        particle[i_particle][2+offset] = max3*rnd3;
        particle[i_particle][3+offset] = max4*rnd4;
        }
    log_exit("uniform_4d_distribution");
    }

