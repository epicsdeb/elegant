/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: analyze.c
 * purpose: Do tracking to find transfer matrix and tunes.
 *          See file analyze.nl for input parameters.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "analyze.h"
#include "matlib.h"

#define CMATRIX_OFFSET 0
#define RMATRIX_OFFSET 6
#define X_BETA_OFFSET RMATRIX_OFFSET+36
#define Y_BETA_OFFSET X_BETA_OFFSET+2
#define DETR_OFFSET Y_BETA_OFFSET+2
#define X_ETA_OFFSET DETR_OFFSET+1
#define Y_ETA_OFFSET X_ETA_OFFSET+2
#define CLORB_ETA_OFFSET Y_ETA_OFFSET+2
#define N_ANALYSIS_COLUMNS CLORB_ETA_OFFSET+4


SDDS_DEFINITION analysis_column[N_ANALYSIS_COLUMNS] = {
    {"C1", "&column name=C1, units=m, symbol=\"C$b1$n\", type=double &end"},
    {"C2", "&column name=C2, symbol=\"C$b2$n\", type=double &end"},
    {"C3", "&column name=C3, units=m, symbol=\"C$b3$n\", type=double &end"},
    {"C4", "&column name=C4, symbol=\"C$b4$n\", type=double &end"},
    {"C5", "&column name=C5, units=m, symbol=\"C$b5$n\", type=double &end"},
    {"C6", "&column name=C6, symbol=\"C$b6$n\", type=double &end"},
    {"R11", "&column name=R11, symbol=\"R$b11$n\", type=double &end"},
    {"R12", "&column name=R12, units=m, symbol=\"R$b12$n\", type=double &end"},
    {"R13", "&column name=R13, symbol=\"R$b13$n\", type=double &end"},
    {"R14", "&column name=R14, units=m, symbol=\"R$b14$n\", type=double &end"},
    {"R15", "&column name=R15, symbol=\"R$b15$n\", type=double &end"},
    {"R16", "&column name=R16, units=m, symbol=\"R$b16$n\", type=double &end"},
    {"R21", "&column name=R21, units=1/m, symbol=\"R$b21$n\", type=double &end"},
    {"R22", "&column name=R22, symbol=\"R$b22$n\", type=double &end"},
    {"R23", "&column name=R23, units=1/m, symbol=\"R$b23$n\", type=double &end"},
    {"R24", "&column name=R24, symbol=\"R$b24$n\", type=double &end"},
    {"R25", "&column name=R25, units=1/m, symbol=\"R$b25$n\", type=double &end"},
    {"R26", "&column name=R26, symbol=\"R$b26$n\", type=double &end"},
    {"R31", "&column name=R31, symbol=\"R$b31$n\", type=double &end"},
    {"R32", "&column name=R32, units=m, symbol=\"R$b32$n\", type=double &end"},
    {"R33", "&column name=R33, symbol=\"R$b33$n\", type=double &end"},
    {"R34", "&column name=R34, units=m, symbol=\"R$b34$n\", type=double &end"},
    {"R35", "&column name=R35, symbol=\"R$b35$n\", type=double &end"},
    {"R36", "&column name=R36, units=m, symbol=\"R$b36$n\", type=double &end"},
    {"R41", "&column name=R41, units=1/m, symbol=\"R$b41$n\", type=double &end"},
    {"R42", "&column name=R42, symbol=\"R$b42$n\", type=double &end"},
    {"R43", "&column name=R43, units=1/m, symbol=\"R$b43$n\", type=double &end"},
    {"R44", "&column name=R44, symbol=\"R$b44$n\", type=double &end"},
    {"R45", "&column name=R45, units=1/m, symbol=\"R$b45$n\", type=double &end"},
    {"R46", "&column name=R46, symbol=\"R$b46$n\", type=double &end"},
    {"R51", "&column name=R51, symbol=\"R$b51$n\", type=double &end"},
    {"R52", "&column name=R52, units=m, symbol=\"R$b52$n\", type=double &end"},
    {"R53", "&column name=R53, symbol=\"R$b53$n\", type=double &end"},
    {"R54", "&column name=R54, units=m, symbol=\"R$b54$n\", type=double &end"},
    {"R55", "&column name=R55, symbol=\"R$b55$n\", type=double &end"},
    {"R56", "&column name=R56, units=m, symbol=\"R$b56$n\", type=double &end"},
    {"R61", "&column name=R61, units=1/m, symbol=\"R$b61$n\", type=double &end"},
    {"R62", "&column name=R62, symbol=\"R$b62$n\", type=double &end"},
    {"R63", "&column name=R63, units=1/m, symbol=\"R$b63$n\", type=double &end"},
    {"R64", "&column name=R64, symbol=\"R$b64$n\", type=double &end"},
    {"R65", "&column name=R65, units=1/m, symbol=\"R$b65$n\", type=double &end"},
    {"R66", "&column name=R66, symbol=\"R$b66$n\", type=double &end"},
    {"betax", "&column name=betax, units=m, symbol=\"$gb$r$bx$n\", type=double &end"},
    {"nux", "&column name=nux, symbol=\"$gn$r$bx$n\", type=double &end"},
    {"betay", "&column name=betay, units=m, symbol=\"$gb$r$by$n\", type=double &end"},
    {"nuy", "&column name=nuy, symbol=\"$gn$r$by$n\", type=double &end"},
    {"detR", "&column name=detR, symbol=\"detR\", type=double &end"},
    {"etax", "&column name=etax, units=m, symbol=\"$gc$r$bx$n\", type=double &end"},
    {"etapx", "&column name=etapx, symbol=\"$gc$r$bx$n$a'$n\", type=double &end"},
    {"etay", "&column name=etay, units=m, symbol=\"$gc$ry\", type=double &end"},
    {"etapy", "&column name=etapy, symbol=\"$gc$r$by$n$a'$n\", type=double &end"},
    {"cetax", "&column name=cetax, units=m, symbol=\"$gc$r$bxc$n\", type=double &end"},
    {"cetapx", "&column name=cetapx, symbol=\"$gc$r$bxc$n$a'$n\", type=double &end"},
    {"cetay", "&column name=cetay, units=m, symbol=\"$gc$r$byc$n\", type=double &end"},
    {"cetapy", "&column name=cetapy, symbol=\"$gc$r$byc$n$a'$n\", type=double &end"},
    } ;

#define IP_STEP 0
#define IP_SUBSTEP 1
#define N_PARAMETERS 2
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"Substep", "&parameter name=Substep, type=long, description=\"Simulation substep\" &end"},
    } ;

static long SDDS_analyze_initialized = 0;
static SDDS_TABLE SDDS_analyze;

void setup_transport_analysis(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control,
    ERRORVAL *errcon
    )
{
    log_entry("setup_transport_analysis");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&analyze_map, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &analyze_map);

    /* check for data errors */
    if (!output)
        bombElegant("no output filename specified", NULL);
    if (n_points!=2 && n_points!=4)
        bombElegant("n_points must be either 2 or 4", NULL);

    output = compose_filename(output, run->rootname);
    SDDS_ElegantOutputSetup(&SDDS_analyze, output, SDDS_BINARY, 1, "transport analysis", 
                            run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                            analysis_column, N_ANALYSIS_COLUMNS, "setup_transport_analysis", 
                            SDDS_EOS_NEWFILE);
    SDDS_analyze_initialized = 1;

    if (!SDDS_DefineSimpleColumns(&SDDS_analyze, control->n_elements_to_vary,
                                  control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumns(&SDDS_analyze, errcon->n_items, errcon->quan_name, errcon->quan_unit, 
                                  SDDS_DOUBLE)) {
        SDDS_SetError("Unable to define additional SDDS columns (setup_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (!SDDS_SaveLayout(&SDDS_analyze) || !SDDS_WriteLayout(&SDDS_analyze)) {
        SDDS_SetError("Unable to write SDDS layout for transport analysis");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    log_exit("setup_transport_analysis");
    }

void do_transport_analysis(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline,
    double *orbit
    )
{
    static double **coord;
    static double *data;
    static double *offset;
    static long initialized = 0;
    static MATRIX *R, *Rc;
    double p_central, beta, tune;
    long n_track, n_trpoint, i, j, effort, index;
    double sin_phi, cos_phi, det;
    static double *orbit_p, *orbit_m;
    TRAJECTORY *clorb=NULL;

    log_entry("do_transport_analysis");
        
    if (center_on_orbit && !orbit)
        bombElegant("you've asked to center the analysis on the closed orbit, but you didn't issue a closed_orbit command", NULL);

    if (!initialized) {
        coord = (double**)czarray_2d(sizeof(**coord), 1+6*n_points, 7);
        data  = (double*)tmalloc(sizeof(*data)*SDDS_analyze.layout.n_columns);
        offset = (double*)tmalloc(sizeof(*offset)*6);
        m_alloc(&R , 6, 6);
        m_alloc(&Rc, 6, 6);
        orbit_p = tmalloc(sizeof(*orbit_p)*6);
        orbit_m = tmalloc(sizeof(*orbit_m)*6);
        initialized = 1;
        }
    if (center_on_orbit)
        clorb = tmalloc(sizeof(*clorb)*(beamline->n_elems+1));

    n_track = n_trpoint = n_points*6+1;
    for (j=0; j<7; j++)
        for (i=0; i<n_track; i++)
            coord[i][j] = (orbit && center_on_orbit?orbit[j]:0);


    /* particles 0 and 1 are for d/dx */
    offset[0] = delta_x ;
    coord[0][0] += delta_x ;
    coord[1][0] -= delta_x ;
    /* particles 2 and 3 are for d/dxp */
    offset[1] = delta_xp;
    coord[2][1] += delta_xp;
    coord[3][1] -= delta_xp;
    /* particles 4 and 5 are for d/dy */
    offset[2] = delta_y ;
    coord[4][2] += delta_y ;
    coord[5][2] -= delta_y ;
    /* particles 6 and 7 are for d/dyp */
    offset[3] = delta_y ;
    coord[6][3] += delta_yp ;
    coord[7][3] -= delta_yp ;
    /* particles 8 and 9 are for d/ds */
    offset[4] = delta_s ;
    coord[8][4] += delta_s ;
    coord[9][4] -= delta_s ;
    /* particles 10 and 11 are for d/dp */
    offset[5] = delta_dp;
    coord[10][5] += delta_dp;
    coord[11][5] -= delta_dp;
    if (n_points==4) {
        /* particles 12 and 13 are for d/dx */
        coord[12][0] += 3*delta_x ;
        coord[13][0] -= 3*delta_x ;
        /* particles 14 and 15 are for d/dxp */
        coord[14][1] += 3*delta_xp;
        coord[15][1] -= 3*delta_xp;
        /* particles 16 and 17 are for d/dy */
        coord[16][2] += 3*delta_y ;
        coord[17][2] -= 3*delta_y ;
        /* particles 18 and 19 are for d/dyp */
        coord[18][3] += 3*delta_yp ;
        coord[19][3] -= 3*delta_yp ;
        /* particles 20 and 21 are for d/ds */
        coord[20][4] += 3*delta_s ;
        coord[21][4] -= 3*delta_s ;
        /* particles 22 and 23 are for d/dp */
        coord[22][5] += 3*delta_dp;
        coord[23][5] -= 3*delta_dp;
        }
    /* particle n_track-1 is the reference particle */

    effort = 0;
    p_central = run->p_central;
    if (do_tracking(NULL, coord, n_trpoint, &effort, beamline, &p_central, 
                    NULL, NULL, NULL, NULL, run, control->i_step, 
		    (control->fiducial_flag&
		     (FIRST_BEAM_IS_FIDUCIAL+RESTRICT_FIDUCIALIZATION))
		    +SILENT_RUNNING+TEST_PARTICLES,
		    control->n_passes, 0,
                    NULL, NULL, NULL, NULL, NULL)!=n_track) {
        fputs("warning: particle(s) lost during transport analysis--continuing with next step", stdout);
        log_exit("do_transport_analysis");
        return;
        }
    if (p_central!=run->p_central && verbosity>0)
      fprintf(stdout, "Central momentum changed from %e to %e\n", run->p_central, p_central);

    if (verbosity>0){
        if (orbit && center_on_orbit) {
            fprintf(stdout, "closed orbit: \n");
            fflush(stdout);
            for (i=0; i<6; i++)
                fprintf(stdout, "%15.8e ", orbit[i]);
            fputc('\n', stdout);
            }
        fprintf(stdout, "final coordinates of refence particle: \n");
        fflush(stdout);
        for (i=0; i<6; i++)
            fprintf(stdout, "%15.8e ", coord[n_track-1][i]);
        fputc('\n', stdout);
        fflush(stdout);
        if (verbosity>1) {
          for (i=0; i<n_track-1; i++) {
            fprintf(stdout, "Particle %ld end coordinates:\n", i);
            for (j=0; j<6; j++)
              fprintf(stdout, "%15.8e%c", coord[i][j], j==5?'\n':' ');
          }
          fflush(stdout);
        }
      }
    

    for (i=0; i<6; i++) 
        data[i+CMATRIX_OFFSET] = coord[n_track-1][i];

    for (i=0; i<6; i++) {
        /* Compute R[i][j] --> data[6*i+j] */
        /* i indexes the dependent quantity */
        for (j=0; j<6; j++) {
            /* j indexes the initial coordinate value */
            if (!offset[j]) 
                R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = (i==j?1:0);
            else {
                if (n_points==2) 
                    R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = (coord[2*j][i]-coord[2*j+1][i])/(2*offset[j]);
                else
                    R->a[i][j] = data[RMATRIX_OFFSET+6*i+j] = 
                        (27*(coord[2*j][i]-coord[2*j+1][i])-(coord[2*j+12][i]-coord[2*j+13][i]))/(48*offset[j]);
                }
            }
        }
    /* find NUx */
    if (fabs(cos_phi = (R->a[0][0]+R->a[1][1])/2)>=1) {
        fprintf(stdout, "warning: beamline unstable for x plane\n");
        fflush(stdout);
        beta = 0;
        tune = 0;
        }
    else {
        beta    = fabs(R->a[0][1]/sin(acos(cos_phi)));
        sin_phi = R->a[0][1]/beta;
        if ((tune = atan2(sin_phi, cos_phi)/PIx2)<0)
            tune += 1;
        }
    data[X_BETA_OFFSET  ] = beta;
    data[X_BETA_OFFSET+1] = tune;
    /* find NUy */
    if (fabs(cos_phi = (R->a[2][2]+R->a[3][3])/2)>=1) {
        fprintf(stdout, "warning: beamline unstable for y plane\n");
        fflush(stdout);
        beta = 0;
        tune = 0;
        }
    else {
        beta    = fabs(R->a[2][3]/sin(acos(cos_phi)));
        sin_phi = R->a[2][3]/beta;
        if ((tune = atan2(sin_phi, cos_phi)/PIx2)<0)
            tune += 1;
        }
    data[Y_BETA_OFFSET  ] = beta;
    data[Y_BETA_OFFSET+1] = tune;
    m_copy(Rc, R);
    data[DETR_OFFSET] = m_det(R);

    /* compute etax and etax' */
    if ((det = (2 - R->a[0][0] - R->a[1][1]))<=0) {
        fprintf(stdout, "error: beamline unstable for x plane--can't match dispersion functions\n");
        fflush(stdout);
        det = 1e-6;
        }
    data[X_ETA_OFFSET  ] = ((1-R->a[1][1])*R->a[0][5]+R->a[0][1]*R->a[1][5])/det;
    data[X_ETA_OFFSET+1] = (R->a[1][0]*R->a[0][5] + (1-R->a[0][0])*R->a[1][5])/det;

    /* compute etay and etay' */
    if ((det = (2 - R->a[2][2] - R->a[3][3]))<=0) {
        fprintf(stdout, "error: beamline unstable for y plane--can't match dispersion functions\n");
        fflush(stdout);
        det = 1e-6;
        }
    data[Y_ETA_OFFSET  ] = ((1-R->a[3][3])*R->a[2][5]+R->a[2][3]*R->a[3][5])/det;
    data[Y_ETA_OFFSET+1] = (R->a[3][2]*R->a[2][5] + (1-R->a[2][2])*R->a[3][5])/det;

    if (delta_dp && center_on_orbit) {
        if (!beamline->matrix)
            beamline->matrix = full_matrix(&(beamline->elem), run, 1);
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, delta_dp, 1, 0, NULL, 1.0,
                          NULL);
        for (i=0; i<6; i++)
            orbit_p[i] = clorb[0].centroid[i];
        find_closed_orbit(clorb, 1e-12, 20, beamline, beamline->matrix, run, -delta_dp, 1, 0, NULL, 1.0,
                          NULL);
        for (i=0; i<6; i++)
            orbit_m[i] = clorb[0].centroid[i];
        for (i=0; i<4; i++)
            data[CLORB_ETA_OFFSET+i] = (orbit_p[i]-orbit_m[i])/(2*delta_dp);
        }
    else
        for (i=0; i<4; i++)
            data[CLORB_ETA_OFFSET+i] = 0;

    index = N_ANALYSIS_COLUMNS;
    for (i=0; i<control->n_elements_to_vary; i++, index++)
        data[index] = control->varied_quan_value[i];
    for (i=0 ; i<errcon->n_items; i++, index++)
        data[index] = errcon->error_value[i];

    if (!SDDS_StartTable(&SDDS_analyze, 1)) {
        fprintf(stdout, "Unable to start SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
        }
    if (!SDDS_SetParameters(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                            IP_STEP, control->i_step, IP_SUBSTEP, control->i_vary, -1)) {
        SDDS_SetError("Unable to set SDDS parameter values (do_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<SDDS_ColumnCount(&SDDS_analyze); i++)
        if (!SDDS_SetRowValues(&SDDS_analyze, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0,
                               i, data[i], -1)) {
            fprintf(stdout, "Unable to set SDDS column %s (do_transport_analysis)\n", 
                    analysis_column[i].name);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exitElegant(1);
            }
    if (!SDDS_WriteTable(&SDDS_analyze)) {
        fprintf(stdout, "Unable to write SDDS table (do_transport_analysis)");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
        }
    if (!inhibitFileSync)
      SDDS_DoFSync(&SDDS_analyze);
    
    if (verbosity>0) {
        fprintf(stdout, "C : ");
	for (i=0; i<6; i++) 
  	    fprintf(stdout, "%20.13e ", coord[n_track-1][i]);
        fputc('\n', stdout);
        for (i=0; i<6; i++) {
            fprintf(stdout, "R%ld: ", i+1);
            fflush(stdout);
            for (j=0; j<6; j++) 
                fprintf(stdout, "%20.13e ", R->a[i][j]);
                fflush(stdout);
            fputc('\n', stdout);
            }
        fprintf(stdout, "horizontal:   tune = %20.13e  beta = %20.13e  eta = %20.13e  eta' = %20.13e\n",
            data[X_BETA_OFFSET+1], data[X_BETA_OFFSET], data[X_ETA_OFFSET], data[X_ETA_OFFSET+1]);
        fflush(stdout);
        fprintf(stdout, "vertical  :   tune = %20.13e  beta = %20.13e  eta = %20.13e  eta' = %20.13e\n",
            data[Y_BETA_OFFSET+1], data[Y_BETA_OFFSET], data[Y_ETA_OFFSET], data[Y_ETA_OFFSET+1]);
        fflush(stdout);
        fprintf(stdout, "determinant of R = 1 + %20.13e\n", data[DETR_OFFSET]-1);
        fflush(stdout);
        fprintf(stdout, "dispersion functions from closed-orbit calculations:\nx: %e m    %e\ny: %e m    %e\n",
            data[CLORB_ETA_OFFSET  ], data[CLORB_ETA_OFFSET+1],
            data[CLORB_ETA_OFFSET+2], data[CLORB_ETA_OFFSET+3]);
        fflush(stdout);
        }

    log_exit("do_transport_analysis");
    }

void finish_transport_analysis(
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    LINE_LIST *beamline
    )
{
    if (SDDS_IsActive(&SDDS_analyze) && !SDDS_Terminate(&SDDS_analyze)) {
        SDDS_SetError("Problem terminating SDDS output (finish_transport_analysis)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }        
    SDDS_analyze_initialized = 0;
    }


VMATRIX *determineMatrix(RUN *run, ELEMENT_LIST *eptr, double *startingCoord, double *stepSize)
{
  double **coord;
  long n_track, i, j;
  VMATRIX *M;
  double **R, *C;
  double defaultStep[6] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5};
  long ltmp1, ltmp2;
  double dgamma, dtmp1, dP[3];
  
  coord = (double**)czarray_2d(sizeof(**coord), 1+6*4, 7);

  if (stepSize==NULL)
    stepSize = defaultStep;
  
  n_track = 4*6+1;
  for (j=0; j<7; j++)
    for (i=0; i<n_track; i++)
      coord[i][j] = startingCoord ? startingCoord[j] : 0;

  /* particles 0 and 1 are for d/dx */
  coord[0][0] += stepSize[0] ;
  coord[1][0] -= stepSize[0] ;
  /* particles 2 and 3 are for d/dxp */
  coord[2][1] += stepSize[1];
  coord[3][1] -= stepSize[1];
  /* particles 4 and 5 are for d/dy */
  coord[4][2] += stepSize[2] ;
  coord[5][2] -= stepSize[2] ;
  /* particles 6 and 7 are for d/dyp */
  coord[6][3] += stepSize[3] ;
  coord[7][3] -= stepSize[3] ;
  /* particles 8 and 9 are for d/ds */
  coord[8][4] += stepSize[4] ;
  coord[9][4] -= stepSize[4] ;
  /* particles 10 and 11 are for d/delta */
  coord[10][5] += stepSize[5];
  coord[11][5] -= stepSize[5];

  /* particles 12 and 13 are for d/dx */
  coord[12][0] += 3*stepSize[0] ;
  coord[13][0] -= 3*stepSize[0] ;
  /* particles 14 and 15 are for d/dxp */
  coord[14][1] += 3*stepSize[1];
  coord[15][1] -= 3*stepSize[1];
  /* particles 16 and 17 are for d/dy */
  coord[16][2] += 3*stepSize[2] ;
  coord[17][2] -= 3*stepSize[2] ;
  /* particles 18 and 19 are for d/dyp */
  coord[18][3] += 3*stepSize[3] ;
  coord[19][3] -= 3*stepSize[3] ;
  /* particles 20 and 21 are for d/ds */
  coord[20][4] += 3*stepSize[4] ;
  coord[21][4] -= 3*stepSize[4] ;
  /* particles 22 and 23 are for d/delta */
  coord[22][5] += 3*stepSize[5];
  coord[23][5] -= 3*stepSize[5];
  /* particle n_track-1 is the reference particle (coordinates set above) */

  switch (eptr->type) {
  case T_CWIGGLER:
    ltmp1 = ((CWIGGLER*)eptr->p_elem)->isr;
    ltmp2 = ((CWIGGLER*)eptr->p_elem)->sr;
    ((CWIGGLER*)eptr->p_elem)->isr = ((CWIGGLER*)eptr->p_elem)->sr = 0;
    GWigSymplecticPass(coord, n_track, run->p_central, (CWIGGLER*)eptr->p_elem);
    ((CWIGGLER*)eptr->p_elem)->isr = ltmp1;
    ((CWIGGLER*)eptr->p_elem)->sr = ltmp2;
    break;
  case T_UKICKMAP:
    ltmp1 = ((UKICKMAP*)eptr->p_elem)->isr;
    ltmp2 = ((UKICKMAP*)eptr->p_elem)->synchRad;
    ((UKICKMAP*)eptr->p_elem)->isr = ((UKICKMAP*)eptr->p_elem)->synchRad = 0;
    if (trackUndulatorKickMap(coord, NULL, n_track, run->p_central, (UKICKMAP*)eptr->p_elem, 0)!=n_track) {
      printf("*** Error: particles lost in determineMatrix call for UKICKMAP\n");
      exitElegant(1);
    }
    ((UKICKMAP*)eptr->p_elem)->isr = ltmp1;
    ((UKICKMAP*)eptr->p_elem)->synchRad = ltmp2;
    break;
  case T_SCRIPT:
    transformBeamWithScript((SCRIPT*)eptr->p_elem, run->p_central, NULL, NULL, coord, n_track, 0,
                            NULL, 0, 2);
    break;
  case T_TWLA:
    motion(coord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    break;
  case T_RFDF:
    /* Don't actually use this */
    track_through_rf_deflector(coord, (RFDF*)eptr->p_elem,
			       coord, n_track, run->p_central, 0, eptr->end_pos, 0);
    break;
  case T_LSRMDLTR:
    ltmp1 = ((LSRMDLTR*)eptr->p_elem)->isr;
    ltmp2 = ((LSRMDLTR*)eptr->p_elem)->synchRad;
    dtmp1 = ((LSRMDLTR*)eptr->p_elem)->laserPeakPower;
    ((LSRMDLTR*)eptr->p_elem)->isr = ((LSRMDLTR*)eptr->p_elem)->synchRad = 0;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = 0;
    motion(coord, n_track, eptr->p_elem, eptr->type, &run->p_central, &dgamma, dP, NULL, 0.0);
    ((LSRMDLTR*)eptr->p_elem)->isr = ltmp1;
    ((LSRMDLTR*)eptr->p_elem)->synchRad = ltmp2;
    ((LSRMDLTR*)eptr->p_elem)->laserPeakPower = dtmp1;
    break;
  case T_FTABLE:
    field_table_tracking(coord, n_track, (FTABLE*)eptr->p_elem, run->p_central, run);
    break;
  default:
    printf("*** Error: determineMatrix called for element that is not supported!\n");
    printf("***        Seek professional help!\n");
    exitElegant(1);
    break;
  }
  
  M = tmalloc(sizeof(*M));
  M->order = 1;
  initialize_matrices(M, M->order);
  R = M->R;
  C = M->C;
  
  for (i=0; i<6; i++) {
    /* i indexes the dependent quantity */

    /* Determine C[i] */
    C[i] = coord[n_track-1][i];

    /* Compute R[i][j] */
    for (j=0; j<6; j++) {
      /* j indexes the initial coordinate value */
      if (n_points==2) 
        R[i][j] = (coord[2*j][i]-coord[2*j+1][i])/(2*stepSize[j]);
      else
        R[i][j] = 
          (27*(coord[2*j][i]-coord[2*j+1][i])-(coord[2*j+12][i]-coord[2*j+13][i]))/(48*stepSize[j]);
    }
  }

  free_czarray_2d((void**)coord, 1+4*6, 7);

  /*
  if (strlen(eptr->name)<900)
    sprintf(s, "\nElement %s#%ld matrix determined from tracking:\n", eptr->name, eptr->occurence);
  else
    sprintf(s, "\nElement matrix determined from tracking:\n");
  print_matrices(stdout, s, M);
  */

  return M;
}

void determineRadiationMatrix(VMATRIX *Mr, RUN *run, ELEMENT_LIST *eptr, double *startingCoord, double *Dr, long nSlices, long order)
{
  CSBEND csbend; CSRCSBEND *csrcsbend; BEND *sbend;
  KQUAD kquad;  QUAD *quad;
  KSEXT ksext; SEXT *sext;
  HCOR hcor; VCOR vcor; HVCOR hvcor;
  double length;
  long i, j, k, slice;
  double *accumD1, *accumD2, *dtmp;
  VMATRIX *M1, *M2, *Ml1, *Mtmp;
  ELEMENT_LIST elem;
  MATRIX *Ms;
  
  /* Accumulated diffusion matrix */
  accumD1 = tmalloc(21*sizeof(*(accumD1)));
  memset(accumD1, 0, 21*sizeof(*(accumD1)));
  accumD2 = tmalloc(21*sizeof(*(accumD2)));
  memset(accumD2, 0, 21*sizeof(*(accumD2)));

  /* Matrices for C, R matrix propagation: */
  initialize_matrices(M1=tmalloc(sizeof(*M1)), order);
  initialize_matrices(M2=tmalloc(sizeof(*M2)), order);
  /* Temporary variable for linear matrix with radiation: */
  initialize_matrices(Ml1=tmalloc(sizeof(*Ml1)), 1);
  for (i=0; i<6; i++) {
    M1->R[i][i] = 1;
    M1->C[i] = startingCoord[i];
  }
  /* Matrix for sigma matrix (D matrix) propagation */
  m_alloc(&Ms, 21, 21);

  elem.end_pos = eptr->end_pos;
  for (slice=0; slice<nSlices; slice++) {
    switch (eptr->type) {
    case T_CSBEND:
      memcpy(&csbend, (CSBEND*)eptr->p_elem, sizeof(CSBEND));
      csbend.isr = 0;
      csbend.angle /= nSlices;
      length = csbend.length /= nSlices;
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1) 
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      break;
    case T_SBEN:
      if (slice!=0)
        csbend.edge1_effects = 0;
      if (slice!=nSlices-1)
        csbend.edge2_effects = 0;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      sbend = (BEND*)eptr->p_elem;
      memset(&csbend, 0, sizeof(csbend));
      csbend.isr = 0;
      csbend.synch_rad = 1;
      length = csbend.length = sbend->length/nSlices;
      csbend.angle = sbend->angle/nSlices;
      csbend.k1 = sbend->k1;
      csbend.e1 = sbend->e1;
      csbend.e2 = sbend->e2;
      csbend.k2 = sbend->k2;
      csbend.h1 = sbend->h1;
      csbend.h2 = sbend->h2;
      csbend.hgap = sbend->hgap;
      csbend.fint = sbend->fint;
      csbend.dx = sbend->dx;
      csbend.dy = sbend->dy;
      csbend.dz = sbend->dz;
      csbend.fse = sbend->fse;
      csbend.tilt = sbend->tilt;
      csbend.etilt = sbend->etilt;
      csbend.edge1_effects = sbend->edge1_effects;
      csbend.edge2_effects = sbend->edge2_effects;
      csbend.edge_order = sbend->edge_order;
      csbend.edgeFlags = sbend->edgeFlags;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1)
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      csbend.k1 = sbend->k1;
      csbend.k2 = sbend->k2;
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      csbend.integration_order = 4;
      break;
    case T_CSRCSBEND:
      if (slice!=0)
        csbend.edge1_effects = 0;
      if (slice!=nSlices-1)
        csbend.edge2_effects = 0;
      elem.type = T_CSBEND;
      elem.p_elem = (void*)&csbend;
      csrcsbend = (CSRCSBEND*)eptr->p_elem;
      memset(&csbend, 0, sizeof(csbend));
      csbend.isr = 0;
      csbend.synch_rad = 1;
      length = csbend.length = csrcsbend->length/nSlices;
      csbend.angle = csrcsbend->angle/nSlices;
      csbend.k1 = csrcsbend->k1;
      csbend.e1 = csrcsbend->e1;
      csbend.e2 = csrcsbend->e2;
      csbend.k2 = csrcsbend->k2;
      csbend.h1 = csrcsbend->h1;
      csbend.h2 = csrcsbend->h2;
      csbend.hgap = csrcsbend->hgap;
      csbend.fint = csrcsbend->fint;
      csbend.dx = csrcsbend->dx;
      csbend.dy = csrcsbend->dy;
      csbend.dz = csrcsbend->dz;
      csbend.fse = csrcsbend->fse;
      csbend.tilt = csrcsbend->tilt;
      csbend.etilt = csrcsbend->etilt;
      csbend.edge1_effects = csrcsbend->edge1_effects;
      csbend.edge2_effects = csrcsbend->edge2_effects;
      csbend.edge_order = csrcsbend->edge_order;
      csbend.edgeFlags = csrcsbend->edgeFlags;
      if (slice!=0)
        csbend.edgeFlags &= ~BEND_EDGE1_EFFECTS;
      if (slice!=nSlices-1)
        csbend.edgeFlags &= ~BEND_EDGE2_EFFECTS;
      csbend.k1 = csrcsbend->k1;
      csbend.k2 = csrcsbend->k2;
      csbend.n_kicks = fabs(csbend.angle/0.005) + 1;
      csbend.integration_order = 4;
      break;
    case T_KQUAD:
      memcpy(&kquad, (KQUAD*)eptr->p_elem, sizeof(KQUAD));
      kquad.isr = 0;
      length = (kquad.length /= nSlices);
      kquad.n_kicks = 4 + (long)(fabs(kquad.k1)*sqr(kquad.length));
      elem.type = T_KQUAD;
      elem.p_elem = (void*)&kquad;
      break;
    case T_QUAD:
      quad = (QUAD*)eptr->p_elem;
      memset(&kquad, 0, sizeof(KQUAD));
      kquad.isr = 0;
      kquad.synch_rad = 1;
      length = (kquad.length = quad->length/nSlices);
      kquad.k1 = quad->k1;
      kquad.tilt = quad->tilt;
      if (quad->ffringe)
        bombElegant("Can't perform radiation matrix calculations when QUAD has nonzero FFRINGE parameter", NULL);
      kquad.dx = quad->dx;
      kquad.dy = quad->dy;
      kquad.dz = quad->dz;
      kquad.fse = quad->fse;
      kquad.xkick = quad->xkick;
      kquad.ykick = quad->ykick;
      kquad.xKickCalibration = quad->xKickCalibration;
      kquad.yKickCalibration = quad->yKickCalibration;
      kquad.n_kicks = 4 + (long)(fabs(kquad.k1)*sqr(kquad.length));
      kquad.integration_order = 4;
      elem.type = T_KQUAD;
      elem.p_elem = (void*)&kquad;
      break;
    case T_KSEXT:
      memcpy(&ksext, (KSEXT*)eptr->p_elem, sizeof(KSEXT));
      ksext.isr = 0;
      ksext.n_kicks = 4;
      length = (ksext.length /= nSlices);
      elem.type = T_KSEXT;
      elem.p_elem = (void*)&ksext;
      break;
    case T_SEXT:
      sext = (SEXT*)eptr->p_elem;
      memset(&ksext, 0, sizeof(KSEXT));
      ksext.isr = 0;
      ksext.synch_rad = 1;
      length = (ksext.length = sext->length/nSlices);
      ksext.k2 = sext->k2;
      ksext.tilt = sext->tilt;
      ksext.dx = sext->dx;
      ksext.dy = sext->dy;
      ksext.dz = sext->dz;
      ksext.fse = sext->fse;
      ksext.n_kicks = 4;
      ksext.integration_order = 4;
      elem.type = T_KSEXT;
      elem.p_elem = (void*)&ksext;
      break;
    case T_RFCA:
      nSlices = 1;
      elem.type = T_RFCA;
      elem.p_elem = eptr->p_elem;
      length = ((RFCA*)eptr->p_elem)->length;
      break;
    case T_TWLA:
      nSlices = 1;
      elem.type = T_TWLA;
      elem.p_elem = eptr->p_elem;
      length = ((TW_LINAC*)eptr->p_elem)->length/nSlices;
      break;
    case T_HCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.p_elem = &hcor;
      memcpy(&hcor, eptr->p_elem, sizeof(hcor));
      length = (hcor.length /= nSlices);
      hcor.lEffRad /= nSlices;
      hcor.kick /= nSlices;
      hcor.isr = 0;
      break;
    case T_VCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.p_elem = &vcor;
      memcpy(&vcor, eptr->p_elem, sizeof(vcor));
      length = (vcor.length /= nSlices);
      vcor.lEffRad /= nSlices;
      vcor.kick /= nSlices;
      vcor.isr = 0;
      break;
    case T_HVCOR:
      memcpy(&elem, eptr, sizeof(elem));
      elem.p_elem = &hvcor;
      memcpy(&hvcor, eptr->p_elem, sizeof(hvcor));
      length = (hvcor.length /= nSlices);
      hvcor.lEffRad /= nSlices;
      hvcor.xkick /= nSlices;
      hvcor.ykick /= nSlices;
      hvcor.isr = 0;
      break;
    default:
      printf("*** Error: determineRadiationMatrix called for element (%s) that is not supported!\n", eptr->name);
      printf("***        Seek professional help!\n");
      exit(1);
      break;
    }

    /* Step 1: determine effective R matrix for this element, as well as the diffusion matrix */
    determineRadiationMatrix1(Ml1, run, &elem, M1->C, accumD2); 
    /*    print_matrices(stdout, "matrix1:", Ml1); */

    /* Step 2: Propagate the diffusion matrix */
    fillSigmaPropagationMatrix(Ms->a, Ml1->R);
    for (i=0; i<21; i++) 
      for (j=0; j<21; j++)
        accumD2[i] += Ms->a[i][j]*accumD1[j];
    dtmp    = accumD1;
    accumD1 = accumD2;
    accumD2 = dtmp;
    
    /* Step 3: Copy the propagated C vector */
    memcpy(M2->C, Ml1->C, 6*sizeof(*(M2->C)));

    /* Step 4: Multiply the R matrices */
    for (i=0; i<6; i++)
      for (j=0; j<6; j++) {
        M2->R[i][j] = 0;
        for (k=0; k<6; k++) 
          M2->R[i][j] += Ml1->R[i][k]*M1->R[k][j];
      }

    Mtmp = M2;
    M2   = M1;
    M1   = Mtmp;

    elem.end_pos += length;
  }
  
  /* Copy the matrix to the caller's buffer */
  for (i=0; i<6; i++) {
    Mr->C[i] = M1->C[i];
    for (j=0; j<6; j++) {
      Mr->R[i][j] = M1->R[i][j];
    }
  }
  for (i=0; i<21; i++)
    Dr[i] = accumD1[i];

  free(accumD1);
  free(accumD2);
  free_matrices(M2); tfree(M2);
  free_matrices(M1); tfree(M1);
  free_matrices(Ml1); tfree(Ml1);
  m_free(&Ms);
}


void determineRadiationMatrix1(VMATRIX *Mr, RUN *run, ELEMENT_LIST *elem, double *startingCoord, double *D)
{
  CSBEND *csbend;
  KQUAD *kquad;
  KSEXT *ksext;
  double **coord, pCentral;
  long n_track, i, j;
  double **R, *C, Cs0;
  double stepSize[6] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-3, 1e-5};
  double sigmaDelta2;
  double dP[3], dgamma;
  VMATRIX *matrix;

  coord = (double**)czarray_2d(sizeof(**coord), 1+6*2, 7);

  Cs0 = startingCoord ? startingCoord[4]: 0;
  n_track = 2*6+1;
  for (j=0; j<6; j++)
    for (i=0; i<n_track; i++)
      coord[i][j] = (startingCoord && j!=4) ? startingCoord[j] : 0;

  /* particles 0 and 1 are for d/dx, 2 and 3 are for d/dxp, etc */
  /* particle n_track-1 is the reference particle (coordinates set above) */
  for (i=0; i<6; i++) {
    coord[2*i+0][i] += stepSize[i] ;
    coord[2*i+1][i] -= stepSize[i] ;
  }
  
  if (D) {
    for (i=0; i<21; i++)
        D[i] = 0;
  }
  sigmaDelta2 = 0;
  
  switch (elem->type) {
  case T_CSBEND:
    csbend = (CSBEND*)elem->p_elem;
    track_through_csbend(coord, n_track, csbend, 0, run->p_central, NULL, elem->end_pos-csbend->length, &sigmaDelta2);
    break;
  case T_SBEN:
    track_particles(coord, elem->matrix, coord, n_track);
    break;
  case T_KQUAD:
    kquad = (KQUAD*)elem->p_elem;
    multipole_tracking2(coord, n_track, elem, 0.0, run->p_central, NULL, elem->end_pos-kquad->length,
                       0.0, 0.0, 0, NULL, &sigmaDelta2);
    break;
  case T_KSEXT:
    ksext = (KSEXT*)elem->p_elem;
    multipole_tracking2(coord, n_track, elem, 0.0, run->p_central, NULL, elem->end_pos-ksext->length,
                       0.0, 0.0, 0, NULL, &sigmaDelta2);
    break;
  case T_RFCA:
    pCentral = run->p_central;
    for (i=0; i<n_track; i++)
      coord[i][4] += Cs0;
    Cs0 = 0;
    simple_rf_cavity(coord, n_track, (RFCA*)elem->p_elem, NULL, &pCentral, elem->end_pos);
    break;
  case T_HCOR:
  case T_VCOR:
  case T_HVCOR:
    matrix = compute_matrix(elem, run, NULL);
    track_particles(coord, matrix, coord, n_track);
    addCorrectorRadiationKick(coord, n_track, elem, elem->type, run->p_central, &sigmaDelta2, 1);
    free_matrices(matrix); free(matrix);
    break;
  case T_TWLA:
    pCentral = run->p_central;
    motion(coord, n_track, elem->p_elem, elem->type, &pCentral, &dgamma, dP, NULL, 0.0);
    break;
  default:
    printf("*** Error: determineRadiationMatrix1 called for element (%s) that is not supported!\n", elem->name);
    printf("***        Seek professional help!\n");
    exitElegant(1);
    break;
  }
  
  R = Mr->R;
  C = Mr->C;
  for (i=0; i<6; i++) {
    /* i indexes the dependent quantity */
    
    /* Determine C[i] */
    C[i] = coord[n_track-1][i];
    
    /* Compute R[i][j] */
    for (j=0; j<6; j++) {
      /* j indexes the initial coordinate value */
      R[i][j] = (coord[2*j][i]-coord[2*j+1][i])/(2*stepSize[j]);
    }
  }

  D[sigmaIndex3[5][5]] = sigmaDelta2;
  
  C[4] += Cs0;
  free_czarray_2d((void**)coord, 1+2*6, 7);

}

