/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: floor.c
 * purpose: computation/output of floor coordinates
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "matlib.h"

static SDDS_TABLE SDDS_floor;

#define IC_S 0
#define IC_X 1
#define IC_Y 2
#define IC_Z 3
#define IC_THETA 4
#define IC_PHI 5
#define IC_PSI 6
#define IC_ELEMENT 7
#define IC_OCCURENCE 8
#define IC_TYPE 9
#define N_COLUMNS 10
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"X", "&column name=X, units=m, type=double, description=\"Transverse survey coordinate\" &end"},
    {"Y", "&column name=Y, units=m, type=double, description=\"Transverse survey coordinate\" &end"},
    {"Z", "&column name=Z, units=m, type=double, description=\"Longitudinal survey coordinate\" &end"},
    {"theta", "&column name=theta, symbol=\"$gq$r\", units=radians, type=double, description=\"Survey angle\" &end"},
    {"phi", "&column name=phi, symbol=\"$gf$r\", units=radians, type=double, description=\"Survey angle\" &end"},
    {"psi", "&column name=psi, symbol=\"$gy$r\", units=radians, type=double, description=\"Roll angle\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
} ;

#include "floor.h"

long advanceFloorCoordinates(MATRIX *V1, MATRIX *W1, MATRIX *V0, MATRIX *W0,
                             double *theta, double *phi, double *psi, double *s,
                             ELEMENT_LIST *elem, ELEMENT_LIST *last_elem, 
                             SDDS_DATASET *SDDS_floor, long row_index);
void computeSurveyAngles(double *theta, double *phi, double *psi, MATRIX *W);
double nearbyAngle(double angle, double reference);

void setupSurveyAngleMatrix(MATRIX *W0, double theta0, double phi0, double psi0) 
{
  MATRIX *temp, *Theta, *Phi, *Psi;

  m_alloc(&Theta, 3, 3);
  m_alloc(&Phi, 3, 3);
  m_alloc(&Psi, 3, 3);
  m_zero(Theta);
  m_zero(Phi);
  m_zero(Psi);
  m_alloc(&temp, 3, 3);

  Theta->a[0][0] = Theta->a[2][2] = cos(theta0);
  Theta->a[0][2] = -(Theta->a[2][0] = -sin(theta0));
  Theta->a[1][1] = 1;

  Phi->a[1][1] = Phi->a[2][2] = cos(phi0);
  Phi->a[1][2] = -(Phi->a[2][1] = -sin(phi0));
  Phi->a[0][0] = 1;

  Psi->a[0][0] = Psi->a[1][1] = cos(psi0);
  Psi->a[0][1] = -(Psi->a[1][0] = sin(psi0));
  Psi->a[2][2] = 1;
  
  m_mult(temp, Theta, Phi);
  m_mult(W0, temp, Psi);

  m_free(&Theta);
  m_free(&Phi);
  m_free(&Psi);
  m_free(&temp);
}

void output_floor_coordinates(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
  ELEMENT_LIST *elem, *last_elem;
  long n_points, row_index;
  MATRIX *V0, *V1;
  MATRIX *W0, *W1;
  double theta, phi, psi, s;
  
  log_entry("output_floor_coordinates");

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&floor_coordinates, nltext);
  if (echoNamelists) print_namelist(stdout, &floor_coordinates);
  
  if (magnet_centers && vertices_only)
    bomb("you can simultaneously request magnet centers and vertices only output", NULL);
  if (filename)
    filename = compose_filename(filename, run->rootname);
  else
    bomb("filename must be given for floor coordinates", NULL);
  
  SDDS_ElegantOutputSetup(&SDDS_floor, filename, SDDS_BINARY, 1, "floor coordinates", 
                          run->runfile, run->lattice, NULL, 0,
                          column_definition, N_COLUMNS, "floor coordinates", 
                          SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

  n_points = beamline->n_elems+1;
  if (vertices_only)
    n_points = 2;
  last_elem = NULL;
  if (include_vertices || vertices_only) {
    elem = &(beamline->elem);
    while (elem) {
      switch (elem->type) {
      case T_RBEN: 
      case T_SBEN: 
      case T_KSBEND: 
      case T_CSBEND:
      case T_CSRCSBEND:
        n_points ++;
        break;
      default:
        break;
      }
      last_elem = elem;
      elem = elem->succ;
    }
  }

  if (!SDDS_StartTable(&SDDS_floor, 2*n_points)) {
    SDDS_SetError("Unable to start SDDS table (output_floor_coordinates)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  row_index = 0;
  if (!SDDS_SetRowValues(&SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
                         IC_S, (double)0.0, IC_X, X0, IC_Y, Y0, IC_Z, Z0, 
                         IC_THETA, theta0, IC_PHI, phi0, IC_PSI, psi0,
                         IC_ELEMENT, "_BEG_", IC_OCCURENCE, (long)1, IC_TYPE, "MARK", -1)) {
    SDDS_SetError("Unable to set SDDS row (output_floor_coordinates)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  elem = &(beamline->elem);

  m_alloc(&V0, 3, 1);
  m_alloc(&V1, 3, 1);
  V0->a[0][0] = X0;
  V0->a[1][0] = Y0;
  V0->a[2][0] = Z0;

  m_alloc(&W0, 3, 3);
  m_alloc(&W1, 3, 3);
  setupSurveyAngleMatrix(W0, theta=theta0, phi=phi0, psi=psi0);

  s = 0;
  while (elem) {
    if (elem->type==T_STRAY) {
      STRAY *stray;
      stray = elem->p_elem;
      if (!stray->WiInitialized) {
        m_alloc((MATRIX**)(&(stray->Wi)), 3, 3);
        stray->WiInitialized = 1;
      }
      m_invert((MATRIX*)stray->Wi, W0);
    }
    row_index = advanceFloorCoordinates(V1, W1, V0, W0, &theta, &phi, &psi, &s, 
                                        elem, last_elem, &SDDS_floor, row_index);
    m_copy(W0, W1);
    m_copy(V0, V1);
    elem = elem->succ;
  }
  if (!SDDS_WriteTable(&SDDS_floor)) {
    SDDS_SetError("Unable to write floor coordinate data (output_floor_coordinates)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_floor);
  if (!SDDS_Terminate(&SDDS_floor)) {
    SDDS_SetError("Unable to terminate SDDS file (output_floor_coordinates)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  m_free(&V0);
  m_free(&V1);
  m_free(&W0);
  m_free(&W1);
}

long advanceFloorCoordinates(MATRIX *V1, MATRIX *W1, MATRIX *V0, MATRIX *W0,
                             double *theta, double *phi, double *psi, double *s,
                             ELEMENT_LIST *elem, ELEMENT_LIST *last_elem, 
                             SDDS_DATASET *SDDS_floor, long row_index)
{
  double dX, dY, dZ, rho=0.0, angle, coord[3], sangle[3], length;
  long is_bend, is_misalignment, is_magnet, is_rotation, i, is_alpha, is_mirror;
  BEND *bend; KSBEND *ksbend; CSBEND *csbend; MALIGN *malign; CSRCSBEND *csrbend;
  ROTATE *rotate; ALPH *alpha;
  char label[200];
  static MATRIX *temp33, *tempV, *R, *S, *T, *TInv;
  static long matricesAllocated = 0;
  double theta0, phi0, psi0, tilt;
  
  if (!matricesAllocated) {
    m_alloc(&temp33, 3, 3);
    m_alloc(&tempV, 3, 1);
    m_alloc(&R, 3, 1);
    m_alloc(&S, 3, 3);
    m_alloc(&T, 3, 3);
    m_alloc(&TInv, 3, 3);
    matricesAllocated = 1;
  }
  
  is_bend = is_magnet = is_rotation = is_misalignment = is_alpha = is_mirror = 0;
  length = dX = dY = dZ = tilt = angle = 0;
  switch (elem->type) {
  case T_RBEN: case T_SBEN:
    is_bend = 1;
    bend = (BEND*)elem->p_elem;
    angle = bend->angle;
    rho = (length=bend->length)/angle;
    tilt = bend->tilt;
    break;
  case T_KSBEND:
    is_bend = 1;
    ksbend = (KSBEND*)elem->p_elem;
    angle = ksbend->angle;
    rho = (length=ksbend->length)/angle;
    tilt = ksbend->tilt;
    break;
  case T_CSBEND:
    is_bend = 1;
    csbend = (CSBEND*)elem->p_elem;
    angle = csbend->angle;
    rho = (length=csbend->length)/angle;
    tilt = csbend->tilt;
    break;
  case T_CSRCSBEND:
    is_bend = 1;
    csrbend = (CSRCSBEND*)elem->p_elem;
    angle = csrbend->angle;
    rho = (length=csrbend->length)/angle;
    tilt = csrbend->tilt;
    break;
  case T_MALIGN:
    malign = (MALIGN*)elem->p_elem;
    dX = malign->dx;
    dY = malign->dy;
    dZ = malign->dz;
    angle = atan(sqrt(sqr(malign->dxp)+sqr(malign->dyp)));
    tilt = atan2(malign->dyp, -malign->dxp);
    is_misalignment = 1;
    break;
  case T_ALPH:
    is_alpha = 1;
    alpha = (ALPH*)elem->p_elem;
    tilt = alpha->tilt;
    switch (alpha->part) {
    case 1:
      dX = alpha->xmax*sin(ALPHA_ANGLE);
      dZ = alpha->xmax*cos(ALPHA_ANGLE);
      angle = -(ALPHA_ANGLE + PI/2);
      break;
    case 2:
      dX = alpha->xmax;
      angle = -(ALPHA_ANGLE + PI/2);
      break;
    default:
      angle = -(2*ALPHA_ANGLE + PI);
      break;
    }
    break;
  case T_ROTATE:
    rotate = (ROTATE*)elem->p_elem;
    tilt = rotate->tilt;
    is_rotation = 1;
    break;
  case T_LMIRROR:
    angle = PI-2*(((LMIRROR*)elem->p_elem)->theta);
    tilt = (((LMIRROR*)elem->p_elem)->tilt);
    length = 0;
    is_mirror = 1;
    break;
  default:
    if (entity_description[elem->type].flags&HAS_LENGTH)
      length = dZ = *((double*)(elem->p_elem));
    if (entity_description[elem->type].flags&IS_MAGNET)
      is_magnet = 1;
    break;
  }
  if (elem->type!=T_FLOORELEMENT) {
    theta0 = *theta;
    phi0 = *phi;
    psi0 = *psi;
    m_identity(S);
    if (is_bend || is_mirror) {
      if (!is_mirror && SDDS_floor && (include_vertices || vertices_only)) {
	/* vertex point is reached by drifting by distance rho*tan(angle/2) */
	R->a[0][0] = R->a[1][0] = 0;
	if (angle && rho)
	  R->a[2][0] = fabs(rho*tan(angle/2));
	else
	  R->a[2][0] = length/2;
	m_mult(tempV, W0, R);
	m_add(V1, tempV, V0);
	sprintf(label, "%s-VP", elem->name);
	if (!SDDS_SetRowValues(SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
			       IC_S, *s+R->a[2][0],
			       IC_X, V1->a[0][0], IC_Y, V1->a[1][0], IC_Z, V1->a[2][0], 
			       IC_THETA, *theta, IC_PHI, *phi, IC_PSI, *psi,
			       IC_ELEMENT, label, IC_OCCURENCE, elem->occurence, IC_TYPE, "VERTEX-POINT", -1)) {
	  SDDS_SetError("Unable to set SDDS row (output_floor_coordinates.1)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      }
      if (angle && !isnan(rho)) {
        if (!is_mirror) {
          R->a[0][0] = rho*(cos(angle)-1);
          R->a[1][0] = 0;
          R->a[2][0] = rho*sin(angle);
        } else {
          R->a[0][0] = R->a[1][0] = R->a[2][0] = 0;
        }
	S->a[0][0] = S->a[2][2] = cos(angle);
	S->a[2][0] = -(S->a[0][2] = -sin(angle));
      } else {
	R->a[0][0] = R->a[1][0] = 0;
	R->a[2][0] = length;
      }
    } else if (is_misalignment || is_alpha) {
      R->a[0][0] = dX;
      R->a[1][0] = dY;
      R->a[2][0] = dZ;
      S->a[0][0] = S->a[2][2] = cos(angle);
      S->a[2][0] = -(S->a[0][2] = -sin(angle));
    } else {
      R->a[0][0] = R->a[1][0] = 0;
      R->a[2][0] = length;
    }

    if (tilt) {
      m_identity(T);
      T->a[0][0] = T->a[1][1] = cos(tilt);
      T->a[1][0] = -(T->a[0][1] = -sin(tilt));
      if (!m_mult(tempV, T, R) ||
	  !m_copy(R, tempV) ||
	  !m_mult(temp33, T, S) ||
	  !m_invert(TInv, T) ||
	  !m_mult(S, temp33, TInv))
	m_error("making tilt transformation");
    }

    m_mult(tempV, W0, R);
    m_add(V1, tempV, V0);
    m_mult(W1, W0, S);
    computeSurveyAngles(theta, phi, psi, W1);
    if ((is_bend || is_magnet) && magnet_centers) {
      for (i=0; i<3; i++)
	coord[i] = V0->a[i][0] + (V1->a[i][0] - V0->a[i][0])/2;
      sangle[0] = (*theta+theta0)/2;
      sangle[1] = (*phi+phi0)/2;
      sangle[2] = (*psi+psi0)/2;
      *s += length/2;
      sprintf(label, "%s-C", elem->name);
    }
    else {
      for (i=0; i<3; i++)
	coord[i] = V1->a[i][0];
      sangle[0] = *theta;
      sangle[1] = *phi;
      sangle[2] = *psi;
      *s += length;
      length = 0;
      strcpy_ss(label, elem->name);
    }
  } else {
    /* floor coordinate reset */
    long i;
    FLOORELEMENT *fep;
    fep = (FLOORELEMENT*)(elem->p_elem);
    for (i=0; i<3; i++) {
      V1->a[i][0] = coord[i] = fep->position[i];
      sangle[i] = fep->angle[i];
    }
    strcpy_ss(label, elem->name);
    setupSurveyAngleMatrix(W1, fep->angle[0], fep->angle[1], fep->angle[2]);
  }
  if (SDDS_floor &&
      (!vertices_only || (!last_elem || elem==last_elem))) {
    if (!SDDS_SetRowValues(SDDS_floor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_index++,
                           IC_S, *s, IC_X, coord[0], IC_Y, coord[1], IC_Z, coord[2],
                           IC_THETA, sangle[0], IC_PHI, sangle[1], IC_PSI, sangle[2],
                           IC_ELEMENT, label, IC_OCCURENCE, elem->occurence, IC_TYPE, 
                           entity_name[elem->type], -1)) {
      SDDS_SetError("Unable to set SDDS row (output_floor_coordinates.2)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  *s += length/2;
  return row_index;
}


void final_floor_coordinates(LINE_LIST *beamline, double *XYZ, double *Angle,
                             double *XYZMin, double *XYZMax)
{
  ELEMENT_LIST *elem;
  double X, Y, Z, theta, phi, psi, s;
  MATRIX *V0, *V1;
  MATRIX *Theta, *Phi, *Psi, *W0, *W1, *temp;
  long i;

  elem = &(beamline->elem);
  X = X0;
  Y = Y0;
  Z = Z0;
  s = 0;
  
  m_alloc(&V0, 3, 1);
  m_alloc(&V1, 3, 1);
  V0->a[0][0] = X0;
  V0->a[1][0] = Y0;
  V0->a[2][0] = Z0;

  m_alloc(&Theta, 3, 3);
  m_alloc(&Phi, 3, 3);
  m_alloc(&Psi, 3, 3);
  m_zero(Theta);
  m_zero(Phi);
  m_zero(Psi);

  theta = theta0;
  Theta->a[0][0] = Theta->a[2][2] = cos(theta0);
  Theta->a[0][2] = -(Theta->a[2][0] = -sin(theta0));
  Theta->a[1][1] = 1;

  phi = phi0;
  Phi->a[1][1] = Phi->a[2][2] = cos(phi0);
  Phi->a[1][2] = -(Phi->a[2][1] = -sin(phi0));
  Phi->a[0][0] = 1;

  psi = psi0;
  Psi->a[0][0] = Psi->a[1][1] = cos(psi0);
  Psi->a[0][1] = -(Psi->a[1][0] = sin(psi0));
  Psi->a[2][2] = 1;
  
  m_alloc(&W0, 3, 3);
  m_alloc(&temp, 3, 3);
  m_mult(temp, Theta, Phi);
  m_mult(W0, temp, Psi);
  m_alloc(&W1, 3, 3);

  if (XYZMin) {
    XYZMin[0] = X0;
    XYZMin[1] = Y0;
    XYZMin[2] = Z0;
  }
  if (XYZMax) {
    XYZMax[0] = X0;
    XYZMax[1] = Y0;
    XYZMax[2] = Z0;
  }
      
  while (elem) {
    advanceFloorCoordinates(V1, W1, V0, W0, &theta, &phi, &psi, &s, elem, NULL, NULL, 0);
    if (elem->type!=T_FLOORELEMENT) {
      long i;
      if (XYZMin)
        for (i=0; i<3; i++)
          XYZMin[i] = MIN(V1->a[i][0], XYZMin[i]);
      if (XYZMax)
        for (i=0; i<3; i++)
          XYZMax[i] = MAX(V1->a[i][0], XYZMax[i]);
    } else {
      long i;
      if (XYZMin)
        for (i=0; i<3; i++)
          XYZMin[i] = V1->a[i][0];
      if (XYZMax)
        for (i=0; i<3; i++)
          XYZMax[i] = V1->a[i][0];
    }
    if (elem->type==T_MARK && ((MARK*)elem->p_elem)->fitpoint) {
      MARK *mark;
      char t[100];
      static char *suffix[7] = {"X", "Y", "Z", "theta", "phi", "psi", "s"};
      mark = (MARK*)(elem->p_elem);
      if (!(mark->init_flags&4)) {
        mark->floor_mem = tmalloc(sizeof(*mark->floor_mem)*7);
        for (i=0; i<7; i++) {
          sprintf(t, "%s#%ld.%s", elem->name, elem->occurence,
                  suffix[i]);
          mark->floor_mem[i] = rpn_create_mem(t, 0);
        }
        mark->init_flags |= 4;
      }
      rpn_store(V1->a[0][0], NULL, mark->floor_mem[0]);
      rpn_store(V1->a[1][0], NULL, mark->floor_mem[1]);
      rpn_store(V1->a[2][0], NULL, mark->floor_mem[2]);
      rpn_store(theta, NULL, mark->floor_mem[3]);
      rpn_store(phi, NULL, mark->floor_mem[4]);
      rpn_store(psi, NULL, mark->floor_mem[5]);
      rpn_store(s, NULL, mark->floor_mem[6]);
    }
    m_copy(V0, V1);
    m_copy(W0, W1);
    elem->end_pos = s;
    elem = elem->succ;
  }
  beamline->revolution_length = s;
  for (i=0; i<3; i++)
    XYZ[i] = V0->a[i][0];
  Angle[0] = theta;
  Angle[1] = phi;
  Angle[2] = psi;
  m_free(&V0);
  m_free(&V1);
  m_free(&Theta);
  m_free(&Phi);
  m_free(&Psi);
  m_free(&W0);
  m_free(&temp);
  m_free(&W1);
}


double nearbyAngle(double angle, double reference)
{
  double minDiff, bestAngle, diff;
  bestAngle = angle;
  minDiff = fabs(reference-angle);

  angle += PIx2;
  diff = fabs(reference-angle);
  if (diff<minDiff) {
    minDiff = diff;
    bestAngle = angle;
  }

  angle -= 2*PIx2;
  diff = fabs(reference-angle);
  if (diff<minDiff) {
    minDiff = diff;
    bestAngle = angle;
  }
  return bestAngle;
}

void computeSurveyAngles(double *theta, double *phi, double *psi, MATRIX *W)
{
  double arg;
  
  arg = sqrt( sqr(W->a[1][0]) + sqr(W->a[1][1]));  /* |cos(phi)| */
  arg = SIGN(cos(*phi))*arg;
  *phi = nearbyAngle(atan2(W->a[1][2], arg), *phi);
  if (fabs(arg)>1e-15) {
    *theta = nearbyAngle(atan2(W->a[0][2], W->a[2][2]), *theta);
    *psi   = nearbyAngle(atan2(W->a[1][0], W->a[1][1]), *psi);
  }
  else {
    *psi = nearbyAngle(atan2(-W->a[0][1], W->a[0][0])-*theta, *psi);
  }
  
}
