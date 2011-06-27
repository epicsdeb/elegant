/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: response.c
 * purpose: trajectory/orbit correction response matrix output for elegant
 * 
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "correctDefs.h"

extern char *correction_mode[2];

typedef struct {
    SDDS_TABLE SDDSout;
    long monitors, correctors;
    char **monitorName;
    long sIndex, monitorNameIndex, *correctorIndex;
    } RESPONSE_OUTPUT;
RESPONSE_OUTPUT xRespOutput, yRespOutput, xInvRespOutput, yInvRespOutput;
RESPONSE_OUTPUT xyRespOutput, yxRespOutput;

#define NORMAL_UNITS 0
#define KNL_UNITS 1
#define BNL_UNITS 2

static CORMON_DATA CMyx, CMxy;

void setup_response_output(RESPONSE_OUTPUT *respOutput,
                           char *filename, char *type, RUN *run, char *beamline_name, CORMON_DATA *CM, STEERING_LIST *SL, 
                           long bpmPlane, long corrPlane, long inverse, long unitsType);
void do_response_output(RESPONSE_OUTPUT *respOutput, CORMON_DATA *CM, STEERING_LIST *SL, long plane,
                   long inverse, long KnL_units, long tune_corrected);
#include "response.h"

void setup_correction_matrix_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CORRECTION *correct,
                                    long *do_response,
                                    long BnLUnitsOK)
{
    long unitsCode;
    log_entry("setup_correction_matrix_output");
    if (correct->mode==-1) {
      printf("Error: you must request orbit/trajectory correction before requesting response matrix output.\n");
      printf("Otherwise, elegant doesn't have the needed information about correctors and monitors.\n");
      exitElegant(1);  
    }
    

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&correction_matrix_output, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &correction_matrix_output);
    *do_response = output_at_each_step;

    unitsCode = KnL_units?KNL_UNITS:(BnL_units?BNL_UNITS:0);
    if (unitsCode==BNL_UNITS && !BnLUnitsOK)
      bombElegant("At present you must give the matrix_output or twiss_output command to use BnL_units=1.  Sorry.", NULL);

    if (coupled) {
      /* set up the coupled response matrix data structures */

      /* CMyx is the vertical response to horizontal correctors */
      /* --- copy corrector information */
      CMyx.ncor = correct->CMx->ncor;
      CMyx.kick_coef = correct->CMx->kick_coef;
      CMyx.sl_index = correct->CMx->sl_index;
      CMyx.ucorr = correct->CMx->ucorr;
      /* --- copy BPM information */      
      CMyx.nmon = correct->CMy->nmon;
      CMyx.umoni = correct->CMy->umoni;
      CMyx.mon_index = correct->CMy->mon_index;
      CMyx.C = NULL;
      CMyx.bpmPlane = 1;
      CMyx.corrPlane = 0;
      
      /* CMxy is the horizontal response to vertical correctors */
      /* --- copy corrector information */
      CMxy.ncor = correct->CMy->ncor;
      CMxy.kick_coef = correct->CMy->kick_coef;
      CMxy.sl_index = correct->CMy->sl_index;
      CMxy.ucorr = correct->CMy->ucorr;
      CMxy.fixed_length = fixed_length;
      
      /* --- copy BPM information */      
      CMxy.nmon = correct->CMx->nmon;
      CMxy.umoni = correct->CMx->umoni;
      CMxy.mon_index = correct->CMx->mon_index;
      CMxy.C = NULL;
      CMxy.bpmPlane = 0;
      CMxy.corrPlane = 1;
      CMxy.fixed_length = fixed_length;
    }
    
#if USE_MPI
    if (isSlave)
      return;
#endif

    if (response[0])
        setup_response_output(&xRespOutput, response[0], correction_mode[correct->mode], run, beamline->name,
                              correct->CMx, &correct->SLx, 0, 0, 0, unitsCode);
    if (response[1])
        setup_response_output(&yRespOutput, response[1], correction_mode[correct->mode], run, beamline->name,
                              correct->CMy, &correct->SLy, 1, 1, 0, unitsCode);

    if (response[2])
        setup_response_output(&yxRespOutput, response[2], correction_mode[correct->mode], run, beamline->name,
                              &CMyx, &correct->SLx, 1, 0, 0, unitsCode);
    if (response[3])
        setup_response_output(&xyRespOutput, response[3], correction_mode[correct->mode], run, beamline->name,
                              &CMxy, &correct->SLy, 0, 1, 0, unitsCode);
      
    if (inverse[0])
        setup_response_output(&xInvRespOutput, inverse[0], correction_mode[correct->mode], run, beamline->name,
                              correct->CMx, &correct->SLx, 0, 0, 1, unitsCode);
    if (inverse[1])
        setup_response_output(&yInvRespOutput, inverse[1], correction_mode[correct->mode], run, beamline->name,
                              correct->CMy, &correct->SLy, 1, 1, 1, unitsCode);

    log_exit("setup_correction_matrix_output");
    }

void setup_response_output(RESPONSE_OUTPUT *respOutput,
                           char *filename, char *type, RUN *run, char *beamline_name, CORMON_DATA *CM, STEERING_LIST *SL, 
                           long bpmPlane, long correctorPlane, long inverse, long unitsCode)
{
    ELEMENT_LIST *eptr;
    static char s[256], t[256], units[32];
    long i, j, *unique_name, sl_index;
    long matrixType;
    char *matrixTypeName[4] = {"horizontal", "vertical", "v->h", "h->v"};

    matrixType = bpmPlane==correctorPlane ? bpmPlane : 2+bpmPlane;
    
    log_entry("setup_response_output");
    filename = compose_filename(filename, run->rootname);
    if (!inverse) {
        sprintf(s, "%s-plane %s %sfixed path-length response matrix for beamline %s of lattice %s", 
                matrixTypeName[matrixType], type,
                fixed_length?"":"non-", beamline_name, run->lattice);
        sprintf(t, "%s-plane %s %sfixed path-length response matrix", matrixTypeName[matrixType], type,
                fixed_length?"":"non-");
        }
    else {
        sprintf(s, "%s-plane %s transposed inverse %sfixed path-length response matrix for beamline %s of lattice %s", 
                matrixTypeName[matrixType], type, 
                fixed_length?"":"non-", beamline_name, run->lattice);
        sprintf(t, "%s-plane %s transposed inverse %sfixed path-length response matrix", matrixTypeName[matrixType],
                type, fixed_length?"":"non-");
        }

    if (!SDDS_InitializeOutput(&respOutput->SDDSout, SDDS_BINARY, 0, s, t, filename)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
        }
#if RW_ASSOCIATES
    SDDS_DefineAssociate(&respOutput->SDDSout, "elegantInput",
                         run->runfile, getenv("PWD"), "elegant input file used to create this file", 
                         "elegant input, parent", 0);
    SDDS_DefineAssociate(&respOutput->SDDSout, "elegantLattice",
                         run->lattice, getenv("PWD"), "elegant lattice file used to create this file", 
                         "elegant lattice, parent", 0);
#endif
    respOutput->monitorNameIndex = 
        SDDS_DefineColumn(&respOutput->SDDSout, "BPMName", NULL, NULL, "beam-position-monitor name", NULL, SDDS_STRING, 0);
    respOutput->sIndex = 
        SDDS_DefineColumn(&respOutput->SDDSout, "s", NULL, "m", "beam-position-monitor location", NULL, SDDS_DOUBLE, 0);
    SDDS_DefineParameter(&respOutput->SDDSout, "CorrectionMatrixType", NULL, NULL, "correction matrix type",
                         NULL, SDDS_STRING, inverse?"Inverse":"Response");
    SDDS_DefineParameter1(&respOutput->SDDSout, "NCorrectors", NULL, NULL, "Number of correctors", NULL, SDDS_LONG, 
                         &CM->ncor);
    SDDS_DefineParameter1(&respOutput->SDDSout, "NMonitors", NULL, NULL, "Number of correctors", NULL, SDDS_LONG, 
                          &CM->nmon);
    SDDS_DefineParameter(&respOutput->SDDSout, "Stage", NULL, NULL, "Simulation stage", NULL, SDDS_STRING, NULL);
    if (SDDS_NumberOfErrors())
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

    respOutput->correctors = CM->ncor;
    respOutput->monitors = CM->nmon;
    respOutput->correctorIndex = tmalloc(sizeof(*respOutput->correctorIndex)*respOutput->correctors);
    respOutput->monitorName = tmalloc(sizeof(*respOutput->monitorName)*respOutput->monitors);

    unique_name = tmalloc(sizeof(*unique_name)*(CM->ncor));
    for (i=0; i<CM->ncor; i++)
        unique_name[i] = 1;
    for (i=0; i<CM->ncor; i++) {
        for (j=0; j<i; j++) {
            if (!unique_name[j])
                continue;
            if (strcmp(CM->ucorr[i]->name, CM->ucorr[j]->name)==0)
                unique_name[i] = unique_name[j] = 0;
            }
        }

    for (i=0; i<CM->ncor; i++) {
        sl_index = CM->sl_index[i];
        eptr = CM->ucorr[i];
        if (!inverse) {
            if (is_blank(entity_description[eptr->type].parameter[SL->param_index[sl_index]].unit))
                strcpy_ss(units, "m");
            else {
                if (unitsCode==KNL_UNITS && (eptr->type==T_HCOR || eptr->type==T_HVCOR || eptr->type==T_VCOR))
                    sprintf(units, "m/K0L");
                else if (unitsCode==BNL_UNITS && (eptr->type==T_HCOR || eptr->type==T_HVCOR || eptr->type==T_VCOR))
                    sprintf(units, "1/T");
                else {
                    sprintf(units, "m/%s", entity_description[eptr->type].parameter[SL->param_index[sl_index]].unit);
                    str_tolower(units);
                    }
                }
            sprintf(s, "response to %s#%ld.%s", eptr->name, eptr->occurence, SL->corr_param[sl_index]);
            }
        else {
            if (is_blank(entity_description[eptr->type].parameter[SL->param_index[sl_index]].unit))
                strcpy_ss(units, "1/m");
            else {
                if (unitsCode==KNL_UNITS 
                    && (eptr->type==T_HCOR || eptr->type==T_HVCOR || eptr->type==T_VCOR))
                  sprintf(units, "K0L/m");
                else if (unitsCode==BNL_UNITS 
                         && (eptr->type==T_HCOR || eptr->type==T_HVCOR || eptr->type==T_VCOR))
                  sprintf(units, "T");
                else {
                    sprintf(units, "%s/m", entity_description[eptr->type].parameter[SL->param_index[sl_index]].unit);
                    str_tolower(units);
                    }
                }
            sprintf(s, "inverse response to %s#%ld.%s", eptr->name, eptr->occurence, SL->corr_param[sl_index]);
            }
        if (unique_name[i])
            respOutput->correctorIndex[i] = 
                SDDS_DefineColumn(&respOutput->SDDSout, CM->ucorr[i]->name, NULL, units, s, NULL, SDDS_DOUBLE, 0);
        else {
            sprintf(t, "%s#%ld", CM->ucorr[i]->name, CM->ucorr[i]->occurence);
            respOutput->correctorIndex[i] = 
                SDDS_DefineColumn(&respOutput->SDDSout, t, NULL, units, s, NULL, SDDS_DOUBLE, 0);
            }
        }
    if (SDDS_NumberOfErrors()) {
        SDDS_PrintErrors(stderr, 1);
        exitElegant(1);
        }
    if (!SDDS_WriteLayout(&respOutput->SDDSout)) {
        SDDS_PrintErrors(stderr, 1);
        exitElegant(1);
        }

    unique_name = trealloc(unique_name, sizeof(*unique_name)*(CM->nmon));
    for (i=0; i<CM->nmon; i++)
        unique_name[i] = 1;
    for (i=0; i<CM->nmon; i++) {
        for (j=0; j<i; j++) {
            if (!unique_name[j])
                continue;
            if (strcmp(CM->umoni[i]->name, CM->umoni[j]->name)==0)
                unique_name[i] = unique_name[j] = 0;
            }
        }

    for (i=0; i<CM->nmon; i++) {
        if (!unique_name[i]) {
            sprintf(s, "%s#%ld", CM->umoni[i]->name, CM->umoni[i]->occurence);
            SDDS_CopyString(&respOutput->monitorName[i], s);
            }
        else
            SDDS_CopyString(&respOutput->monitorName[i], CM->umoni[i]->name);
        }
    set_namelist_processing_flags(0);
    set_print_namelist_flags(0);
    log_exit("setup_response_output");
    }

/* This routine is only run by the optimizer.  It results in values being placed in
 * rpn memories for use in optimization.
 */

void update_response(RUN *run, LINE_LIST *beamline, CORRECTION *correct) 
{
    MAT *Cx, *Cy, *Tx, *Ty;

    /* Copy the matrices from the correction structure so we can put them back when we are done. 
     * We have to do this because the correction command may have different settings (e.g., fixed-length
     * constraint or non-perturbed matrix.
     */

    Cx = correct->CMx->C;
    Tx = correct->CMx->T;
    correct->CMx->C = correct->CMx->T = NULL;

    Cy = correct->CMy->C;
    Ty = correct->CMy->T;
    correct->CMy->C = correct->CMy->T = NULL;

    if (correct->mode==TRAJECTORY_CORRECTION) {
      compute_trajcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, COMPUTE_RESPONSE_SILENT);
      compute_trajcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, COMPUTE_RESPONSE_SILENT);
      if (coupled) {
        compute_trajcor_matrices(&CMxy, &correct->SLy, 0, run, beamline, COMPUTE_RESPONSE_SILENT);
        compute_trajcor_matrices(&CMyx, &correct->SLx, 2, run, beamline, COMPUTE_RESPONSE_SILENT);
      }
    }
    else if (correct->mode==ORBIT_CORRECTION) {
      compute_orbcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, COMPUTE_RESPONSE_SILENT|(fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
      compute_orbcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, COMPUTE_RESPONSE_SILENT|(fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
      if (coupled) {
        compute_orbcor_matrices1(&CMxy, &correct->SLy, 0, run, beamline, COMPUTE_RESPONSE_SILENT|(fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
        compute_orbcor_matrices1(&CMyx, &correct->SLx, 2, run, beamline, COMPUTE_RESPONSE_SILENT|(fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
      }
    }
    else
      bombElegant("bad correction mode (update_response)", NULL);

    /* copy matrices back to the correction structure and free memory */
    matrix_free(correct->CMx->C);
    matrix_free(correct->CMx->T);
    correct->CMx->C = Cx;
    correct->CMx->T = Tx;

    matrix_free(correct->CMy->C);
    matrix_free(correct->CMy->T);
    correct->CMy->C = Cy;
    correct->CMy->T = Ty;

  }



void run_response_output(RUN *run, LINE_LIST *beamline, CORRECTION *correct, long tune_corrected)
{
    long unitsCode;
    MAT *Cx, *Cy, *Tx, *Ty;

    unitsCode = KnL_units?KNL_UNITS:(BnL_units?BNL_UNITS:0);
    if (tune_corrected==0 && !output_before_tune_correction)
        return;

    /* Copy the matrices from the correction structure so we can put them back when we are done. 
     * We have to do this because the correction command may have different settings (e.g., fixed-length
     * constraint or non-perturbed matrix.
     */

    Cx = correct->CMx->C;
    Tx = correct->CMx->T;
    correct->CMx->C = correct->CMx->T = NULL;

    Cy = correct->CMy->C;
    Ty = correct->CMy->T;
    correct->CMy->C = correct->CMy->T = NULL;

    if (correct->mode==TRAJECTORY_CORRECTION) {
      printf("Computing trajectory correction matrices for output.\n");
      fflush(stdout);
      compute_trajcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 
                               (!(inverse[0]==NULL || SDDS_StringIsBlank(inverse[0])) ? COMPUTE_RESPONSE_INVERT : 0)|COMPUTE_RESPONSE_SILENT);
      compute_trajcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, 
                               (!(inverse[1]==NULL || SDDS_StringIsBlank(inverse[1])) ? COMPUTE_RESPONSE_INVERT : 0)|COMPUTE_RESPONSE_SILENT);

      if (coupled) {
        compute_trajcor_matrices(&CMxy, &correct->SLy, 0, run, beamline, COMPUTE_RESPONSE_SILENT);
        compute_trajcor_matrices(&CMyx, &correct->SLx, 2, run, beamline, COMPUTE_RESPONSE_SILENT);
      }
    }
    else if (correct->mode==ORBIT_CORRECTION) {
      printf("Computing orbit correction matrices for output, with %s length.\n",
             fixed_length?"fixed":"variable");
      fflush(stdout);
      compute_orbcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 
                              (!(inverse[0]==NULL || SDDS_StringIsBlank(inverse[0])) ? COMPUTE_RESPONSE_INVERT : 0)|COMPUTE_RESPONSE_SILENT|
                              (fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
      compute_orbcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, 
                              (!(inverse[1]==NULL || SDDS_StringIsBlank(inverse[1])) ? COMPUTE_RESPONSE_INVERT : 0)|COMPUTE_RESPONSE_SILENT|
                              (fixed_length?COMPUTE_RESPONSE_FIXEDLENGTH:0));
      if (coupled) {
        compute_orbcor_matrices1(&CMxy, &correct->SLy, 0, run, beamline, COMPUTE_RESPONSE_SILENT);
        compute_orbcor_matrices1(&CMyx, &correct->SLx, 2, run, beamline, COMPUTE_RESPONSE_SILENT);
      }
    }
    else
      bombElegant("bad correction mode (run_response_output)", NULL);

#if USE_MPI
    if (!isSlave) {
#endif

      if (response[0])
        do_response_output(&xRespOutput, correct->CMx, &correct->SLx, 0, 0, unitsCode, tune_corrected);
      if (response[1])
        do_response_output(&yRespOutput, correct->CMy, &correct->SLy, 1, 0, unitsCode, tune_corrected);
      if (response[2])
        do_response_output(&yxRespOutput, &CMyx, &correct->SLx, 0, 0, unitsCode, tune_corrected);
      if (response[3])
        do_response_output(&xyRespOutput, &CMxy, &correct->SLy, 1, 0, unitsCode, tune_corrected);
      if (inverse[0])
        do_response_output(&xInvRespOutput, correct->CMx, &correct->SLx, 0, 1, unitsCode, tune_corrected);
      if (inverse[1])
        do_response_output(&yInvRespOutput, correct->CMy, &correct->SLy, 1, 1, unitsCode, tune_corrected);

#if USE_MPI
  }
#endif
    
    /* copy matrices back to the correction structure and free memory */
    matrix_free(correct->CMx->C);
    matrix_free(correct->CMx->T);
    correct->CMx->C = Cx;
    correct->CMx->T = Tx;

    matrix_free(correct->CMy->C);
    matrix_free(correct->CMy->T);
    correct->CMy->C = Cy;
    correct->CMy->T = Ty;

    }

void do_response_output(RESPONSE_OUTPUT *respOutput, CORMON_DATA *CM, STEERING_LIST *SL, long corrPlane,
                   long inverse, long unitsCode, long tune_corrected)
{
    long i, j;
    ELEMENT_LIST *eptr;
    double value;

    if (!SDDS_StartTable(&respOutput->SDDSout, CM->nmon) ||
        !SDDS_SetParameters(&respOutput->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                            "Stage", 
                            tune_corrected==0?"tune uncorrected":
                            (tune_corrected==1?"tune corrected":"input lattice response"), NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<CM->nmon; i++) {
        if (!SDDS_SetRowValues(&respOutput->SDDSout, SDDS_PASS_BY_VALUE|SDDS_SET_BY_INDEX, i,
                               respOutput->monitorNameIndex, respOutput->monitorName[i],
                               respOutput->sIndex, CM->umoni[i]->end_pos, -1))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (unitsCode==BNL_UNITS) 
            for (j=0; j<CM->ncor; j++) {
                eptr = CM->ucorr[j];
                value = (inverse?Mij(CM->T, j, i):Mij(CM->C, i, j));
                if (eptr->type==T_HCOR || eptr->type==T_HVCOR || eptr->type==T_VCOR) {
                    value *= (inverse?eptr->Pref_output/586.679:586.679/(eptr->Pref_output+1e-10));
                  }
                if (!SDDS_SetRowValues(&respOutput->SDDSout, SDDS_PASS_BY_VALUE|SDDS_SET_BY_INDEX, i,
                                       respOutput->correctorIndex[j], value, -1))
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        else if (unitsCode==KNL_UNITS && corrPlane==0) 
            for (j=0; j<CM->ncor; j++) {
                eptr = CM->ucorr[j];
                value = (inverse?Mij(CM->T, j, i):Mij(CM->C, i, j));
                if (eptr->type==T_HCOR || eptr->type==T_HVCOR)
                    value = -value;
                if (!SDDS_SetRowValues(&respOutput->SDDSout, SDDS_PASS_BY_VALUE|SDDS_SET_BY_INDEX, i,
                                       respOutput->correctorIndex[j], value, -1))
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        else
            for (j=0; j<CM->ncor; j++) {
                value = (inverse?Mij(CM->T, j, i):Mij(CM->C, i, j));
                if (!SDDS_SetRowValues(&respOutput->SDDSout, SDDS_PASS_BY_VALUE|SDDS_SET_BY_INDEX, i,
                                       respOutput->correctorIndex[j], value, -1))
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                }
        }
    if (!SDDS_WriteTable(&respOutput->SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

void finish_response_output(void)
{
    log_entry("finish_response_output");
    if (response[0] && !SDDS_Terminate(&xRespOutput.SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (response[1] && !SDDS_Terminate(&yRespOutput.SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (inverse[0] && !SDDS_Terminate(&xInvRespOutput.SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (inverse[1] && !SDDS_Terminate(&yInvRespOutput.SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    log_exit("finish_response_output");
    }
