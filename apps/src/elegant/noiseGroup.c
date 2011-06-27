/*************************************************************************\
* Copyright (c) 2006 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2006 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"

static long *idValue = NULL;
static short *valueValid = NULL;
static double *randomValue = NULL;
static long groups = 0, maxGroups = 0;

long DefineNoiseGroup(long groupId) 
{
  long i;
  for (i=0; i<groups; i++) {
    if (groupId==idValue[i])
      break;
  }
  if (i==groups) {
    groups ++;
    printf("Defining noise group %ld\n", groupId);
    fflush(stdout);
    if (groups>=maxGroups) {
      maxGroups += 10;
      printf("Resizing to maximum of %ld noise groups\n", maxGroups);
      fflush(stdout);
      if (!(idValue = SDDS_Realloc(idValue, sizeof(*idValue)*maxGroups)) ||
	  !(valueValid = SDDS_Realloc(valueValid, sizeof(*randomValue)*maxGroups)) ||
	  !(randomValue = SDDS_Realloc(randomValue, sizeof(*randomValue)*maxGroups))) {
	return -1;
      }
    }
    idValue[i] = groupId;
    randomValue[i] = DBL_MAX;
    valueValid[i] = 0;
    return 1;
  }
  return 0;
}

long ResetNoiseGroupValues()
{
  long i;
  for (i=0; i<groups; i++) {
    valueValid[i] = 0;
  }
  return groups;
}

double GetNoiseGroupValue(long groupId) 
{
  long i;
  TRACKING_CONTEXT tContext;
  for (i=0; i<groups; i++) {
    if (groupId==idValue[i]) {
      getTrackingContext(&tContext);
      if (!valueValid[i]) {
	randomValue[i] = gauss_rn_lim(0.0, 1.0, 2, random_3);
	printf("Noise value for group %ld is %e\n", idValue[i], randomValue[i]);
      }
      printf("Returning noise value %e for %s\n", randomValue[i], tContext.elementName);
      fflush(stdout);
      valueValid[i] = 1;
      return randomValue[i];
    }
  }
  getTrackingContext(&tContext);
  printf("Error: noise group value requested for nonexistent group %ld for %s\n", groupId,
	 tContext.elementName);
  return DBL_MAX;
}
