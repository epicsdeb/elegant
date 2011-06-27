/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: subprocess.c
 *
 * Michael Borland, 1994
 */
#include <stdio.h>
#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif
#include "mdb.h"
#include "track.h"
#include "subprocess.h"
#include <signal.h>


/* dummy signal handler for use with sigpause */
void subprocess_sigusr1()
{
}

void run_subprocess(NAMELIST_TEXT *nltext, RUN *run)
{
  static char buffer[1024];
  char *ptr, *ptr0;

  log_entry("run_subprocess");

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&subprocess, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &subprocess);

  if (command) {
    buffer[0] = 0;
    ptr0 = command;
    while ((ptr=strstr(ptr0, "%s"))) {
      if (ptr!=command && *(ptr-1)=='%') {
        *(ptr-1) = 0;
        strcat(buffer, ptr0);
        strcat(buffer, "%s");
        ptr += 2;
        ptr0 = ptr;
      }
      else {
        if (!run || !run->rootname)
          bombElegant("rootname must be initialized prior to subprocess execution if \%s substitution is used", NULL);
        *ptr = 0;
        ptr += 2;
        strcat(buffer, ptr0);
        strcat(buffer, run->rootname);
        ptr0 = ptr;
      }
    }
    strcat(buffer, ptr0);
    fprintf(stdout, "%s\n", buffer);
    fflush(stdout);
    executeCshCommand(buffer);
  }

  log_exit("run_subprocess");
}

void executeCshCommand(char *cmd)
{
  FILE *fp;
  char *filename;
  char cmd2[1000];
  
#if defined(CONDOR_COMPILE)
  _condor_ckpt_disable();
#endif

  filename = tmpname(NULL);
  fp = fopen(filename, "w");
  fprintf(fp, "set nonomatch\n");
  fprintf(fp, "%s\n", cmd);
  fclose(fp);
  sprintf(cmd2, "csh %s\n", filename);
  system(cmd2);
  remove(filename);

#if defined(CONDOR_COMPILE)
  _condor_ckpt_enable();
#endif
}

