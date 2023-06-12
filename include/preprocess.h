#ifndef __PREPROCESS__
#define __PREPROCESS__

#include <nics_config.h>
#include <nicslu.h>
#include "type.h"

#ifdef __cplusplus
extern "C" {
#endif

int preprocess(char *matrixName, SNicsLU *nicslu, double **ax,
               unsigned int **ai, unsigned int **ap);

#ifdef __cplusplus
}
#endif

#endif
