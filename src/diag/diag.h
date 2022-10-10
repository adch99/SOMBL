#ifndef MBL_DIAG_H
#define MBL_DIAG_H

#include "../constants.h"

int diag_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals);
int diag_get_eigh(CDTYPE * matrix, int size, DTYPE * eigvals);

#endif