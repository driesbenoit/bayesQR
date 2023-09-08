#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(qrb_al_mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(qrb_mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(qrc_al_mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(qrc_mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"qrb_al_mcmc", (DL_FUNC) &F77_NAME(qrb_al_mcmc), 10},
    {"qrb_mcmc",    (DL_FUNC) &F77_NAME(qrb_mcmc),    10},
    {"qrc_al_mcmc", (DL_FUNC) &F77_NAME(qrc_al_mcmc), 14},
    {"qrc_mcmc",    (DL_FUNC) &F77_NAME(qrc_mcmc),    14},
    {NULL, NULL, 0}
};

void R_init_bayesQR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
