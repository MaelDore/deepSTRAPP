#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void setrecursivesequence(void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP cache_descendants(SEXP);

static const R_CMethodDef CEntries[] = {
    {"setrecursivesequence", (DL_FUNC) &setrecursivesequence, 6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"cache_descendants",   (DL_FUNC) &cache_descendants,   1},
    {NULL, NULL, 0}
};

void R_init_deepSTRAPP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}