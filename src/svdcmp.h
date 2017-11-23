#pragma once
#ifndef _SVDCMP_H
#define _SVDCMP_H
#ifdef __cplusplus
extern "C" {
#endif //  _cplusplus
#include "point3D.h"
#include "matrix.h"

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

double **dmatrix(int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);
double pythag(double a, double b);
void svdcmp(double **a, int m, int n, double w[], double **v);
void rigid_transform3D(point_xyz_set *A, point_xyz_set *B, matrix *R, matrix *trans);

#ifdef __cplusplus
}
#endif // _cplusplus

#endif // !_SVDCMP_H
