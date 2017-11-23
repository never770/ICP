//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/* 
 * June 2001, Johnny Park
 * 
 * June 2003,
 *
*/

/*	
 *	Collection of useful and often used functions.
 *
 *	Started Jan 2001, Johnny Park
 */
#ifdef __cplusplus
extern "C" {
#endif
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <f2c.h>
#include <clapack.h>
#include "point3D.h"
#include "Datatype.h"


#ifndef   SQ
#define   SQ(x)   ((x) * (x))
#endif

#ifndef   PI
#define   PI      3.1415926535897932   
#endif  

/* 
 * Data definitions
 */
typedef struct vector {
  int n;
  double *entry;
} vector;

typedef struct matrix {
  int r;
  int c;
  double **entry;
} matrix;
typedef struct point_set_matrix {
	matrix *point_set;
} point_set_matrix;

/*
 * Function Declaration
 */
 //打印矩阵
void PrintMatrix(matrix *M);
//打印向量
void PrintVector(vector *p);
//计算行列式
double Determinant33(matrix *A);
//矩阵转置
void Transpose(matrix *M1, matrix *M2);
void SetMatrix(matrix *M, double *v);
void GetMatrix(matrix *M, double *v);
void SetMatrixRowMajor(matrix *M, double *v);
void GetMatrixRowMajor(matrix *M, double *v);
matrix *AllocateMatrix(int row, int col);
vector *AllocateVector(int num);
point_set_matrix *Allocatpoint_set_matrix(int row, int col);
void FreePoint_set_matrix(point_set_matrix *p);
void PutPointSetIntoMatrix(const point_xyz_set *point_set_A, point_set_matrix *point_set_Matrix);
void PutPointIntoMatrix(const point_xyz *A, matrix *pointmatrix);
void FreeMatrix(matrix *M);
void FreeVector(vector *v);
void CopyMatrix(matrix *M1, matrix *M2);
void CopyVector(vector *v1, vector *v2);
void MultMatrix(matrix *M1, matrix *M2, matrix *M3);
void AddMatrix(matrix *M1, matrix *M2, matrix *M3);
void AddVector(vector *p1, vector *p2, vector *p3);
void SubMatrix(matrix *M1, matrix *M2, matrix *M3);
void SubVector(vector *p1, vector *p2, vector *p3);
void CrossProduct(vector *p1, vector *p2, vector *p3);
double DotProduct(vector *p1, vector *p2);
void UnitVector(vector *p1, vector *p2);
void MultConstantToMatrix(double x, matrix *M1, matrix *M2);
double Trace(matrix *M);
void RotationAboutVector(double theta, double rx, double ry, double rz,
	matrix *R);
void AntiSymmetric(matrix *M1, matrix *M2);
void RotationQuaternion(vector *p, matrix *R);
void QuaternionFromRotation(matrix *R, vector *p);
int  SymmetricEigens(matrix *A, matrix *Q, vector *v);
void SymmetricLargestEigens(matrix *A, vector *evec, double *eval);
void SymmetricSmallestEigens(matrix *A, vector *evec, double *eval);
void InverseMatrix(matrix *M1, matrix *M2);
void initia_matrix4f(Matrix4f center_A, char *type);
void Initia_matrix(matrix *mat, char *type);
void PrintMatrix4f(Matrix4f a);

//4*4转移矩阵和点集相乘
void Matrix4fMultPoint_xyz_set(const Matrix4f trans, const point_xyz_set *ps, point_xyz_set *MultResult);

/*
 * clapack functions
 */
 
typedef long int integer;
typedef double doublereal;

int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	         integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	         integer *info);
	         
int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
	          doublereal *a, integer *lda, doublereal *vl, doublereal *vu, 
	          integer *il, integer *iu, doublereal *abstol, integer *m, 
	          doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	          integer *lwork, integer *iwork, integer *ifail, integer *info);
	
int dgesv_(integer *n, integer *nrhs, doublereal *a, integer *lda, 
           integer *ipiv, doublereal *b, integer *ldb, integer *info);



#endif

#ifdef __cplusplus
}
#endif