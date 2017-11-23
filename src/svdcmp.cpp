/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "svdcmp.h"
#define dubug 0

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	/*
	nrl   行的起始位置
	nrh   行的终止位置
	ncl  列的起始位置
	nch  列的终止位置
	*/

	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;
	/* allocate pointers to rows */    //申请行的
	m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (double *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
	return v - nl + NR_END;
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.0 + (absb / absa)*(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + (absa / absb)*(absa / absb)));
}

/******************************************************************************/
void svdcmp(double **a, int m, int n, double w[], double **v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	rv1 = dvector(1, n);
	g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale *g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s*rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = DMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) { /* Accumulation of right-hand transformations. */
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++) /* Double division to avoid possible underflow. */
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) v[k][j] += s*v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--) { /* Accumulation of left-hand transformations. */
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i])*g;
				for (k = i; k <= m; k++) a[k][j] += f*a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k = n; k >= 1; k--) { /* Diagonalization of the bidiagonal form. */
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) { /* Test for splitting. */
				nm = l - 1; /* Note that rv1[1] is always zero. */
				if ((double)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((double)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0; /* Cancellation of rv1[l], if l > 1. */
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((double)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g*h;
					s = -f*h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c + z*s;
						a[j][i] = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x = w[l]; /* Shift from bottom 2-by-2 minor. */
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h*((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0; /* Next QR transformation: */
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c + z*s;
					v[jj][i] = z*c - x*s;
				}
				z = pythag(f, h);
				w[j] = z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z = 1.0 / z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c + z*s;
					a[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_dvector(rv1, 1, n);
}
//对A进行变换，使之与B重合
void rigid_transform3D(point_xyz_set *A, point_xyz_set *B, matrix *R, matrix *trans)
{

	//printf("----------------------p_base---------------\n");
	//PrintPointXYZSet(A);
	//printf("-------------------------------------\n");
	//PrintPointXYZSet(B);
    //求点云集的质心
	point_xyz centroid_A, centroid_B;
	centroid_A = CenterOfMass(A->head, A->num);
	centroid_B = CenterOfMass(B->head, B->num);
	matrix *centroid_A_matrix;                       
	centroid_A_matrix = AllocateMatrix(1, 3);
	PutPointIntoMatrix(&centroid_A, centroid_A_matrix);
	matrix *centroid_B_matrix;                                 
	centroid_B_matrix = AllocateMatrix(1, 3);
	PutPointIntoMatrix(&centroid_B, centroid_B_matrix);
    int N = A->num;

    //将点云集A 质心平移至原点
    point_xyz_set *center_point_set_A;
    center_point_set_A = AllocatePoint_xyz_set(A->num);          
    Matrix4f tran_center_A;
    char eye[] = "eye";
    initia_matrix4f(tran_center_A, eye);
    tran_center_A[0][3] = -1 * centroid_A.x;
    tran_center_A[1][3] = -1 * centroid_A.y;
    tran_center_A[2][3] = -1 * centroid_A.z;
	
    //打印点集A的变换矩阵
    if (dubug) 
	{
		printf("-------------A的变换矩阵----------------\n");
	    PrintMatrix4f(tran_center_A);

     }
	//进行点集变换
	Matrix4fMultPoint_xyz_set(tran_center_A, A, center_point_set_A);
	//打印变换后的点集
	if (dubug)
	{
		printf("----------------A中心化  点集------------------\n");
		PrintPointXYZSet(center_point_set_A);
	}
	//将中心化后的点集A，放入到矩阵中
	point_set_matrix *A_point_set_Matrix;                       //已释放
	A_point_set_Matrix = Allocatpoint_set_matrix(center_point_set_A->num, 3);   //未释放
	PutPointSetIntoMatrix(center_point_set_A, A_point_set_Matrix);
	if (dubug)
	{
		printf("---------------A移至中心  矩阵----------------\n");
		PrintMatrix(A_point_set_Matrix->point_set);
	}
	//对B点集进行变换
    point_xyz_set *center_point_set_B;
    center_point_set_B = AllocatePoint_xyz_set(B->num);        
    Matrix4f tran_center_B;
    initia_matrix4f(tran_center_B, eye);
    tran_center_B[0][3] = -1 * centroid_B.x;
    tran_center_B[1][3] = -1 * centroid_B.y;
    tran_center_B[2][3] = -1 * centroid_B.z;

	//对点集B进行变换
    Matrix4fMultPoint_xyz_set(tran_center_B, B, center_point_set_B);
	
	//打印B点集
    if (dubug)
    {
		printf("-------------------将B移至中心，点集----------------\n");
	    PrintPointXYZSet(center_point_set_B);
     }
	//将点集B放入到矩阵中
	point_set_matrix *B_point_set_Matrix;                      //已释放
	B_point_set_Matrix = Allocatpoint_set_matrix(center_point_set_B->num, 3);
	PutPointSetIntoMatrix(center_point_set_B, B_point_set_Matrix);
	if (dubug)
	{
		printf("------------------B 移至中心 矩阵-------------------\n");
		PrintMatrix(B_point_set_Matrix->point_set);
	}


    //将A移至中心的表示矩阵，进行转置
    matrix *A_Transpose;                                       //已释放
    A_Transpose = AllocateMatrix(3, A_point_set_Matrix->point_set->r);
    Transpose(A_point_set_Matrix->point_set, A_Transpose);

    if (dubug)
    {
		printf("---------------A_Transpose，点集A的矩阵的转置------------\n");
	    PrintMatrix(A_Transpose);
     }

	matrix *H;
	H = AllocateMatrix(3, 3);
    MultMatrix(A_Transpose, B_point_set_Matrix->point_set, H);
    if (dubug) 
	{
		printf("----------------H矩阵-----------------\n" );
	    PrintMatrix(H);
     }
     double **u;
     double *w, **v;
     u = dmatrix(1, 3, 1, 3);
     w = dvector(1, 3);
     v = dmatrix(1, 3, 1, 3);
	 int i, j;
	 for (i = 1; i <= 3; i++)
	 {
		 for (j = 1; j <= 3; j++)
			 u[i][j] = H->entry[i - 1][j - 1];
	 }
	 svdcmp(u, 3, 3, w, v);
	 if (dubug)
	 {
		 //打印函数
		 printf("----------------------u矩阵-----------------\n");
		 for (i = 1; i <= 3; i++)
		 {
			 printf("%le %le %le\n", u[i][1], u[i][2], u[i][3]);
		 }
		 printf("----------------------v矩阵-----------------\n");
		 for (i = 1; i <= 3; i++)
		 {
			 printf("%le %le %le\n", v[i][1], v[i][2], v[i][3]);
		 }
	 }

	 matrix *U, *V;                              
	 U = AllocateMatrix(3, 3);
	 V = AllocateMatrix(3, 3);
	 for (int i = 0; i < 3; i++)
		 for (int j = 0; j < 3; j++)
		 {
			 U->entry[i][j] = u[i + 1][j + 1];
			 V->entry[i][j] = v[i + 1][j + 1];
		 }
	 matrix *UT;                                   //已释放
	 UT = AllocateMatrix(U->c, U->r);
	 Transpose(U, UT);
	 MultMatrix(V, UT, R);
	 if (Determinant33(R) < 0)
	 {
		 for (int i = 0; i < 3; i++)
		 {
			 V->entry[0][1] *= -1;
			 V->entry[1][1] *= -1;
			 V->entry[2][1] *= -1;
		 }
		 MultMatrix(V, UT, R);
	 }
	 if (dubug)
	 {
		 printf("-----------------R------------------\n");
		 PrintMatrix(R);
	 }

	 matrix *centroid_A_Tran;                      //已释放
	 centroid_A_Tran = AllocateMatrix(3, 1);
	 centroid_A_Tran->entry[0][0] = centroid_A.x;
	 centroid_A_Tran->entry[1][0] = centroid_A.y;
	 centroid_A_Tran->entry[2][0] = centroid_A.z;
	 matrix *tempRmultA_T;                          //已释放
	 tempRmultA_T = AllocateMatrix(3, 1);
	 MultMatrix(R, centroid_A_Tran, tempRmultA_T);

	 matrix *centroid_B_Tran;                        //已释放
	 centroid_B_Tran = AllocateMatrix(3, 1);
	 centroid_B_Tran->entry[0][0] = centroid_B.x;
	 centroid_B_Tran->entry[1][0] = centroid_B.y;
	 centroid_B_Tran->entry[2][0] = centroid_B.z;

	 SubMatrix(centroid_B_Tran, tempRmultA_T, trans);
	 if (true)
	 {
		 printf("---------------------T--------------------\n");
		 PrintMatrix(trans);
	 }
	 FreeMatrix(centroid_B_Tran);
	 FreeMatrix(tempRmultA_T);
	 FreeMatrix(UT);
	 FreeMatrix(U);
	 FreeMatrix(V);
	 free_dvector(w, 1, 3);
	 FreePoint_set_matrix(B_point_set_Matrix);
	 FreeMatrix(A_Transpose);
	 FreePoint_set_matrix(A_point_set_Matrix);
	 FreeMatrix(centroid_B_matrix);
	 FreeMatrix(centroid_A_matrix);
}