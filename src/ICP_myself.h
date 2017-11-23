#ifndef _ICP_MYSELF_H
#define _ICP_MYSELF_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "point3D.h"
#include "time_keeper.h"

/*
* linked list pointing members of a bin
*/
typedef struct member {
	struct member *next;    // pointer to next member
	int    thispt;          // array index of corresponding point
} member;

typedef struct bin {
	struct member *first;
	struct member *current;
} bin;
//定义多箱子结合体，用来盛放点云
typedef struct Bins {
	bin ***qbin;
	int numbineachdimen;
	//x_size, y_size, z_size, x_minboun, y_minboun, z_minboun
	double size_minboundary[6];
}Bins;

//定义函数
void PointCenterOfMass(point_xyz_set *q, point_xyz_set *q_centermass);
void FindClosePointViolent(point_xyz_set *p_base, point_xyz_set *q_new, int *closest_pt);
void FindClosePointViaElias(Bins *q, point_xyz_set *qst, point_xyz_set *p_new, int *closest_pt);
void ICPalgorithm_myself(matrix *Translation, matrix *fR, vector *ft, point_xyz_set *p_base, point_xyz_set *q);

//void ComputeRotationAndTranslation(void);
void ComputeRotationAndTranslation_SVD(point_xyz_set *p_base, point_xyz_set *q_new, int *closest_pt, matrix *fR, vector *ft);
void ApplyRotationAndTranslation(point_xyz_set *q_new, matrix *fR, vector *ft);

//Bins
void CreateBinsHashmap(Bins * q_bins, point_xyz_set *q);
bin ***Allocatebin(int num);
void Freebin(bin ***qbin, int num_bin);
void FreeBins(Bins* q_bins);
void FreePointXyzSet(point_xyz_set *pointSet);

#ifdef __cplusplus
}
#endif
#endif





