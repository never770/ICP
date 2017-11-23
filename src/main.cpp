#pragma warning(disable:4996)
#include <string.h>
#include<stdio.h>
#include"rdata_vis.h"
//#include "..\ICP-master\ICP.h"
#include "matrix.h"
#include "Datatype.h"
#include "Dataprocess.h"
#include "transfor_matrix2pos_pot.h"
#include "svdcmp.h"
#include ".\ICP_myself.h"
int Rentry2dcm(double **p, Matrix3f &dcm);
int main()
{	
	matrix      *R;
	vector      *t;
	double      Mf[16];
	R = AllocateMatrix(3, 3);
	t = AllocateVector(3);
	point_xyz_set *ps, *qs;

	// 读取基准模型
	int base_num;
	base_num = coutrow("Points_Moving.txt");
	ps = AllocatePoint_xyz_set(base_num);
	readtxt2pointset("Points_Moving.txt", ps, base_num);
	//printf("-------------Points_Moving.txt----------\n");
	//PrintPointXYZSet(ps);

	//读取变换后的点云
	int change_num;
	change_num = coutrow("Points_Static.txt");
	qs = AllocatePoint_xyz_set(change_num);
	readtxt2pointset("Points_Static.txt", qs, change_num);
	//printf("-------------Points_Static.txt----------\n");
	//PrintPointXYZSet(qs);

	matrix *Translation;
	Translation = AllocateMatrix(4, 4);
	char eye[] = "eye";
	Initia_matrix(Translation,eye);
	ICPalgorithm_myself(Translation, R, t, ps, qs);

	/*
	1. 解算出来的R复制到dcm中
	2.将dcm解算出欧拉角
	*/
	Matrix3f dcm;
	vector3f angle;
	Rentry2dcm(R->entry, dcm);
	dcm2angle(dcm, angle);
	printf("\n %f %f %f\n", angle[0], angle[1], angle[2]);

	FreeMatrix(R); FreeVector(t);
	FreePointXyzSet(ps);
	FreePointXyzSet(qs);

	return 0;
}
//int readpointxyz(char *filename, point_xyz *read, int &num)
//{
//	FILE *fp_txt;
//	if ((fp_txt = fopen(filename, "r")) == NULL) {
//		printf("Cannot read %s\n", filename);
//		return -1;
//	}
//	int i = 0;
//	float a, b, c;
//	while (((fscanf(fp_txt, "%f %f %f", &a, &b, &c))!= EOF) && i < num)
//	{
//		//&(read[i].x), &(read[i].y), &(read[i].z)
//		//printf("%f %f %f\n", a, b, c) ;
//		read[i].x = a;
//		read[i].y = b;
//		read[i].z = c;
//		i++;
//	}
//}
int Rentry2dcm(double **p, Matrix3f &dcm)
{
	int dcm_size = 3;
	int i = 0,
		j = 0;

	for (i = 0; i < dcm_size; i++)
		for(j = 0; j < dcm_size; j++)
	{
			dcm[i][j] = p[i][j];
	}
	return 0;
}