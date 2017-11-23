#include <stdio.h>
#include <math.h>
#include<string.h>
#include "Datatype.h"
#include "Dataprocess.h"
#include "transfor_matrix2pos_pot.h"


//根据旋转矩阵，计算欧拉角
int dcm2angle(const Matrix3f dcm, vector3f angle)
{
	//angle[0] is the rotation angles of z axis;
	//angle[1] is the rotation angles of y axis;
	//angle[2] is the rotation angles of x axis;
	if (dcm[2][0] != 1 && dcm[2][0] != -1)
	{
		double cos_theata, cos_theata1;
		angle[1] = -1 * asin(dcm[2][0]);
		cos_theata = cos(angle[1]);
		angle[2] = atan2(dcm[2][1] / cos_theata, dcm[2][2]/ cos_theata);
		angle[0] = atan2(dcm[1][0] / cos_theata, dcm[0][0] / cos_theata);
	}
	else
	{
		angle[0] = 0;
		if (dcm[2][0] == -1)
		{
			angle[1] = PI / 2.0;
			angle[2] = atan2(dcm[0][1], dcm[0][2]);
		}
		else
		{
			angle[1] = -1 * PI / 2.0;
			angle[2] = atan2(-1 * dcm[0][1], -1 * dcm[0][2]);
		}
		
	}
	int i;
	for (i = 0; i < 3; i++)
	{
		angle[i] = angle[i] * 180.0 / PI;
	}
	

	return 0;
}


//根据欧拉角，计算3*3旋转矩阵
int angle2dcm(const vector3f angle, Matrix3f dcm)
{
	//zyx
	/*
	 cy*cz    sy*sx*cz-sz*cx    sy*cx*cz+sz*sx
	 cy*sz    sy*sx*sz+cz*cx    sy*cx*sz-cz*sx
	 -sy      cy*sx             cy*cx
	*/
	//default angle[0]——x axis_psi, angle[1]——y axis_theta, angle[2]—z axis_phi  
	float cx = cos(angle[0]);
	float cy = cos(angle[1]);
	float cz = cos(angle[2]);
	float sx = sin(angle[0]);
	float sy = sin(angle[1]);
	float sz = sin(angle[2]);
	dcm[0][0] = cy*cz;
	dcm[0][1] = sy*sx*cz - sz*cx;
	dcm[0][2] = sy*cx*cz + sz*sx;

	dcm[1][0] = cy*sz;
	dcm[1][1] = sy*sx*sz + cz*cx;
	dcm[1][2] = sy*cx*sz - cz*sx;

	dcm[2][0] = -sy;
	dcm[2][1] = cy*sx;
	dcm[2][2] = cy*cx;
	return 0;
}

//根据变换矩阵，计算位置和姿态
int trans_mat2pos_pot_struct(const Matrix4f &transf, pos_pot_strcut &pos_pot)
{
	Matrix3f dcm;
	//copy rotation matrix from the transform matrix
	//memcpy(dcm[0], transf[0], 3);
	//memcpy(dcm[1], transf[1], 3);
	//memcpy(dcm[2], transf[2], 3);
	int i, j;
	for(i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			dcm[i][j] = transf[i][j];
		}
	dcm2angle(dcm, pos_pot.pos1);
	int vector_len = 3;
	for (i = 0; i < vector_len; i++)
	{
		pos_pot.pot[i] = transf[i][3];
	}
	return 0;
}

//根据位置和旋转角，计算4*4的变换矩阵
int pos_pot_struct2transmat(const pos_pot_strcut & pos_pot, Matrix4f & transf)
{
	//positon 
	transf[0][0] = pos_pot.pot[0];
	transf[0][1] = pos_pot.pot[1];
	transf[0][2] = pos_pot.pot[2];
	//calculate the rotation matrix
	Matrix3f dcm;
	angle2dcm(pos_pot.pos1, dcm);
	int i, j;
	int vector_len = 3;
	for(i = 0; i < vector_len; i++)
		for (int j = 0; j < vector_len; j++)
		{
			transf[i][j] = dcm[i][j];
		}
	return 0;
}
