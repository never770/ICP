#pragma once
#ifndef DATATYPE_H
#define DATATYPE_H
#define PI 3.1415926
typedef  float Matrix3f[3][3];
typedef  float vector3f[3];
typedef  float Matrix4f[4][4];
//define the struct which contain the posture and position
typedef struct pos_pot_strcut{
	vector3f pos1;
	vector3f pot;
} pos_pot_strcut;
#endif