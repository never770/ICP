#ifndef TRANSFOR_MATRIX2POS_POT_H
#define TRANSFOR_MATRIX2POS_POT_H
#include "Datatype.h"

/*
 Calculates the set of rotation angles from direction cosine matrix
 The default rotation sequence is 'ZYX'.
 angle[0] is the rotation angle of z axis;
 angle[1] is the rotation angle of y axis;
 angle[2] is the rotation angle of x axis;
*/
int dcm2angle(const Matrix3f dcm, vector3f angle);
int angle2dcm(const vector3f angle, Matrix3f dcm);

//Calculates the posture and position from The transformation matrix
int trans_mat2pos_pot_struct(const Matrix4f &transf, pos_pot_strcut &pos_pot);
int pos_pot_struct2transmat(const pos_pot_strcut &pos_pot, Matrix4f &transf);
#endif