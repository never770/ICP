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
 * 3D points
 * 
 * June 2001, Johnny Park
 */

#ifndef _POINT3D_H
#define _POINT3D_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#ifndef SQ
#define SQ(x)   ((x) * (x))
#endif

#ifndef PI
#define PI      3.1415927
#endif

/* 
 * Structures
 */
typedef struct point_xyz {        /* (x, y, z) position of a point*/
  double x, y, z;
} point_xyz;
typedef struct point_xyz_set {            /* set of point*/
	point_xyz *head;
	int num;
}point_xyz_set;
typedef struct point_normal {     /* unit normal of a point */ 
  double nx, ny, nz;
} point_normal;

typedef struct point_nor {     /* unit normal of a point */ 
  float nx, ny, nz;
} point_nor;

typedef struct point_rgb {        /* color RGB values of a point */
  double r, g, b;
} point_rgb;


typedef struct point3D {
  double x, y, z;      
  double nx, ny, nz;   
  double r, g, b;      
  double w;            
} point3D;


/*
 * Function Declaration
 */
point_xyz_set *AllocatePoint_xyz_set(int num);
double E_distance(point_xyz *a, point_xyz *b);
double Angle(point_normal *a, point_normal *b);
point_xyz CenterOfMass(point_xyz *a, int num_point);
int coutrow(char *filename);
void readtxt2pointset(char *filename, point_xyz_set *set, int n);
void PrintPointXYZSet(point_xyz_set *ps);
void CopyPointXYZ(point_xyz *dst, point_xyz *src);
#ifdef __cplusplus
}
#endif 
#endif
