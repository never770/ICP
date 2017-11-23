//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/*
 * Some useful functions on 3D points
 * 
 * June 2001, Johnny Park
 */
#pragma warning(disable:4996)
#include "point3D.h"
#include "ICP_myself.h"

//为点云申请空间
point_xyz_set *AllocatePoint_xyz_set(int num)
{
	
	point_xyz_set *ps;
	ps = (point_xyz_set *)malloc(sizeof(point_xyz_set));
	ps->num = num;
	ps->head = (point_xyz *)malloc(num * sizeof(point_xyz));
	return ps;
}

void FreePointXyzSet(point_xyz_set *pointSet)
{
	free(pointSet->head);
	free(pointSet);
}
//释放点云空间
void FreePoint_xyz_set(point_xyz_set *ps)
{
	free(ps->head);
	free(ps);
}

/*
 * Eucleadian distance between two point_xyz's
 */
double E_distance(point_xyz *a, point_xyz *b)
{ 
  return(sqrt(SQ(a->x - b->x) + SQ(a->y - b->y) + SQ(a->z - b->z)));
}

/*
 * Angle between two unit normals
 */
double Angle(point_normal *a, point_normal *b)
{
  return( (180/PI) * acos((a->nx)*(b->nx)+(a->ny)*(b->ny)+(a->nz)*(b->nz)) );
}

/*
 * Finds center of mass of point_xyz's
 */
point_xyz CenterOfMass(point_xyz *a, int num_point)
{
  register int i;
  point_xyz center;

  center.x = 0; center.y = 0; center.z = 0;

  for (i=0; i<num_point; i++) {
    center.x += a[i].x;
    center.y += a[i].y;
    center.z += a[i].z;
  }
  center.x /= (double)num_point;
  center.y /= (double)num_point;
  center.z /= (double)num_point;

  return (center);
}



//test  github upload
//读取文件中的数组到点云集中
//添加一个函数，判断一个txt文件中保存了多少行
void readtxt2pointset(char *filename, point_xyz_set *set, int n)
{
	if (set->num != n)
	{
		printf("点集个数与读取点集个数不相等\n");
		exit(-1);
	}
	FILE *fp_txt;
	if ((fp_txt = fopen(filename, "r")) == NULL) {
		printf("Cannot read %s\n", filename);
		return -1;
	}
	int i = 0;
	float a, b, c;
	while (((fscanf(fp_txt, "%f %f %f", &a, &b, &c)) != EOF) && i < n)
	{
		//&(read[i].x), &(read[i].y), &(read[i].z)
		//printf("%f %f %f\n", a, b, c) ;
		set->head[i].x = a;
		set->head[i].y = b;
		set->head[i].z = c;
		i++;
	}
}

void PrintPointXYZSet(point_xyz_set * ps)
{
	int size = ps->num;
	for (int i = 0; i < size; i++)
	{
		printf(" %f %f %f\n", ps->head[i].x, ps->head[i].y, ps->head[i].z);
	}
}

void CopyPointXYZ(point_xyz * dst, point_xyz *src)
{
	dst->x = src->x;
	dst->y = src->y;
	dst->z = src->z;
}

int coutrow(char *filename)//输入工程目录下的文件名，或者其他目录下绝对路径名例如：c:\\1.txt;
{
	char c;
	int h = 0;
	FILE *fp;
	fp = fopen(filename, "r");
	if (fp == NULL)
		return -1;//表示文件打开错误
	while ((c = fgetc(fp)) != EOF)
	{
		if (c == '\n')
			h++;
		else
		{
			c = fgetc(fp);//这是处理最后一行可能没有换行标志，但是确文件结束。
			if (c == EOF)
			{
				h++;
				break;
			}
		}
	}
	return h;
}