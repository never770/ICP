#pragma warning(disable:4996)
#ifdef __cplusplus
extern "C" {
#endif
#include "svdcmp.h"
#include "ICP_myself.h"
#include "point3D.h"
#include<math.h>
#define num_p_t 5333
#define debug 0
double Dmax;
#define MAX   99999999;
//p_base  �����任�����뵽q��ȥ
void ICPalgorithm_myself(matrix *Translation, matrix *fR, vector *ft, point_xyz_set *p_base, point_xyz_set *q)
{
	//��ʼ������
	//������ʱ��������������������̵�p_base
	point_xyz_set *p_new;
	p_new = AllocatePoint_xyz_set(p_base->num);                                                   //δ�ͷ�      ���ͷ�
	
	double previous_error = 99999,    /* error in the previous iteration */
		delta_error = 99999;       /* current error - previous error */
	
	for (int i = 0; i < p_base->num; i++)
	{
		p_new->head[i].x = p_base->head[i].x;
		p_new->head[i].y = p_base->head[i].y;
		p_new->head[i].z = p_base->head[i].z;
	}

	//�������飬��������p_newÿ���㣬��q�е������
	int *closest_pt;
	closest_pt = (int *)malloc(p_new->num * sizeof(int));                                      //δ�ͷ��ڴ�             ���ͷ��ڴ�
	//��ʼ������
	for (int i = 0; i < p_new->num; i++)
	{
		closest_pt[i] = -1;
	}
	int iteration_count = 0;

	//�������Ļ���
	point_xyz_set *q_centermass;
	q_centermass = AllocatePoint_xyz_set(q->num);                                            //δ�ͷ��ڴ�    //���ͷ��ڴ�
	point_xyz_set *p_new_centermass;
	p_new_centermass = AllocatePoint_xyz_set(p_new->num);                                    //δ�ͷ��ڴ�     //���ͷ��ڴ�


	//���������㷨
	//����������ṹ������Hashmap
	Bins *q_bins;
	q_bins = (Bins *)malloc(sizeof(Bins));                                                 //δ�ͷ��ڴ�   //���ͷ�
	PointCenterOfMass(q, q_centermass);
	CreateBinsHashmap(q_bins, q_centermass);


	//meanDmax[0]  ��ƽ�����
	//meanDmax[1]   ��������Χ
	double meanDmax[2];
	meanDmax[0] = MAX;
	meanDmax[1] = MAX;

	/*����ICPֹͣ������
	1. ������������
	2.����ƽ�����С�ڹ涨���                          mean_error > low_bound_mean_error
	3. �������С�ڹ涨�Ĳ������                       delta_error > low_bound_delta_error
	*/

	double low_bound_mean_error = 0.0000000005;
	double low_bound_delta_error = 0.000000002;
	
	while (iteration_count < 100 /*&&
		(meanDmax[0] > low_bound_mean_error) &&
		(delta_error > low_bound_delta_error)*/
		)
	{
		printf("----------��%d�ε���---------------\n", iteration_count);

		PointCenterOfMass(p_new, p_new_centermass);
		FindClosePointViaElias(q_bins, q_centermass, p_new_centermass, closest_pt, meanDmax);

		//����µĿ��Ʊ���
		delta_error = fabs(meanDmax[0] - previous_error);
		previous_error = meanDmax[0];

		ComputeRotationAndTranslation_SVD(q, p_new, closest_pt, fR, ft);
		ApplyRotationAndTranslation(p_new, fR, ft);
		iteration_count++;
	}
	ComputeRotationAndTranslation_SVD(q, p_base, closest_pt, fR, ft);
	ApplyRotationAndTranslation(p_base, fR, ft);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Translation->entry[i][j] = fR->entry[i][j];
		}
	}
	Translation->entry[0][3] = ft->entry[0];
	Translation->entry[1][3] = ft->entry[1];
	Translation->entry[2][3] = ft->entry[2];

	//�ͷ��ڴ�
	FreeBins(q_bins);
	FreePointXyzSet(p_new_centermass);
	FreePointXyzSet(q_centermass);
	FreePointXyzSet(p_new);
	free(closest_pt);

}
/*
��p_new�е�ÿһ������q��Ѱ������㣬
������closest_pt��
*/
void FindClosePointViolent(point_xyz_set *q, point_xyz_set *p_new, int *closest_pt)
{
	double x_p_new, y_p_new, z_p_new;
	int minthispt;
	double distance, current_distance;
	for (int i = 0; i < p_new->num; i++)
	{
		minthispt = 0;
		distance = E_distance(&(p_new->head[i]), &(q->head[0]));
		for (int j = 1; j < q->num; j++)
		{
			current_distance = E_distance(&(p_new->head[i]), &(q->head[j]));
			if (distance > current_distance)
			{
				minthispt = j;
				distance = current_distance;
			}
		}
		closest_pt[i] = minthispt;
		printf("%d-%d\n", i, minthispt);
	}
}

void FindClosePointViaElias(Bins *q, point_xyz_set *qst, point_xyz_set *p_new, int *closest_pt, double *meanDmax)
{
	//x_size, y_size, z_size, x_minboun, y_minboun, z_minboun
	int num_p = p_new->num;
	int sub_sample = 1;
	int bin_i = -1, bin_j = -1, bin_k = -1;  /* index of bins */

	//printf("----------------Ѱ����������---------------\n");
	//printf("----------------q�㼯------------------------\n");
	//PrintPointXYZSet(qst);
	//printf("---------------p_new-------------------------\n");
	//PrintPointXYZSet(p_new);

	//int *closet_pt_temp;
	//closet_pt_temp = (int *)calloc(p_new->num, sizeof(int));
	//FindClosePointViolent(qst, qst, closet_pt_temp);

	double binsize_x = q->size_minboundary[0];
	double binsize_y = q->size_minboundary[1];
	double binsize_z = q->size_minboundary[2];
	double min_halfbin_x = q->size_minboundary[3];
	double min_halfbin_y = q->size_minboundary[4];
	double min_halfbin_z = q->size_minboundary[5];
	double max_halfbin_x = min_halfbin_x + q->numbineachdimen * binsize_x;
	double max_halfbin_y = min_halfbin_y + q->numbineachdimen * binsize_y;
	double max_halfbin_z = min_halfbin_z + q->numbineachdimen * binsize_z;
	int bin_num_erverdimen = q->numbineachdimen;

	double init_min_distance = sqrt(SQ(max_halfbin_x - max_halfbin_x) +
		SQ(max_halfbin_y - min_halfbin_y) +
		SQ(max_halfbin_z - min_halfbin_z));

	double current_distance;
	int min_distance_pt;
	double min_distance;
	double mean_error = 0;
	int current_q = -1;
	int start_i, end_i,       /* start and end of neighbor index */
		start_j, end_j,
		start_k, end_k;
	int num_closest_pt = 0;
	double *dist_pq;
	dist_pq = (double *)malloc(p_new->num * sizeof(double));
	Dmax = meanDmax[1];

	for (int m = 0; m < num_p; m++)
	{
		/*printf("m   %f  %f  %f\n", qst->head[m].x, qst->head[m].y, qst->head[m].z);*/
		//x_size, y_size, z_size, x_boun, y_boun, z_boun
		//bin_i = (int)((qst->head[m].x - min_halfbin_x) / binsize_x);
		//bin_j = (int)((qst->head[m].y - min_halfbin_y) / binsize_y);
		//bin_k = (int)((qst->head[m].z - min_halfbin_z) / binsize_z);


		bin_i = (int)((p_new->head[m].x - q->size_minboundary[3]) / q->size_minboundary[0]);
		bin_j = (int)((p_new->head[m].y - q->size_minboundary[4]) / q->size_minboundary[1]);
		bin_k = (int)((p_new->head[m].z - q->size_minboundary[5]) / q->size_minboundary[2]);

		bin_i = (bin_i < 0) ? 0 : bin_i;
		bin_i = (bin_i > (q->numbineachdimen)) ? q->numbineachdimen: bin_i;

		bin_j = (bin_j < 0) ? 0 : bin_j;
		bin_j = (bin_j > (q->numbineachdimen)) ? q->numbineachdimen : bin_j;

		bin_k = (bin_k < 0) ? 0 : bin_k;
		bin_k = (bin_k > (q->numbineachdimen)) ? q->numbineachdimen : bin_k;

		{
			min_distance = init_min_distance;
			min_distance_pt = -1;
			//��������    �Ź������������  
			int search_radius = 1;

			//��֤������Χ���ڵ�����
			while (search_radius < q->numbineachdimen && 
				   search_radius * binsize_x * 0.2 < Dmax  &&
				   search_radius * binsize_y * 0.2 < Dmax  &&
				   search_radius * binsize_z * 0.2 < Dmax)
			{

				start_i = (bin_i - search_radius < 0) ? 0 : bin_i - search_radius;
				end_i = (bin_i + search_radius > bin_num_erverdimen - 1) ?
					bin_num_erverdimen - 1 : bin_i + search_radius;

				start_j = (bin_j - search_radius < 0) ? 0 : bin_j - search_radius;
				end_j = (bin_j + search_radius > bin_num_erverdimen - 1) ?
					bin_num_erverdimen - 1 : bin_j + search_radius;

				start_k = (bin_k - search_radius < 0) ? 0 : bin_k - search_radius;
				end_k = (bin_k + search_radius > bin_num_erverdimen - 1) ?
					bin_num_erverdimen - 1 : bin_k + search_radius;

				for (int i = start_i; i <= end_i; i++)
					for (int j = start_j; j <= end_j; j++)
						for (int k = start_k; k <= end_k; k++) {
							//�����һ������������ʱ������
							if (i > bin_i - (search_radius - 1) && i < bin_i + (search_radius - 1) &&
								j > bin_j - (search_radius - 1) && j < bin_j - (search_radius - 1) &&
								k > bin_k - (search_radius - 1) && k < bin_k - (search_radius - 1))
							{
								printf("����һ��\n");
								continue;

							}
							//���i, j, k���ڵ����ھ�����

							/* at least one member is in this bin */
							if (q->qbin[i][j][k].first != NULL) {
								/* set a pointer to the first member of this bin */
								q->qbin[i][j][k].current = q->qbin[i][j][k].first;

								/* find the member with the closest E-distance */
								while (q->qbin[i][j][k].current != NULL) {
									current_q = q->qbin[i][j][k].current->thispt;

									current_distance = E_distance(&(p_new->head[m]), &(qst->head[current_q]));
 									if (current_distance < min_distance) {
										min_distance = current_distance;
										min_distance_pt = current_q;
									}
									q->qbin[i][j][k].current = q->qbin[i][j][k].current->next;
									//printf("\tcurrent_distance:%f\n", current_distance);
									//printf("\tmin_distance:%f\n", min_distance);
								}
							}
						}

				//����ҵ�����㣬��ֹͣ��������������������Χ
				if (min_distance_pt != -1)
				{
					closest_pt[m] = min_distance_pt;
					mean_error += min_distance;
					dist_pq[m] = min_distance;
					num_closest_pt++;
					//printf("%d-----%d \n ", m, min_distance_pt);
					break;
				}
				else
				{
					search_radius++;

				}
			}


			printf("%d--hashmap--%d\n", m, closest_pt[m]);

		}

	}

	mean_error /= (double)num_closest_pt;
	double std_error = 0;
	for (int i = 0; i < p_new->num; i++) 
	{
		if(closest_pt[i] != -1)
		{
		    std_error += (dist_pq[i] - mean_error) * (dist_pq[i] - mean_error);
		}
	}
	std_error /= (double)num_closest_pt;
	std_error = sqrt(std_error);
	Dmax = 3.0*std_error + mean_error;
	meanDmax[0] = mean_error;
	meanDmax[1] = Dmax;

}

//��p_new���б任��ʹ֮��q�غ�
void ComputeRotationAndTranslation_SVD(point_xyz_set *q, point_xyz_set *p_new, int *closest_pt, matrix *fR, vector *ft)
{
	//��������������������Լ�
	matrix *trans_tmp;
	trans_tmp = AllocateMatrix(3, 1);                                                         //δ�ͷ��ڴ�         //���ͷ�
	point_xyz_set *p_new_close_pt, *q_close_pt;

	//ͳ���ҵ��˼��������
	int close_pt_num = 0;
	for (int i = 0; i < p_new->num; i++)
	{
		if (closest_pt[i] != -1)
		{
			close_pt_num++;
		}
	}
	p_new_close_pt = AllocatePoint_xyz_set(close_pt_num);                                       //δ�ͷ��ڴ�         //���ͷ�
	q_close_pt = AllocatePoint_xyz_set(close_pt_num);                                            //δ�ͷ��ڴ�        //���ͷ�
	int p_new_close_pt_current = 0;
	for (int i = 0; i < p_new->num; i++)
	{
		if (closest_pt[i] != -1)
		{
			CopyPointXYZ(&(p_new_close_pt->head[p_new_close_pt_current]), &(p_new->head[i]));
			CopyPointXYZ(&(q_close_pt->head[p_new_close_pt_current]), &(q->head[closest_pt[i]]));
			p_new_close_pt_current++;
		}
	}
	//��p_new_close_pt���б任��ʹ֮��q_close_pt�غ�


	//printf("--------------p_new_close_pt-------------\n");
	//PrintPointXYZSet(p_new_close_pt);
	//printf("--------------q_new_close_pt--------------\n");
	//PrintPointXYZSet(q_close_pt);

	rigid_transform3D(p_new_close_pt, q_close_pt, fR, trans_tmp);
	printf("-------------------p_new_close_pt----------------\n");
	PrintPointXYZSet(p_new_close_pt);
	printf("-------------------q_close_pt----------------\n");
	PrintPointXYZSet(q_close_pt);
	for (int j = 0; j < 3; j++)
	{
		ft->entry[0] = trans_tmp->entry[0][0];
		ft->entry[1] = trans_tmp->entry[1][0];
		ft->entry[2] = trans_tmp->entry[2][0];
	}

	FreeMatrix(trans_tmp);
	FreePointXyzSet(p_new_close_pt);
	FreePointXyzSet(q_close_pt);
}
void ApplyRotationAndTranslation(point_xyz_set *p_new, matrix *fR, vector *ft)
{
	if (true)
	{
		printf("---------��ת����-----------\n");
		PrintMatrix(fR);
		printf("------------ƽ������--------\n");
		PrintVector(ft);
	}
	double x_new, y_new, z_new;
	for (int i = 0; i < p_new->num; i++) {
		x_new = (fR->entry[0][0] * p_new->head[i].x) +
			(fR->entry[0][1] * p_new->head[i].y) +
			(fR->entry[0][2] * p_new->head[i].z) + ft->entry[0];

		y_new = (fR->entry[1][0] * p_new->head[i].x) +
			(fR->entry[1][1] * p_new->head[i].y) +
			(fR->entry[1][2] * p_new->head[i].z) + ft->entry[1];

		z_new = (fR->entry[2][0] * p_new->head[i].x) +
			(fR->entry[2][1] * p_new->head[i].y) +
			(fR->entry[2][2] * p_new->head[i].z) + ft->entry[2];
		p_new->head[i].x = x_new;
		p_new->head[i].y = y_new;
		p_new->head[i].z = z_new;
	}
}
void SimultaneousMaxMin(double *maxminxyz, const point_xyz_set *q_set)
{
	int num_q = q_set->num;
	point_xyz *q = q_set->head;
	int i;
	//printf("finding max and min of (x,y,z) values...\n"); fflush(stdout);
	/* initialize values */
	double max_x = -9999999, max_y = -9999999, max_z = -9999999, min_x = 9999999, min_y = 9999999, min_z = 9999999;
	for (i = 0; i < floor(num_q / 2); i++) {
		//��x����ֵ
		if (q[2 * i].x > q[2 * i + 1].x) {
			if (max_x < q[2 * i].x)    max_x = q[2 * i].x;
			if (min_x > q[2 * i + 1].x)  min_x = q[2 * i + 1].x;
		}
		else {
			if (max_x < q[2 * i + 1].x)  max_x = q[2 * i + 1].x;
			if (min_x > q[2 * i].x)    min_x = q[2 * i].x;
		}
		//��y����ֵ
		if (q[2 * i].y > q[2 * i + 1].y) {
			if (max_y < q[2 * i].y)    max_y = q[2 * i].y;
			if (min_y > q[2 * i + 1].y)  min_y = q[2 * i + 1].y;
		}
		else {
			if (max_y < q[2 * i + 1].y)  max_y = q[2 * i + 1].y;
			if (min_y > q[2 * i].y)    min_y = q[2 * i].y;
		}
		//��z����ֵ                                ��ֵȡֵ�����
		if (q[2 * i].z > q[2 * i + 1].z) {
			if (max_z < q[2 * i].z)    max_z = q[2 * i].z;
			if (min_z > q[2 * i + 1].z)  min_z = q[2 * i + 1].z;
		}
		else {
			if (max_z < q[2 * i + 1].z)  max_z = q[2 * i + 1].z;
			if (min_z > q[2 * i].z)    min_z = q[2 * i].z;
		}
	}

	/* check the last element in case num_p is an odd number*/
	if (max_x < q[num_q - 1].x)  max_x = q[num_q - 1].x;
	if (min_x > q[num_q - 1].x)  min_x = q[num_q - 1].x;
	if (max_y < q[num_q - 1].y)  max_y = q[num_q - 1].y;
	if (min_y > q[num_q - 1].y)  min_y = q[num_q - 1].y;
	if (max_z < q[num_q - 1].z)  max_z = q[num_q - 1].z;
	if (min_z > q[num_q - 1].z)  min_z = q[num_q - 1].z;

	printf("\tx: (%f, %f)\n\ty: (%f, %f)\n\tz: (%f, %f)\n\n",
		min_x, max_x, min_y, max_y, min_z, max_z);
	maxminxyz[0] = max_x; maxminxyz[1] = min_x;
	maxminxyz[2] = max_y; maxminxyz[3] = min_y;
	maxminxyz[4] = max_z; maxminxyz[5] = min_z;
	PrintTime();  fflush(stdout);
}
void CreateBinsHashmap(Bins *q_bins, point_xyz_set *q)
{
	int num_q = q->num;
	int num_bin = (int)(pow(q->num, 1.0 / 3)) + 1;
	q_bins->numbineachdimen = num_bin;

	//�������鲢�ҳ����Ƽ��������Сֵ
	//����0��Ԫ��  x�����ֵ�� 1��Ԫ�� x����Сֵ
	double maxminxyz[6];
	SimultaneousMaxMin(maxminxyz, q);
	//��������
	//Ӧ���ǰѵ���뵽������
	int i, j;
	int bin_i, bin_j, bin_k;      /* index of bins */

									//x_size, y_size, z_size, x_boun, y_boun, z_boun
	q_bins->size_minboundary[0] = (double)(maxminxyz[0] - maxminxyz[1] + 1.0) / (double)(q_bins->numbineachdimen);
	q_bins->size_minboundary[1] = (double)(maxminxyz[2] - maxminxyz[3] + 1.0) / (double)(q_bins->numbineachdimen);
	q_bins->size_minboundary[2] = (double)(maxminxyz[4] - maxminxyz[5] + 1.0) / (double)(q_bins->numbineachdimen);

	q_bins->size_minboundary[3] = maxminxyz[1];// - (q_bins->size_minboundary[0] / 2.0);
	q_bins->size_minboundary[4] = maxminxyz[3]; //- (q_bins->size_minboundary[1] / 2.0);
	q_bins->size_minboundary[5] = maxminxyz[5]; //- (q_bins->size_minboundary[2] / 2.0);


	q_bins->qbin = Allocatebin(q_bins->numbineachdimen);
	


	for (i = 0; i < num_q; i++) {
		//����õ�λ�ڵ��Ǹ�������
		//x_size, y_size, z_size, x_boun, y_boun, z_boun
		printf("i %f, %f, %f \t", q->head[i].x, q->head[i].y, q->head[i].z);
		bin_i = (int)((q->head[i].x - q_bins->size_minboundary[3]) / q_bins->size_minboundary[0]);
		bin_j = (int)((q->head[i].y - q_bins->size_minboundary[4]) / q_bins->size_minboundary[1]);
		bin_k = (int)((q->head[i].z - q_bins->size_minboundary[5]) / q_bins->size_minboundary[2]);
		printf("%d %d %d\n", bin_i, bin_j, bin_k);
		/*һ�����ӵ�������������1.��һ���ڵ�  2. ��ǰ�ڵ�*/
		//��ָ��ĸ���
		/* first member */
		if (q_bins->qbin[bin_i][bin_j][bin_k].first == NULL) {
			printf("�ýڵ��һ�ν���\n");
			q_bins->qbin[bin_i][bin_j][bin_k].first = (member *)calloc(1, sizeof(member));
			q_bins->qbin[bin_i][bin_j][bin_k].current = q_bins->qbin[bin_i][bin_j][bin_k].first;
		}
		/* not a first member */
		else {
			printf("�ýڵ㲻�ǵ�һ�ν���\n");
			q_bins->qbin[bin_i][bin_j][bin_k].current->next =
				(member *)calloc(1, sizeof(member));
			q_bins->qbin[bin_i][bin_j][bin_k].current =
				q_bins->qbin[bin_i][bin_j][bin_k].current->next;
		}
		q_bins->qbin[bin_i][bin_j][bin_k].current->thispt = i;
	}
}

//Ϊ���������ڴ�
bin ***Allocatebin(int num_bin)
{
	bin ***qbin;
	qbin = (bin ***)malloc(sizeof(bin **) * num_bin);
	for (int i = 0; i < num_bin; i++) {
		qbin[i] = (bin **)malloc(sizeof(bin *) * num_bin);
		for (int j = 0; j < num_bin; j++) {
			qbin[i][j] = (bin *)calloc(num_bin, sizeof(bin));
		}
	}
	return qbin;
}
void Freebin(bin ***qbin, int num_bin)
{
	for (int i = 0; i < num_bin; i++)
	{
		for (int j = 0; j < num_bin; j++)
		{
			free(qbin[i][j]);
		}
		free(qbin[i]);

	}
}

void FreeBins(Bins* q_bins)
{
	Freebin(q_bins->qbin,q_bins->numbineachdimen);
	free(q_bins);
}

//�ͷ��ڴ�

void PointCenterOfMass(point_xyz_set *q, point_xyz_set *q_centermass)
{
	point_xyz centroid_q;
	centroid_q = CenterOfMass(q->head, q->num);
	int N = q->num;
	//�����Ƽ�q ����ƽ����ԭ��
	Matrix4f tran_center_q;
	char eye[] = "eye";
	initia_matrix4f(tran_center_q, eye);
	tran_center_q[0][3] = -1 * centroid_q.x;
	tran_center_q[1][3] = -1 * centroid_q.y;
	tran_center_q[2][3] = -1 * centroid_q.z;
	//���е㼯�任
	Matrix4fMultPoint_xyz_set(tran_center_q, q, q_centermass);
}

#ifdef __cplusplus
}
#endif

