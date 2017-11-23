#include "dataprocess.h"

int vector3fcopy(vector3f dest, const vector3f src)
{
	int vector_len = 3;
	int i;
	for (int i = 0; i < vector_len; i++)
	{
		dest[i] = src[i];
	}
	return 0;
}
