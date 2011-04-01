#include "itkImage.h"
#include <iostream>

#include <CIS_Array_Image3D.h>


int main()
{
	CIS_Array_Image3D_float *img3D = new CIS_Array_Image3D_float();
	img3D->Load("c:\imagedata\mr10_092_13p.i0344_53.hdr");

	system("pause");
	return 0;
}