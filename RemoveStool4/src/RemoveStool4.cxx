#include "itkImage.h"
#include <iostream>
#include <utils.h>

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/wr1_019_10s.i0462/dcm");

	ImageType::DirectionType direction = input->GetDirection();

	std::cout << direction[0][0] << " " << direction[1][1] << " " << direction[2][2] << std::endl;

	WriteITK <ImageType> (input, "input.nii");

	system("pause");
	return 0;
}