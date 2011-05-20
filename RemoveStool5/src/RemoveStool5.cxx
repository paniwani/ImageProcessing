#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include "itkColonSegmentationFilter.h"

int main()
{
	// load image
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10_092_13p.i0344/dcm",85,90);
	WriteITK <ImageType> (input,"input.nii");

	// segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonFilter = ColonSegmentationFilterType::New();
	colonFilter->SetInput(input);
	colonFilter->SetOutputForegroundValue(255);
	colonFilter->SetOutputBackgroundValue(0);
	colonFilter->SetRemoveBoneLung(false);
	colonFilter->Update();
	ByteImageType::Pointer colon = colonFilter->GetOutput();
	WriteITK <ByteImageType> (colon,"colon.nii");

	system("pause");
	return 0;
}