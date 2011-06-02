#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkMorphologicalDistanceTransformImageFilter.h>

int main()
{
	// make image
	ImageType2D::Pointer input = ImageType2D::New();

	ImageType2D::RegionType region;
	
	ImageType2D::IndexType index;
	index[0] = 0;
	index[1] = 0;

	ImageType2D::SizeType size;
	size[0] = 5;
	size[1] = 5;

	region.SetIndex(index);
	region.SetSize(size);

	input->SetRegions(region);

	ImageType2D::SpacingType spacing;
	spacing[0] = 0.5;
	spacing[1] = 0.5;

	input->SetSpacing(spacing);

	input->Allocate();
	input->FillBuffer(0);
	
	ImageType2D::IndexType idx;
	idx[0] = 2;
	idx[1] = 2;

	input->SetPixel(idx,255);

	WriteITK <ImageType2D> (input,"input.nii");

	// get distance map
	typedef itk::MorphologicalDistanceTransformImageFilter<ImageType2D,ImageType2D> DistanceType;
	DistanceType::Pointer distanceFilter = DistanceType::New();
	distanceFilter->SetInput(input);
	distanceFilter->SetOutsideValue(255);
	distanceFilter->Update();

	WriteITK <ImageType2D> (distanceFilter->GetOutput(),"distance.nii");

	system("pause");
	return 0;
}