#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkLocalRangeImageFilter.h>

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85,90);
	WriteITK <ImageType> (input, "input.nii");

	typedef itk::LocalRangeImageFilter<ImageType,ImageType> LocalRangeImageFilterType;
	LocalRangeImageFilterType::Pointer rangeFilter = LocalRangeImageFilterType::New();
	rangeFilter->SetInput(input);
	
	ImageType::SizeType radius;
	radius[0] = 1; radius[1] = 1; radius[2] = 0;

	rangeFilter->SetSize(radius);
	rangeFilter->Update();

	WriteITK <ImageType> (rangeFilter->GetOutput(),"range.nii");

	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(radius,input,input->GetLargestPossibleRegion());

	for (nit.GoToBegin(); !nit.IsAtEnd(); ++nit)
	{
		PixelType sum=0;

		for (unsigned int i=0; i<nit.Size(); i++)
		{
			sum += nit.GetPixel(i);
		}

		PixelType mean = sum / nit.Size();


	}


	system("pause");
	return 0;
}