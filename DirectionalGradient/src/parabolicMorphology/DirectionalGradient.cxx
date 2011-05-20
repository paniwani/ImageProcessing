#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkDirectionalGradientImageFilter.h>

int main()
{
	// load image
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("C:/ImageData/mr10_092_13p.i0344_100.hdr");
	WriteITK <ImageType2D> (input,"input.nii");

	// make mask
	/*ImageType2D::RegionType maskRegion;
	
	ImageType2D::IndexType maskIndex;
	maskIndex[0] = 256;
	maskIndex[1] = 256;
	
	ImageType2D::SizeType maskSize;
	maskSize[0] = 5;
	maskSize[1] = 5;

	maskRegion.SetIndex(maskIndex);
	maskRegion.SetSize(maskSize);*/

	ByteImageType2D::Pointer mask = ByteImageType2D::New();
	mask->SetRegions(input->GetLargestPossibleRegion());
	mask->SetSpacing(input->GetSpacing());
	mask->Allocate();
	mask->FillBuffer(0);
	
	ByteIteratorType2D maskIt(mask,input->GetLargestPossibleRegion());
	IteratorType2D inputIt(input,input->GetLargestPossibleRegion());

	for (inputIt.GoToBegin(), maskIt.GoToBegin(); !maskIt.IsAtEnd();++maskIt, ++inputIt)
	{
		if (inputIt.Get() > -200 && inputIt.Get() < 200)
			maskIt.Set(1);
	}

	WriteITK <ByteImageType2D> (mask,"mask.nii");
	
	// directional gradient
	typedef itk::DirectionalGradientImageFilter<ImageType2D,ByteImageType2D,ImageType2D> DirectionalGradientImageFilterType;
	DirectionalGradientImageFilterType::Pointer dgFilter = DirectionalGradientImageFilterType::New();
	dgFilter->SetInput(input);
	dgFilter->SetMaskImage(mask);
	dgFilter->SetOutsideValue(1);
	dgFilter->SetSigma(1);
	dgFilter->Update();
	WriteITK <ImageType2D> (dgFilter->GetOutput(),"dg.nii");




	system("pause");
	return 0;
}