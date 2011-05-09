#include "itkImage.h"
#include <iostream>
#include <utils.h>

#include <itkScalarImageToTextureFeaturesFilter.h>
#include <itkNeighborhoodIterator.h>

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr9/mr2_002_13p.i0393/dcm",5,8);

	WriteITK <ImageType> (input, "input.nii");

	ImageType::SizeType radius;
	radius.Fill(1);

	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(radius,input,input->GetLargestPossibleRegion());

	ImageType::Pointer mask = ImageType::New();
	mask->SetRegions(input->GetLargestPossibleRegion());
	mask->SetSpacing(input->GetSpacing());
	mask->Allocate();

	ImageType::Pointer texture = ImageType::New();
	texture->SetRegions(input->GetLargestPossibleRegion());
	texture->SetSpacing(input->GetSpacing());
	texture->Allocate();
	texture->FillBuffer(0);

	IteratorType texture_iter(texture,input->GetLargestPossibleRegion());

	for (nit.GoToBegin(),texture_iter.GoToBegin(); !nit.IsAtEnd(); ++nit,++texture_iter)
	{
		mask->FillBuffer(0);

		for (int i=0; i<nit.Size()-1; i++)
		{
			bool isInBounds;
			PixelType val = nit.GetPixel(i,isInBounds);

			if (isInBounds)
				mask->SetPixel(nit.GetIndex(i),1);
		}	

		typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<ImageType> ScalarImageToTextureFeaturesFilterType;
		ScalarImageToTextureFeaturesFilterType::Pointer textureFilter = ScalarImageToTextureFeaturesFilterType::New();

		textureFilter->SetInput(input);
		textureFilter->SetMaskImage(mask);
		
		ScalarImageToTextureFeaturesFilterType::FeatureNameVectorPointer requestedFeatures = ScalarImageToTextureFeaturesFilterType::FeatureNameVector::New();
		requestedFeatures->push_back(ScalarImageToTextureFeaturesFilterType::TextureFeaturesFilterType::Energy);
		textureFilter->SetRequestedFeatures(requestedFeatures);
		
		textureFilter->FastCalculationsOn();
		textureFilter->Update();

		//WriteITK <ImageType> (textureFilter->GetOutput(0),"one.nii");

		const ScalarImageToTextureFeaturesFilterType::FeatureValueVectorDataObjectType * mean = textureFilter->GetFeatureMeansOutput();
		const ScalarImageToTextureFeaturesFilterType::FeatureValueVector * _mean = mean->Get();
		//std::cout << _mean->ElementAt(0) << std::endl;

		//std::cout << "stop" << std::endl;

		texture_iter.Set( _mean->ElementAt(0) );


	}

	WriteITK <ImageType> (texture,"texture.nii");

	system("pause");
	return 0;
}