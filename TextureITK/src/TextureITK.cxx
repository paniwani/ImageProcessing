#include "itkImage.h"
#include <iostream>
#include <utils.h>

#include <itkScalarImageToTextureFeaturesFilter.h>
#include <itkNeighborhoodIterator.h>

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr9/mr2_002_13p.i0393/dcm",5,6);

	WriteITK <ImageType> (input, "input.nii");

	/*ImageType::SizeType radius;
	radius.Fill(1);

	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(radius,input,input->GetLargestPossibleRegion());

	ImageType::Pointer mask = ImageType::New();
	mask->SetRegions(input->GetLargestPossibleRegion());
	mask->SetSpacing(input->GetSpacing());
	mask->Allocate();

	for (nit.GoToBegin(); !nit.IsAtEnd(); ++nit)
	{*/
		/*mask->FillBuffer(0);

		for (int i=0; i<nit.Size()-1; i++)
		{
			mask->SetPixel(nit.GetIndex(i),1);
		}*/

		typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<ImageType> ScalarImageToTextureFeaturesFilterType;
		ScalarImageToTextureFeaturesFilterType::Pointer textureFilter = ScalarImageToTextureFeaturesFilterType::New();

		textureFilter->SetInput(input);
		//textureFilter->SetMaskImage(mask);
		
		ScalarImageToTextureFeaturesFilterType::FeatureNameVectorPointer requestedFeatures = ScalarImageToTextureFeaturesFilterType::FeatureNameVector::New();
		requestedFeatures->push_back(ScalarImageToTextureFeaturesFilterType::TextureFeaturesFilterType::Energy);
		textureFilter->SetRequestedFeatures(requestedFeatures);
		
		textureFilter->FastCalculationsOn();
		textureFilter->Update();

		//WriteITK <ImageType> (textureFilter->GetOutput(0),"one.nii");

		const ScalarImageToTextureFeaturesFilterType::FeatureValueVectorDataObjectType * mean = textureFilter->GetFeatureMeansOutput();

		std::cout << "stop" << std::endl;


	//}

	system("pause");
	return 0;
}