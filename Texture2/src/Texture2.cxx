#include "itkImage.h"
#include <iostream>
#include <utils.h>

#include <itkScalarImageToCooccurrenceMatrixFilter2.h>
#include <itkHistogramToTextureFeaturesFilter.h>

int main()
{
		
	// Load image
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("c:/imagedata/Modified_mr10_092_13p.i0344_85_100_40_slice13.hdr");
	WriteITK <ImageType2D> (input,"input.nii");
	
	ImageType2D::RegionType region = input->GetLargestPossibleRegion();
	IteratorType2D input_iter(input,region);

	const unsigned int dim = ImageType2D::ImageDimension;

	// Allocate texture image
	ImageType2D::Pointer textureImage = ImageType2D::New();
	textureImage->SetSpacing(input->GetSpacing());
	textureImage->SetRegions(region);
	textureImage->Allocate();
	textureImage->FillBuffer(0);
	IteratorType2D textureImage_iter(textureImage,region);

	// Set radius
	unsigned int radius = 1;

	ImageType2D::SizeType localSize;
	for (unsigned int i=0; i<dim; i++)
		localSize[i] = 2*radius + 1;

	// Compute texture locally
	int count=0;

	for (input_iter.GoToBegin(), textureImage_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++textureImage_iter)
	{
		// Get local region
		ImageType2D::IndexType idx = input_iter.GetIndex();
		
		//for (unsigned int i=0; i<dim; i++)
			//idx[i] -= radius;

		if (idx[0]+localSize[0] <= region.GetIndex()[0]+region.GetSize()[0]-1 && 
			idx[1]+localSize[1] <= region.GetIndex()[1]+region.GetSize()[1]-1 )
		{
			ImageType2D::RegionType localRegion;
			localRegion.SetIndex(idx);
			localRegion.SetSize(localSize);

			// Get histogram
			typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter2<ImageType2D> GLCMGeneratorType;
			GLCMGeneratorType::Pointer GLCMGenerator = GLCMGeneratorType::New();
			GLCMGenerator->SetLocalRegion(localRegion);
			GLCMGenerator->SetNormalize(true);

			ImageType2D::OffsetType offset;
			offset[0]=1;
			offset[1]=0;

			GLCMGenerator->SetOffset(offset);
			GLCMGenerator->SetInput(input);
			GLCMGenerator->SetPixelValueMinMax(-1500,1500);
			GLCMGenerator->Update();

			//std::cout << "Done GLCM: " << count << std::endl;

			//// Compute texture feature
			//typedef itk::Statistics::HistogramToTextureFeaturesFilter< GLCMGeneratorType::HistogramType > TexturesFilterType;
			//TexturesFilterType::Pointer textureFilter = TexturesFilterType::New();
			//
			//textureFilter->SetInput(GLCMGenerator->GetOutput());
			//textureFilter->Update();

			//textureImage_iter.Set(textureFilter->GetEnergy());

			//std::cout << "Done Texture: " << count << std::endl;
			//TexturesFilterType::MeasurementType measurement = textureFilter->GetEnergy();
		}
		
		count++;
	}

	//WriteITK <ImageType2D> (energy,"energy.nii");


	system("pause");
	return 0;
}