#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <otbScalarImageToHigherOrderTexturesFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>

int main()
{
	// load input
	ImageType2D::Pointer input = ReadDicom <ImageType2D> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85);
	WriteITK <ImageType2D> (input,"input.nii");

	// crop
	ImageType2D::RegionType cRegion;
	
	ImageType2D::IndexType cIndex;
	cIndex[0] = 350; cIndex[1] = 230;

	ImageType2D::SizeType cSize;
	cSize[0] = 50; cSize[1] = 50;

	cRegion.SetIndex(cIndex);
	cRegion.SetSize(cSize);
	
	typedef itk::RegionOfInterestImageFilter<ImageType2D,ImageType2D> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(cRegion);
	cropper->Update();
	input = cropper->GetOutput();
	WriteITK <ImageType2D> (input,"inputCropped.nii");
	
	// rescale
	typedef itk::RescaleIntensityImageFilter<ImageType2D> RescaleIntensityImageFilterType;
	RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
	rescaler->SetInput(input);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	input = rescaler->GetOutput();
	WriteITK <ImageType2D> (input,"inputRescaled.nii");

	// compute textures
	typedef otb::ScalarImageToHigherOrderTexturesFilter<ImageType2D,ImageType2D> ScalarImageToHigherOrderTexturesFilterType;
	ScalarImageToHigherOrderTexturesFilterType::Pointer textureFilter = ScalarImageToHigherOrderTexturesFilterType::New();
	textureFilter->SetInput(input);
	textureFilter->SetInputImageMinimum(0);
	textureFilter->SetInputImageMaximum(255);

	ImageType2D::SizeType radius;
	radius.Fill(1);
	
	textureFilter->SetRadius(radius);

	ImageType2D::OffsetType offset;
	offset[0] = 0;
	offset[1] = 1;

	textureFilter->SetOffset(offset);

	textureFilter->Update();

	// get texture outputs
	for (int i=0; i<11; i++)
	{
		/*typedef itk::RescaleIntensityImageFilter<ImageType2D,ByteImageType2D> RescaleToByteType;
		RescaleToByteType::Pointer rbFilter = RescaleToByteType::New();
		rbFilter->SetInput(textureFilter->GetOutput(i));
		rbFilter->SetOutputMaximum(255);
		rbFilter->SetOutputMinimum(0);
		rbFilter->Update();*/

		std::stringstream ss;
		ss << "texture" << i << ".nii";
		WriteITK <ImageType2D> (textureFilter->GetOutput(i),ss.str());
	}

	return 0;
}