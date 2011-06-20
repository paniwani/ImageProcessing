const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkConfidenceConnectedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkNoiseImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBlackTopHatImageFilter.h>
#include <itkConvolutionImageFilter.h>
#include <otbScalarImageToTexturesMaskFilter.h>
 												
int main(int argc, char * argv[])				
{ 					
	// Load image
	ImageType::Pointer input = ReadITK <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm/mr10_092_13p_i0091.dcm");
	ImageType::RegionType region = input->GetLargestPossibleRegion();
		
	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ImageType> (input,colonRegion);
	Mask <ImageType,ByteImageType> (input,colon,-1025);

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	// Setup 8 outputs
	std::vector<FloatImageType::Pointer> outVector;
	outVector.resize(8);

	for (int i=0; i<8; i++)
	{
		outVector[i] = FloatImageType::New();
		outVector[i]->SetRegions(region);
		outVector[i]->CopyInformation(input);
		outVector[i]->Allocate();
		outVector[i]->FillBuffer(0);
	}

	// Setup offsets
	typedef ImageType::OffsetType OffsetType;

	std::vector<OffsetType> offsetVector;
	offsetVector.resize(2);
	offsetVector[0][0] = 0; offsetVector[0][1] = 1;
	offsetVector[1][0] = 1; offsetVector[1][1] = 0;

	// Setup texture filter
	typedef otb::ScalarImageToTexturesMaskFilter<FloatImageType,FloatImageType> TextureFilterType;
	TextureFilterType::Pointer textureFilter = TextureFilterType::New();
	textureFilter->SetInput( Rescale <ImageType,FloatImageType> (input,0,255) );
	textureFilter->SetInputImageMinimum(0);
	textureFilter->SetInputImageMaximum(255);
	textureFilter->SetNumberOfThreads(1);
	textureFilter->SetMaskImage( Cast <ByteImageType,FloatImageType> (colon) );

	ImageType::SizeType radius;
	radius.Fill(1);
	textureFilter->SetRadius(radius);

	for (int offcount = 0; offcount < 2; offcount++)
	{
		// Get textures for a particular offset
		textureFilter->SetOffset(offsetVector[offcount]);
		textureFilter->Update();
		
		// Sum textures across offsets
		for (int i=0; i<8; i++)
		{
			FloatImageType::Pointer texture = textureFilter->GetOutput(i);		
			outVector[i] = Add <FloatImageType> (outVector[i],texture);

		}
	}

	// Average offsets
	for (int i=0; i<8; i++)
	{
		DivideByConstant <FloatImageType> (outVector[i],offsetVector.size());

		std::stringstream ss;
		ss << "texture" << i << ".nii";
		WriteITK <ByteImageType> ( Rescale <FloatImageType,ByteImageType> (outVector[i],0,255) ,ss.str());
	}

	// Show gradient magnitude
	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeFilterType;
	GradientMagnitudeFilterType::Pointer gmFilter = GradientMagnitudeFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	WriteITK <FloatImageType> (gmFilter->GetOutput(),"gm.nii");

	// Show standard deviation
	typedef itk::NoiseImageFilter<ImageType,FloatImageType> NoiseFilterType;
	NoiseFilterType::Pointer noiseFilter = NoiseFilterType::New();
	noiseFilter->SetInput(input);
	noiseFilter->SetRadius(radius);
	noiseFilter->Update();
	WriteITK <FloatImageType> (noiseFilter->GetOutput(),"noise.nii");

	return 0;
}
