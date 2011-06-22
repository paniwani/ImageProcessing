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
#include <otbScalarImageToTexturesFilter.h>
#include <otbScalarImageToTexturesMaskFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkRelabelLabelMapFilter.h>
 												
int main(int argc, char * argv[])				
{ 					
	if( argc < 3 )
    {
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImage radius numBins allOffsets" << std::endl;
		system("pause");
		return 1;
    }

	// Load image
	ImageType::Pointer inputOriginal = ReadITK <ImageType> (argv[1]);
	
	// Convert to unsigned char [0-255]
	ByteImageType::Pointer input = Rescale <ImageType,ByteImageType> (inputOriginal,0,255);
	
	ByteImageType::RegionType region = input->GetLargestPossibleRegion();
	WriteITK <ByteImageType> (input,"input.nii");

	unsigned int Radius = atoi( argv[2] );
	unsigned int numBins = atoi( argv[3] );
	unsigned int allOffsets = atoi( argv[4] );

	std::stringstream ss;
	
	/*
	// Segment colon
	typedef itk::ColonSegmentationFilter<ByteImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ByteImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ByteImageType> (input,colonRegion);
	Mask <ImageType,ByteImageType> (input,colon,-1025);
	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ByteImageType::SizeType size = region.GetSize();
	*/

	//// Get gradient magnitude
	//typedef itk::GradientMagnitudeImageFilter<ByteImageType,FloatImageType> GradientMagnitudeFilterType;
	//GradientMagnitudeFilterType::Pointer gmFilter = GradientMagnitudeFilterType::New();
	//gmFilter->SetInput(input);
	//gmFilter->Update();
	//WriteITK <FloatImageType> (gmFilter->GetOutput(),"gm.nii");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Haralick Textures
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Setup 8 outputs
	std::vector<DoubleImageType::Pointer> outVector;
	outVector.resize(8);

	for (int i=0; i<8; i++)
	{
		outVector[i] = DoubleImageType::New();
		outVector[i]->SetRegions(region);
		outVector[i]->CopyInformation(input);
		outVector[i]->Allocate();
		outVector[i]->FillBuffer(0);
	}

	// Setup offsets
	typedef DoubleImageType::OffsetType OffsetType;

	std::vector<OffsetType> offsetVector;

	if (allOffsets)
	{
		offsetVector.resize(4);
		offsetVector[0][0] = -1; offsetVector[0][1] = 0;
		offsetVector[1][0] = -1; offsetVector[1][1] = -1;
		offsetVector[0][0] = 0; offsetVector[0][1] = -1;
		offsetVector[1][0] = 1; offsetVector[1][1] = -1;
	} else {
		offsetVector.resize(2);
		offsetVector[0][0] = 0; offsetVector[0][1] = 1;
		offsetVector[1][0] = 1; offsetVector[1][1] = 0;
	}

	// Setup texture filter
	typedef otb::ScalarImageToTexturesFilter<DoubleImageType,DoubleImageType> TextureFilterType;
	TextureFilterType::Pointer textureFilter = TextureFilterType::New();
	textureFilter->SetInput( Cast <ByteImageType, DoubleImageType> (input) ); // use double type input
	textureFilter->SetInputImageMinimum(0);
	textureFilter->SetInputImageMaximum(255);
	textureFilter->SetNumberOfThreads(1);
	textureFilter->SetNumberOfBinsPerAxis( numBins );
	//textureFilter->SetMaskImage( Cast <ByteImageType,FloatImageType> (colon) );

	DoubleImageType::SizeType radius;
	radius.Fill(Radius);
	textureFilter->SetRadius(radius);

	for (int offcount = 0; offcount < offsetVector.size(); offcount++)
	{
		// Get textures for a particular offset
		textureFilter->SetOffset(offsetVector[offcount]);
		textureFilter->Update();
		
		// Sum textures across offsets
		for (int i=0; i<8; i++)
		{
			DoubleImageType::Pointer texture = textureFilter->GetOutput(i);		
			outVector[i] = Add <DoubleImageType> (outVector[i],texture);

			std::stringstream ss;
			ss << "texture" << i << "_offset_" << offcount << ".nii";
			//WriteITK <FloatImageType> (texture,ss.str());

		}
	}

	// Average offsets
	for (int i=0; i<8; i++)
	{
		DivideByConstant <DoubleImageType> (outVector[i],offsetVector.size());

		std::stringstream ss;
		ss << "textureAveragedOffsets" << i << "_radius_" << Radius << "_numBins_" << numBins << "_allOffsets_" << allOffsets << ".nii";
		WriteITK <ByteImageType> ( Rescale <DoubleImageType,ByteImageType> (outVector[i],0,255) ,ss.str());
	}

	return 0;
}	