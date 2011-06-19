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
#include <itkSimpleFuzzyConnectednessScalarImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkBlackTopHatImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkFastBilateralImageFilter.h>
 												
int main(int argc, char * argv[])				
{ 					
	// Load image
	ImageType::Pointer input = ReadITK <ImageType> ("C:/ImageData/mr10_092_13p.i0344/dcm/mr10_092_13p_i0091.dcm");
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

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	// Smooth input
	typedef itk::FastBilateralImageFilter<ImageType,ImageType> BilateralFilterType;
	BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
	bilateralFilter->SetInput(input);
	bilateralFilter->SetRangeSigma(200);
	bilateralFilter->SetDomainSigma(5);
	bilateralFilter->Update();
	input = bilateralFilter->GetOutput();

	WriteITK <ImageType> (input,"inputSmoothed.nii");

	Mask <ImageType,ByteImageType> (input,colon,-1025);

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	// Get tagged mask
	typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> ThresholdFilterType;
	ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
	thresholder->SetInput(input);
	thresholder->SetInsideValue(255);
	thresholder->SetOutsideValue(0);
	thresholder->SetLowerThreshold(180);
	thresholder->Update();
	ByteImageType::Pointer tag = thresholder->GetOutput();
	WriteITK <ByteImageType> (tag,"tag.nii");

	// Find local minima at several radii
	unsigned int radius[4]  = {1,2,3,4};

	for (int i=0; i<4; i++)
	{
		typedef itk::BinaryBallStructuringElement<BytePixelType,Dimension> StructuringElementType;  
		StructuringElementType se;
		
		ByteImageType::SizeType rad;
		rad.Fill(0);
		rad[0] = radius[i];
		rad[1] = radius[i];

		se.SetRadius( rad );
		se.CreateStructuringElement();

		typedef itk::BlackTopHatImageFilter<ImageType,ImageType,StructuringElementType> BlackTopHatFilterType;
		BlackTopHatFilterType::Pointer blackTopHat = BlackTopHatFilterType::New();
		blackTopHat->SetInput(input);
		blackTopHat->SetKernel(se);
		blackTopHat->Update();

		ImageType::Pointer bth = blackTopHat->GetOutput();

		// Mask with tagged
		//Mask <ImageType,ByteImageType> (bth,tag);

		std::stringstream ss;
		ss << "bth_rad" << radius[i] << ".nii";
		WriteITK <ImageType> (bth,ss.str());
	}

	return 0;
}