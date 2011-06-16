const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkAdaptiveOtsuThresholdImageFilterModified.h>
 												
int main(int argc, char * argv[])				
{ 					
	// Load image
	ImageType::Pointer input = ReadITK <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm/mr10_092_13p_i0140.dcm");
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
	//Mask <ImageType,ByteImageType> (input,colon,-1025);

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	typedef itk::AdaptiveOtsuThresholdImageFilterModified<ImageType,ByteImageType> AdaptiveThresholdFilterType;
	AdaptiveThresholdFilterType::Pointer thresholder = AdaptiveThresholdFilterType::New();
	thresholder->SetInput(input);
	thresholder->SetInsideValue(255);

	/*
	thresholder->SetOutsideValue(0);
	thresholder->SetNumberOfHistogramBins(255);
	thresholder->SetNumberOfControlPoints(30);
	thresholder->SetNumberOfLevels(4);
	thresholder->SetNumberOfSamples(5000);
	*/

	thresholder->Update();

	WriteITK <ByteImageType> (thresholder->GetOutput(),"threshold.nii");
	
	return 0;
}