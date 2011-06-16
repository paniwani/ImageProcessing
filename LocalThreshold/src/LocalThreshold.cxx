const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>
 												
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

	// Compute global otsu threshold
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1200);
	otsuCalculator->Compute();
	PixelType globalOtsu = otsuCalculator->GetThreshold();

	// Allocate threshold image



	// Divide image into several subregions
	

	typedef itk::ImageRegionMultidimensionalSplitter<Dimension> SplitterType;
	SplitterType::Pointer splitter = SplitterType::New();
	
	unsigned int requestedNumDivisions = 40;
	unsigned int numDivisions = GetNumberOfSplits(region,requestedNumDivisions);

	for (unsigned int i=0; i < numDivisions; i++)
	{
		// Get local region
		ImageType::RegionType localRegion = splitter->GetSplit(i,numDivisions,region);
		
		// Use global threshold if no stool is present
		bool hasStool = false;
		
		IteratorType it(input,region);
		
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (it.Get() >= 200)
			{
				hasStool = true;
				break;
			}
		}
		
		if (hasStool)


		// Compute local threshold

	}


	return 0;
}