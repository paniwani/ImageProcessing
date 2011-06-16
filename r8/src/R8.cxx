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

	typedef itk::ConnectedThresholdImageFilter<ImageType,ByteImageType> ConnectedThresholdImageFilterType;
	ConnectedThresholdImageFilterType::Pointer connectedThresholdFilter = ConnectedThresholdImageFilterType::New();
	connectedThresholdFilter->SetInput(input);
	connectedThresholdFilter->SetReplaceValue(255);
	connectedThresholdFilter->SetLower(-250);
	connectedThresholdFilter->SetUpper(200);

	// Place seed in center of image
	ImageType::IndexType centerIndex;
	for (int i=0; i<Dimension; i++)
		centerIndex[i] = (long) size[i]/2;

	connectedThresholdFilter->SetSeed(centerIndex);
	connectedThresholdFilter->Update();

	ByteImageType::Pointer connectedThresholdImage = connectedThresholdFilter->GetOutput();
	WriteITK <ByteImageType> (connectedThresholdImage, "connectedThreshold.nii");

	// Remove components with gradient too high
	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	FloatImageType::Pointer gm = gmFilter->GetOutput();
	WriteITK <FloatImageType> (gm,"gm.nii");

	typedef itk::BinaryThresholdImageFilter<FloatImageType,ByteImageType> BinaryThresholdImageFilterType;
	BinaryThresholdImageFilterType::Pointer thresholder = BinaryThresholdImageFilterType::New();
	thresholder->SetInput(gm);
	thresholder->SetInsideValue(0);
	thresholder->SetOutsideValue(255);
	thresholder->SetLowerThreshold(300);
	thresholder->Update();
	ByteImageType::Pointer gmLow = thresholder->GetOutput();
	WriteITK <ByteImageType> (gmLow,"gmLow.nii");

	Mask <ByteImageType,ByteImageType> (connectedThresholdImage, gmLow);
	WriteITK <ByteImageType> (connectedThresholdImage,"ctLowGm.nii");

	// Compute mean and variance of connected tissue
	float mean=0, variance=0;
	unsigned int count=0;

	IteratorType it(input,region);
	ByteIteratorType mit(connectedThresholdImage,region);

	for (it.GoToBegin(), mit.GoToBegin(); !it.IsAtEnd(); ++it, ++mit)
	{
		if (mit.Get() != 0)
		{
			mean += it.Get();
			count++;
		}	
	}

	mean /= count;

	for (it.GoToBegin(), mit.GoToBegin(); !it.IsAtEnd(); ++it, ++mit)
	{
		if (mit.Get() != 0)
		{
			variance += (it.Get() - mean)*(it.Get() - mean);
		}	
	}

	variance /= count;

	std::cout << "Mean: " << mean << std::endl;
	std::cout << "Std: " << sqrt(variance) << std::endl;

	// Set input to air where gradient is too high
	Mask <ImageType, ByteImageType> (input,gmLow,-1024);
	WriteITK <ImageType> (input,"inputLowGm.nii");

	// Fuzzy with several standard deviations
	typedef itk::SimpleFuzzyConnectednessScalarImageFilter<ImageType,ByteImageType> FuzzySegmentationFilterType;
	typedef FuzzySegmentationFilterType::FuzzySceneType FuzzySceneType; 
	
	FuzzySegmentationFilterType::Pointer fuzzySegmenter = FuzzySegmentationFilterType::New();
	fuzzySegmenter->SetInput(input);
	fuzzySegmenter->SetObjectSeed( centerIndex );
	fuzzySegmenter->SetMean( 0 );
	fuzzySegmenter->SetThreshold( 0.5 );

	float stdArray[5] = {100,200,300,400,500};

	for (int i=0; i<5; i++)
	{
		float std = stdArray[i];

		fuzzySegmenter->SetVariance( std*std );
		
		fuzzySegmenter->Update();

		const FuzzySceneType * fuzzyScene = fuzzySegmenter->GetFuzzyScene();

		//typedef itk::MaskImageFilter<FuzzySceneType,ByteImageType> MaskerType;
		//MaskerType::Pointer masker = MaskerType::New();
		//masker->SetInput1(fuzzyScene);
		//masker->SetInput2(gmLow);

		// Rescale fuzzy scene
		typedef itk::RescaleIntensityImageFilter<FuzzySceneType,ByteImageType> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(fuzzyScene);
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();

		std::stringstream ss;
		ss << "fuzzyScene_std_" << stdArray[i] << ".nii";

		WriteITK <ByteImageType> (rescaler->GetOutput(),ss.str());
	}
	
	/*typedef itk::ImageFileWriter<FuzzySceneType> FuzzyWriterType;
	FuzzyWriterType::Pointer fwriter = FuzzyWriterType::New();
	fwriter->SetInput( fuzzySegmenter->GetFuzzyScene() );
	fwriter->SetFileName("fuzzyScene.nii");
	fwriter->Update();*/

	//// Compute centroid of tissue intensity
	//ImageType::IndexType centroidIndex;

	//float centroid[2] = {0,0};
	//float sumIntensity = 0;

	//IteratorType it(input,region);
	//for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	//{
	//	if ( it.Get() > -250 && it.Get() < 200 )
	//	{
	//		ImageType::IndexType idx = it.GetIndex();
	//		
	//		for (int i=0; i<2; i++)
	//		{
	//			centroid[i] += it.Get()*idx[i];
	//		}
	//		
	//		sumIntensity += it.Get();
	//	}
	//}	

	//for (int i=0; i<2; i++)
	//{
	//	centroid[i] /= sumIntensity;
	//	centroidIndex[i] = (long) centroid[i];
	//}

	//std::cout << "Centroid X: " << centroidIndex[0] << std::endl;
	//std::cout << "Centroid Y: " << centroidIndex[1] << std::endl;

	//// Allocate entire confidence image
	//ByteImageType::Pointer confidenceImage = ByteImageType::New();
	//confidenceImage->SetRegions(region);
	//confidenceImage->CopyInformation(input);
	//confidenceImage->Allocate();
	//confidenceImage->FillBuffer(0);

	//// Confidence connected on each colon segment
	//typedef unsigned short LabelType;
	//typedef itk::ShapeLabelObject<LabelType,Dimension> LabelObjectType;
	//typedef itk::LabelMap<LabelObjectType> LabelMapType;	
	//
	//typedef itk::BinaryImageToLabelMapFilter<ByteImageType,LabelMapType> BinaryImageToLabelMapFilterType;

	//BinaryImageToLabelMapFilterType::Pointer converter = BinaryImageToLabelMapFilterType::New();
	//converter->SetInput(colon);
	//converter->SetInputForegroundValue(255);
	//converter->SetOutputBackgroundValue(0);
	//converter->Update();
	//LabelMapType::Pointer colonMap = converter->GetOutput();

	//for (unsigned int label=1; label<=colonMap->GetNumberOfLabelObjects(); label++)
	//{
	//	// Get label object
	//	const LabelObjectType * labelObject = labelMap->GetLabelObject( label );
	//	
	//	// Create sub image containing this component only
	//	ImageType::RegionType labelRegion = labelObject->GetRegion();
	//	ImageType::Pointer inputSub = CropByRegion <ImageType> (input,labelRegion);	

	//	// Perform confidence connected on this region
	//	// Select first voxel with tissue-like intensity as seed		
	//	typedef itk::ConfidenceConnectedImageFilter<ImageType,ByteImageType> ConfidenceConnectedImageFilterType;
	//	ConfidenceConnectedImageFilterType::Pointer confidenceConnectedFilter = ConfidenceConnectedImageFilterType::New();
	//	confidenceConnectedFilter->SetInput(input);
	//	confidenceConnectedFilter->SetMultiplier(2.5);
	//	confidenceConnectedFilter->SetNumberOfIterations(2);
	//	confidenceConnectedFilter->SetReplaceValue(255);

	//	IteratorType sit(inputSub,inputSub->GetLargestPossibleRegion);
	//	for (sit.GoToBegin(); !sit.IsAtEnd(); ++sit)
	//	{
	//		if (sit.Get() > -250 && sit.Get() < 180)
	//		{
	//			confidenceConnectedFilter->AddSeed( sit.GetIndex() );
	//			break;
	//		}
	//	}
	//	
	//	



	//	confidenceConnectedFilter->Update();
	//	WriteITK <ByteImageType> (confidenceConnectedFilter->GetOutput(), "confidenceConnected.nii");
	//}







	
				
	system("pause");
	return 0; 									
} 												
