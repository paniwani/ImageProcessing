const unsigned int Dimension = 3;
#include <itkImage.h> 						
#include <iostream> 							
#include <utils2.h> 
#include <itkColonSegmentationFilter.h>
#include <itkMaskImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkOrImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkFastBilateralImageFilter.h>

int main(int argc, char * argv[])				
{ 		
	// Load image
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",140,160);
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Segment colon
	typedef itk::ColonSegmentationFilter< ImageType, ByteImageType > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( input );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();
	
	ImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ImageType> (input,colonRegion);

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	// Smooth input
	typedef itk::FastBilateralImageFilter<ImageType,ImageType> BilateralImageFilterType;
	BilateralImageFilterType::Pointer bilateralFilter = BilateralImageFilterType::New();
	bilateralFilter->SetInput(input);
	bilateralFilter->SetDomainSigma(3);
	bilateralFilter->SetRangeSigma(50);
	bilateralFilter->Update();
	input = bilateralFilter->GetOutput();
	WriteITK <ImageType> (input,"inputSmoothed.nii");

	// Mask input with colon
	Mask <ImageType,ByteImageType> (input,colon,-1025);

	//// Find all air
	//typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> BinarizerType;
	//BinarizerType::Pointer binarizer = BinarizerType::New();
	//binarizer->SetInput(input);
	//binarizer->SetOutsideValue(0);
	//binarizer->SetInsideValue(255);
	//binarizer->SetUpperThreshold(-600);
	//binarizer->Update();
	//ByteImageType::Pointer air = binarizer->GetOutput();
	//WriteITK <ByteImageType> (air,"air.nii");

	//// Subtract background
	//typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType> BinaryKeeperType;
	//BinaryKeeperType::Pointer keeper = BinaryKeeperType::New();
	//keeper->SetInput(air);
	//keeper->SetAttribute("Size");
	//keeper->SetNumberOfObjects(1);
	//keeper->SetForegroundValue(255);
	//keeper->SetBackgroundValue(0);

	//typedef itk::SubtractImageFilter<ByteImageType> SubtracterType;
	//SubtracterType::Pointer subtracter = SubtracterType::New();
	//subtracter->SetInput1(air);
	//subtracter->SetInput2(keeper->GetOutput());
	//subtracter->Update();
	//air = subtracter->GetOutput();
	//WriteITK <ByteImageType> (air,"airNoBkg.nii");

	PixelType intialThreshold = 180;

	// Find tagged
	typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> BinarizerType;
	BinarizerType::Pointer binarizer = BinarizerType::New();
	binarizer->SetInput(input);
	binarizer->SetOutsideValue(0);
	binarizer->SetInsideValue(255);
	binarizer->SetLowerThreshold(180);
	binarizer->Update();
	WriteITK <ByteImageType> (binarizer->GetOutput(),"allTagged.nii");
	
	// Remove small components
	typedef itk::BinaryShapeOpeningImageFilter<ByteImageType> BinaryOpeningType;
	BinaryOpeningType::Pointer opener = BinaryOpeningType::New();
	opener->SetInput(binarizer->GetOutput());
	opener->SetAttribute("Size");
	opener->SetForegroundValue(255);
	opener->SetBackgroundValue(0);
	opener->SetLambda(5);
	opener->Update();
	WriteITK <ByteImageType> (opener->GetOutput(),"allTaggedRemovedSmall.nii");

	ByteImageType::Pointer tag = opener->GetOutput();

	/*
	// Fill holes
	typedef itk::BinaryBallStructuringElement<unsigned char, 2> StructuringElementType;
	StructuringElementType se;
	se.SetRadius(1);
	se.CreateStructuringElement();

	typedef itk::BinaryMorphologicalClosingImageFilter<ByteImageType,ByteImageType,StructuringElementType> ClosingType;
	ClosingType::Pointer closer = ClosingType::New();
	closer->SetInput(opener->GetOutput());
	closer->SetKernel(se);
	closer->SetForegroundValue(255);
	closer->Update();

	ByteImageType::Pointer tag = closer->GetOutput();
	WriteITK <ByteImageType> (closer->GetOutput(),"allTaggedRemovedSmall.nii");
	*/

	//// Overlay air and tagged
	//typedef itk::OrImageFilter<ByteImageType> OrImageFilterType;
	//OrImageFilterType::Pointer orFilter = OrImageFilterType::New();
	//orFilter->SetInput1(air);
	//orFilter->SetInput2(tag);

	//// Fill holes
	//closer = ClosingType::New();
	//closer->SetInput(orFilter->GetOutput());
	//closer->SetKernel(se);
	//closer->SetForegroundValue(255);
	//closer->Update();

	//ByteImageType::Pointer at = closer->GetOutput();
	//WriteITK <ByteImageType> (at,"at.nii");

	// Setup otsu
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(input);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1200);
	otsuCalculator->SetNumberOfHistogramBins(128);

	// Setup threshold image
	ImageType::Pointer thresholdImage = ImageType::New();
	thresholdImage->SetRegions(region);
	thresholdImage->CopyInformation(input);
	thresholdImage->Allocate();

	// Set initial threshold
	thresholdImage->FillBuffer(intialThreshold);

	// Get label attributes
	typedef itk::BinaryImageToShapeLabelMapFilter<ByteImageType> BinaryImageToShapeLabelMapFilterType;
	typedef BinaryImageToShapeLabelMapFilterType::OutputImageType LabelMapType;
	typedef LabelMapType::LabelObjectType LabelObjectType;

	BinaryImageToShapeLabelMapFilterType::Pointer labeler = BinaryImageToShapeLabelMapFilterType::New();
	labeler->SetInput(tag);
	labeler->SetInputForegroundValue(255);
	labeler->SetOutputBackgroundValue(0);
	labeler->Update();
	LabelMapType::Pointer labelMap = labeler->GetOutput();

	// Convert label map to label image for viewing
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,IntImageType> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer converter = LabelMapToLabelImageFilterType::New();
	converter->SetInput(labelMap);
	converter->Update();
	WriteITK <IntImageType> (converter->GetOutput(),"labelImage.nii");

	unsigned int padding = 5;

	for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
	{
		const LabelObjectType * labelObject = labelMap->GetLabelObject( label );
		
		// Get and label region
		ByteImageType::RegionType labelRegion = labelObject->GetRegion();
		ByteImageType::IndexType labelIndex = labelRegion.GetIndex();
		ByteImageType::SizeType labelSize = labelRegion.GetSize();

		// Expand region in XY plane
		
		for ( unsigned int j=0; j<2; j++ )
		{
			labelIndex[j] = labelIndex[j] - padding > 0 ? labelIndex[j] - padding : 0;

			labelSize[j] += 2*padding;

			int over = labelIndex[j] + labelSize[j] > region.GetSize()[j];

			if (over > 0)
				labelSize[j] -= over;
		}

		labelRegion.SetIndex(labelIndex);
		labelRegion.SetSize(labelSize);

		// Get otsu for sub region
		otsuCalculator->SetRegion(labelRegion);

		std::stringstream ss;
		ss << label << "_intensity.csv";

		//otsuCalculator->SetPrintHistogram(ss.str());

		otsuCalculator->Compute();

		PixelType ot = otsuCalculator->GetThreshold();

		std::cout << "Label: " << label << "\tOtsu: " << ot << std::endl;

		// Store largest threshold in image
		IteratorType it(thresholdImage,labelRegion);

		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if (ot > it.Get())
			{
				it.Set(ot);
			}
		}

	}

	WriteITK <ImageType> (thresholdImage,"thresholdImage.nii");

	// Get gradient magnitude
	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	FloatImageType::Pointer gm = gmFilter->GetOutput();
	WriteITK <FloatImageType> (gm,"gm.nii");

	// Allocate output classified image
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// Show classification with hard threshold at 200
	FloatIteratorType gmIt(gm,region);
	IteratorType inputIt(input,region);
	IteratorType thIt(thresholdImage,region);
	ByteIteratorType colonIt(colon,region);
	ByteIteratorType outIt(out,region);

	for (gmIt.GoToBegin(), inputIt.GoToBegin(), colonIt.GoToBegin(), outIt.GoToBegin(), thIt.GoToBegin(); !gmIt.IsAtEnd(); ++gmIt, ++inputIt, ++colonIt, ++outIt, ++thIt)
	{
		if (colonIt.Get() != 0)
		{
			PixelType I = inputIt.Get();
			PixelType T = 434;
			float G = gmIt.Get();

			if ( ( I >= T && G < 0.8*I ) || I > 1000 )
			{
				outIt.Set(1);
			} else if ( I <= -700 ) {
				outIt.Set(2);
			} else if ( I < T && I > -300 && G <= 300 ) {
				outIt.Set(3);
			}
		}
	}

	// Show classification with adaptive threshold
	WriteITK <ByteImageType> (out,"outPre.nii");

	for (gmIt.GoToBegin(), inputIt.GoToBegin(), colonIt.GoToBegin(), outIt.GoToBegin(), thIt.GoToBegin(); !gmIt.IsAtEnd(); ++gmIt, ++inputIt, ++colonIt, ++outIt, ++thIt)
	{
		if (colonIt.Get() != 0)
		{
			PixelType I = inputIt.Get();
			PixelType T = thIt.Get();
			float G = gmIt.Get();

			if ( ( I >= T && G < 0.8*I ) || I > 1000 )
			{
				outIt.Set(1);
			} else if ( I <= -700 ) {
				outIt.Set(2);
			} else if ( I < T && I > -300 && G <= 300 ) {
				outIt.Set(3);
			}
		}
	}

	WriteITK <ByteImageType> (out,"outPost.nii");

	system("pause"); 							
	return 0; 									
}