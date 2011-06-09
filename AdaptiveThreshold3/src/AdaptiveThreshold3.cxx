#include <itkImage.h> 						
#include <iostream> 							
#include <utils.h> 
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


int main(int argc, char * argv[])				
{ 		
	// Load image
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm/mr10_092_13p_i0246.dcm");
	WriteITK <ImageType2D> (input,"input.nii");

	ImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// Segment colon
	typedef itk::ColonSegmentationFilter< ImageType2D, ByteImageType2D > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( input );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType2D::Pointer colon = colonSegmenter->GetOutput();
	WriteITK <ByteImageType2D> (colon,"colon.nii");

	// Mask input with colon
	typedef itk::MaskImageFilter<ImageType2D,ByteImageType2D> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(input);
	masker->SetInput2(colon);
	masker->SetOutsideValue(-1025);
	masker->Update();
	input = masker->GetOutput();
	WriteITK <ImageType2D> (input,"inputMasked.nii");

	//// Find all air
	//typedef itk::BinaryThresholdImageFilter<ImageType2D,ByteImageType2D> BinarizerType;
	//BinarizerType::Pointer binarizer = BinarizerType::New();
	//binarizer->SetInput(input);
	//binarizer->SetOutsideValue(0);
	//binarizer->SetInsideValue(255);
	//binarizer->SetUpperThreshold(-600);
	//binarizer->Update();
	//ByteImageType2D::Pointer air = binarizer->GetOutput();
	//WriteITK <ByteImageType2D> (air,"air.nii");

	//// Subtract background
	//typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType2D> BinaryKeeperType;
	//BinaryKeeperType::Pointer keeper = BinaryKeeperType::New();
	//keeper->SetInput(air);
	//keeper->SetAttribute("Size");
	//keeper->SetNumberOfObjects(1);
	//keeper->SetForegroundValue(255);
	//keeper->SetBackgroundValue(0);

	//typedef itk::SubtractImageFilter<ByteImageType2D> SubtracterType;
	//SubtracterType::Pointer subtracter = SubtracterType::New();
	//subtracter->SetInput1(air);
	//subtracter->SetInput2(keeper->GetOutput());
	//subtracter->Update();
	//air = subtracter->GetOutput();
	//WriteITK <ByteImageType2D> (air,"airNoBkg.nii");

	// Find tagged
	typedef itk::BinaryThresholdImageFilter<ImageType2D,ByteImageType2D> BinarizerType;
	BinarizerType::Pointer binarizer = BinarizerType::New();
	binarizer->SetInput(input);
	binarizer->SetOutsideValue(0);
	binarizer->SetInsideValue(255);
	binarizer->SetLowerThreshold(200);
	
	// Remove small components
	typedef itk::BinaryShapeOpeningImageFilter<ByteImageType2D> BinaryOpeningType;
	BinaryOpeningType::Pointer opener = BinaryOpeningType::New();
	opener->SetInput(binarizer->GetOutput());
	opener->SetAttribute("Size");
	opener->SetForegroundValue(255);
	opener->SetBackgroundValue(0);
	opener->SetLambda(5);

	// Fill holes
	typedef itk::BinaryBallStructuringElement<unsigned char, 2> StructuringElementType;
	StructuringElementType se;
	se.SetRadius(1);
	se.CreateStructuringElement();

	typedef itk::BinaryMorphologicalClosingImageFilter<ByteImageType2D,ByteImageType2D,StructuringElementType> ClosingType;
	ClosingType::Pointer closer = ClosingType::New();
	closer->SetInput(opener->GetOutput());
	closer->SetKernel(se);
	closer->SetForegroundValue(255);
	closer->Update();

	ByteImageType2D::Pointer tag = closer->GetOutput();
	WriteITK <ByteImageType2D> (tag,"tag.nii");

	//// Overlay air and tagged
	//typedef itk::OrImageFilter<ByteImageType2D> OrImageFilterType;
	//OrImageFilterType::Pointer orFilter = OrImageFilterType::New();
	//orFilter->SetInput1(air);
	//orFilter->SetInput2(tag);

	//// Fill holes
	//closer = ClosingType::New();
	//closer->SetInput(orFilter->GetOutput());
	//closer->SetKernel(se);
	//closer->SetForegroundValue(255);
	//closer->Update();

	//ByteImageType2D::Pointer at = closer->GetOutput();
	//WriteITK <ByteImageType2D> (at,"at.nii");

	// Setup otsu
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType2D > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(input);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1500);

	// Setup thresholdImage image
	ImageType2D::Pointer thresholdImage = ImageType2D::New();
	thresholdImage->SetRegions(region);
	thresholdImage->CopyInformation(input);
	thresholdImage->Allocate();
	thresholdImage->FillBuffer(0);

	// Get label attributes
	typedef itk::BinaryImageToShapeLabelMapFilter<ByteImageType2D> BinaryImageToShapeLabelMapFilterType;
	typedef BinaryImageToShapeLabelMapFilterType::OutputImageType LabelMapType;
	typedef LabelMapType::LabelObjectType LabelObjectType;

	BinaryImageToShapeLabelMapFilterType::Pointer labeler = BinaryImageToShapeLabelMapFilterType::New();
	labeler->SetInput(tag);
	labeler->SetInputForegroundValue(255);
	labeler->SetOutputBackgroundValue(0);
	labeler->Update();
	LabelMapType::Pointer labelMap = labeler->GetOutput();

	// Convert label map to label image for viewing
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,IntImageType2D> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer converter = LabelMapToLabelImageFilterType::New();
	converter->SetInput(labelMap);
	converter->Update();
	WriteITK <IntImageType2D> (converter->GetOutput(),"labelImage.nii");

	unsigned int padding = 10;

	for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
	{
		const LabelObjectType * labelObject = labelMap->GetLabelObject( label );
		
		// Get and expand label region
		ByteImageType2D::RegionType labelRegion = labelObject->GetRegion();
		ByteImageType2D::IndexType labelIndex = labelRegion.GetIndex();
		ByteImageType2D::SizeType labelSize = labelRegion.GetSize();
		
		for ( unsigned int j=0; j<ByteImageType2D::ImageDimension; j++ )
		{
			labelIndex[j] = labelIndex[j] - padding > 0 ? labelIndex[j] - 10 : 0;

			labelSize[j] += 2*padding;

			int over = labelIndex[j] + labelSize[j] > region.GetSize()[j];

			if (over > 0)
				labelSize[j] -= over;
		}

		labelRegion.SetIndex(labelIndex);
		labelRegion.SetSize(labelSize);

		// Get otsu for that region
		otsuCalculator->SetRegion(labelRegion);

		std::stringstream ss;
		ss << label << "_intensity.csv";

		otsuCalculator->SetPrintHistogram(ss.str());

		otsuCalculator->Compute();

		// Store threshold in image
		IteratorType2D it(thresholdImage,labelRegion);

		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
		it.Set(otsuCalculator->GetThreshold());
		}

	}

	WriteITK <ImageType2D> (thresholdImage,"thresholdImage.nii");

	system("pause"); 							
	return 0; 									
}