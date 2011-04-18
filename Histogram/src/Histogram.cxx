#include <itkImage.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkColonSegmentationFilter.h>
#include <itkMaskImageFilter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <utils.h>

int main()
{
	// Load image
	const ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr9/mr2_002_13p.i0393/dcm");
	WriteITK <ImageType> (input,"input.hdr");

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilter;
	ColonSegmentationFilter::Pointer segmenter = ColonSegmentationFilter::New();
	segmenter->SetInput( input );
	segmenter->Update();
	ByteImageType::Pointer colon = segmenter->GetOutput();
	ByteIteratorType colon_iter(colon,region);

	WriteITK <ByteImageType> (colon,"colon_segmentation.hdr");

	// Mask the input with colon segmentation
	typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( input );
	masker->SetInput2( colon );
	masker->SetOutsideValue( -1500 ); // arbitrary
	masker->Update();
	ImageType::Pointer input_mask = masker->GetOutput();

	WriteITK <ImageType> (input_mask,"input_masked.hdr");

	// Compute otsu threshold for tissue and stool intensity
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType >  OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-250);
	otsuCalculator->SetHistogramMax(1400);
	otsuCalculator->SetPrintHistogram("intensity.csv");
	otsuCalculator->Compute();

	std::cout << "Otsu Intensity: " << otsuCalculator->GetThreshold() << std::endl;

	// Get gradient magnitude
	typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gradientFilter = GradientMagnitudeImageFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->Update();
	ImageType::Pointer gradient = gradientFilter->GetOutput();

	WriteITK <ImageType> (gradient, "gradient.hdr");

	// Mask gradient magnitude
	masker = MaskImageFilterType::New();
	masker->SetInput1( gradient );
	masker->SetInput2( colon );
	masker->SetOutsideValue( -10 ); // arbitrary
	masker->Update();
	ImageType::Pointer gradient_mask = masker->GetOutput();
	WriteITK <ImageType> (gradient_mask,"gradient_masked.hdr");

	// Get gradient max
	typedef itk::MinimumMaximumImageCalculator<ImageType> RangeCalculator;
	RangeCalculator::Pointer rangeCalculator = RangeCalculator::New();
	rangeCalculator->SetImage( gradient_mask );
	rangeCalculator->Compute();
	float gradient_max = rangeCalculator->GetMaximum();

	std::cout << "Gradient Max: " << gradient_max << std::endl;

	// Compute otsu threshold for gradient magnitude
	otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(gradient_mask);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(0);
	otsuCalculator->SetHistogramMax(gradient_max);
	otsuCalculator->SetPrintHistogram("gradient.csv");
	otsuCalculator->Compute();

	std::cout << "Otsu Gradient: " << otsuCalculator->GetThreshold() << std::endl;


	//// Threshold to only include Tissue and Stool regions
	//float hMin = -250;
	//float hMax = 1400;

	//double totalPixels = 0;

	//IteratorType input_iter(input,region);
	//for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	//{
	//	if ( input_iter.Get() < hMin+1 || input_iter.Get() > hMax )
	//	{
	//		input_iter.Set( -1500 );
	//	} else {
	//		totalPixels++;
	//	}
	//}

	//WriteITK <ImageType> (input,"input_threshold.hdr");

	//std::cout << "Total number of pixels: " << totalPixels << std::endl;

	//// Get histogram
	//typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > HistogramGeneratorType;
	//HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
	//histogramGenerator->SetInput( input );
	//histogramGenerator->SetNumberOfBins( hMax-hMin );
	//histogramGenerator->SetMarginalScale( 10.0 );
	//histogramGenerator->SetHistogramMin( hMin-0.5 );
	//histogramGenerator->SetHistogramMax( hMax+0.5 );
	//histogramGenerator->Compute();

	//typedef HistogramGeneratorType::HistogramType  HistogramType;
	//const HistogramType * histogram = histogramGenerator->GetOutput();
	//const unsigned int histogramSize = histogram->Size();
	//std::cout << "Histogram size " << histogramSize << std::endl;

	//// Write histogram to text file
	//std::ofstream file;
	//file.open("histogram.csv");
	//file << "Bin,Frequency\n";

	//int bin;
	//for( bin=1; bin < histogramSize; bin++ )
	//{
	//	file << hMin+bin << "," << histogram->GetFrequency( bin, 0 ) << "\n";
	//	//std::cout << "bin = " << bin << " frequency = ";
	//	//std::cout << histogram->GetFrequency( bin, 0 ) << std::endl;
	//}

	//
	//typedef itk::OtsuThresholdImageFilter< ImageType, ByteImageType >  OtsuThresholdImageFilterType;
	//OtsuThresholdImageFilterType::Pointer otsuFilter = OtsuThresholdImageFilterType::New();
	//otsuFilter->SetOutsideValue( 0 );
	//otsuFilter->SetInsideValue( 1 );
	//otsuFilter->SetNumberOfHistogramBins( hMax-hMin );
	//// Software Guide : EndCodeSnippet


	////  Software Guide : BeginLatex
	////  
	////  The execution of the filter is triggered by invoking the \code{Update()}
	////  method.   If the filter's output has been passed as input to subsequent
	////  filters, the \code{Update()} call on any posterior filters in the
	////  pipeline will indirectly trigger the update of this filter.
	////
	////  Software Guide : EndLatex 

	//// Software Guide : BeginCodeSnippet
	//filter->Update();
	//// Software Guide : EndCodeSnippet



	////  Software Guide : BeginLatex
	////  
	////  We print out here the Threshold value that was computed internally by the
	////  filter. For this we invoke the \code{GetThreshold} method.
	////
	////  Software Guide : EndLatex 

	//// Software Guide : BeginCodeSnippet
	//int threshold = filter->GetThreshold();
	//std::cout << "Threshold = " << threshold << std::endl;
	//// Software Guide : EndCodeSnippet



	//file.close();
	system("pause");
	return 0;
}