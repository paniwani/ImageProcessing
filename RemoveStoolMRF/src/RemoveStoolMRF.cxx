#include "itkImage.h"
#include "itkFixedArray.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkScalarToArrayCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMRFImageFilter.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageClassifierBase.h"
#include "itkBinaryShapeKeepNObjectsImageFilter.h"

template <typename T>
typename T::Pointer ReadITK(char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< T > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );

	try {
		reader->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error reading image: " << err << std::endl;
		return NULL;
	}

	return reader->GetOutput();
}

template <typename T>
void WriteITK(typename T::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< T >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	
	try {
		writer->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error writing image: " << err << std::endl;
	}
}


int main( int argc, char * argv [] )
{
  /*if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0];
    std::cerr << " inputScalarImage inputLabeledImage";
    std::cerr << " outputLabeledImage numberOfIterations";
    std::cerr << " smoothingFactor numberOfClasses";
    std::cerr << " mean1 mean2 ... meanN " << std::endl;
    return EXIT_FAILURE;
    }*/

	// typedefs
	typedef signed short        PixelType;
	const unsigned int          Dimension = 2;
	typedef itk::Image<PixelType, Dimension > ImageType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	typedef unsigned char       LabelPixelType;
	typedef itk::Image<LabelPixelType, Dimension > LabelImageType;
	typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;

	typedef itk::FixedArray<LabelPixelType,1>  ArrayPixelType;
	typedef itk::Image< ArrayPixelType, Dimension > ArrayImageType;

	typedef itk::ScalarToArrayCastImageFilter< ImageType, ArrayImageType > ScalarToArrayFilterType;

	typedef itk::MRFImageFilter< ArrayImageType, LabelImageType > MRFFilterType;

	// load inputs
	ImageType::Pointer input = ReadITK <ImageType> ("input.nii");
	LabelImageType::Pointer map = ReadITK <LabelImageType> ("vmapPreQR.nii");

	// rescale input
	typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleInputType;
	RescaleInputType::Pointer inputRescaler = RescaleInputType::New();
	inputRescaler->SetInput(input);
	inputRescaler->SetOutputMaximum(255);
	inputRescaler->SetOutputMinimum(0);
	inputRescaler->Update();
	input = inputRescaler->GetOutput();

	WriteITK <ImageType> (input,"input.nii");

	// convert outer unclassified to air class
	IteratorType inputIt(input,input->GetLargestPossibleRegion());
	LabelIteratorType mapIt(map,map->GetLargestPossibleRegion());

	LabelImageType::Pointer unclass = LabelImageType::New();
	unclass->SetRegions(map->GetLargestPossibleRegion());
	unclass->Allocate();
	unclass->FillBuffer(0);
	LabelIteratorType unclassIt(unclass,unclass->GetLargestPossibleRegion());

	for (mapIt.GoToBegin(), unclassIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt, ++unclassIt)
	{
		if (mapIt.Get() == 4)
		{
			unclassIt.Set(255);
		} else {
			unclassIt.Set(0);
		}
	}

	typedef itk::BinaryShapeKeepNObjectsImageFilter<LabelImageType> KeeperFilterType;
	KeeperFilterType::Pointer keeper = KeeperFilterType::New();
	keeper->SetInput(unclass);
	keeper->SetAttribute("Size");
	keeper->SetNumberOfObjects(1);
	keeper->SetBackgroundValue(0);
	keeper->SetForegroundValue(255);
	keeper->Update();
	unclass = keeper->GetOutput();
	unclassIt = LabelIteratorType(unclass,unclass->GetLargestPossibleRegion());

	for (mapIt.GoToBegin(), unclassIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt, ++unclassIt)
	{
		if (unclassIt.Get() == 255)
		{
			mapIt.Set(2);
		}
	}

	for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
	{
		/*if (mapIt.Get() == 4)
			mapIt.Set(2);
		if (mapIt.Get() == 5)
			mapIt.Set(4);
		if (mapIt.Get() == 6)
			mapIt.Set(5);
		if (mapIt.Get() == 7)
			mapIt.Set(6);*/

		mapIt.Set(mapIt.Get()-1);
	}

	WriteITK <LabelImageType> (map,"vmap2.nii");

	const unsigned int numOfClasses = 4;

	// get means of all classes
	double mean[4] = {0,0,0,0};
	int count[4] = {0,0,0,0};

	for (inputIt.GoToBegin(), mapIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++mapIt)
	{
		count[ mapIt.Get() ]++;
		mean [ mapIt.Get() ] += inputIt.Get(); 
	}

	for (int i=0; i<numOfClasses; i++)
	{
		mean[i] /= (double) count[i];
		std::cout << "Class " << i*255/(numOfClasses-1) << " mean: " << mean[i] << std::endl;
	}

	mean[3] = 255;

	// convert scalar input to array image
	ScalarToArrayFilterType::Pointer scalarToArrayFilter = ScalarToArrayFilterType::New();
	scalarToArrayFilter->SetInput( input );	

	// setup mrf filter
	MRFFilterType::Pointer mrfFilter = MRFFilterType::New();
	mrfFilter->SetInput( scalarToArrayFilter->GetOutput() );
	mrfFilter->SetNumberOfClasses( numOfClasses );
	mrfFilter->SetMaximumNumberOfIterations( 50 );
	mrfFilter->SetErrorTolerance( 1e-7 );
	mrfFilter->SetSmoothingFactor( 3 );

	// setup classifier
	typedef itk::ImageClassifierBase<ArrayImageType, LabelImageType > SupervisedClassifierType;
	SupervisedClassifierType::Pointer classifier = SupervisedClassifierType::New();

	// setup decision rule
	typedef itk::MinimumDecisionRule DecisionRuleType;
	DecisionRuleType::Pointer  classifierDecisionRule = DecisionRuleType::New();
	classifier->SetDecisionRule( classifierDecisionRule.GetPointer() );

	// setup membership function
	typedef itk::Statistics::DistanceToCentroidMembershipFunction< ArrayPixelType > MembershipFunctionType;
	typedef MembershipFunctionType::Pointer MembershipFunctionPointer;
	
	double meanDistance = 0;

	itk::Array<double> centroid(1);

	for( unsigned int i=0; i < numOfClasses; i++ )
	{
		MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();

		centroid[0] = mean[i]; 

		membershipFunction->SetCentroid( centroid );

		classifier->AddMembershipFunction( membershipFunction );
		meanDistance += static_cast< double > (centroid[0]);
	}
	meanDistance /= numOfClasses;

	mrfFilter->SetNeighborhoodRadius( 1 );

	std::vector< double > weights;
	weights.push_back(1.5);
	weights.push_back(2.0);
	weights.push_back(1.5);
	weights.push_back(2.0);
	weights.push_back(0.0); // This is the central pixel
	weights.push_back(2.0);
	weights.push_back(1.5);
	weights.push_back(2.0);
	weights.push_back(1.5);

	// Software Guide : EndCodeSnippet


	// Software Guide : BeginLatex
	// We now scale weights so that the smoothing function and the image fidelity
	// functions have comparable value. This is necessary since the label
	// image and the input image can have different dynamic ranges. The fidelity
	// function is usually computed using a distance function, such as the
	// \doxygen{DistanceToCentroidMembershipFunction} or one of the other 
	// membership functions. They tend to have values in the order of the means
	// specified. 
	// Software Guide : EndLatex 

	// Software Guide : BeginCodeSnippet

	double totalWeight = 0;
	for(std::vector< double >::const_iterator wcIt = weights.begin(); 
	  wcIt != weights.end(); ++wcIt )
	{
	totalWeight += *wcIt;
	}
	for(std::vector< double >::iterator wIt = weights.begin(); 
	  wIt != weights.end(); wIt++ )
	{
	*wIt = static_cast< double > ( (*wIt) * meanDistance / (2 * totalWeight));
	}

	mrfFilter->SetMRFNeighborhoodWeight( weights );
	mrfFilter->SetClassifier( classifier );

	typedef MRFFilterType::OutputImageType  OutputImageType;

	// Rescale outputs to the dynamic range of the display
	typedef itk::Image< unsigned char, Dimension > RescaledOutputImageType;
	typedef itk::RescaleIntensityImageFilter< 
			 OutputImageType, RescaledOutputImageType >   RescalerType;

	RescalerType::Pointer intensityRescaler = RescalerType::New();
	intensityRescaler->SetOutputMinimum(   0 );
	intensityRescaler->SetOutputMaximum( 255 );
	intensityRescaler->SetInput( mrfFilter->GetOutput() );  

	// Software Guide : BeginCodeSnippet
	typedef itk::ImageFileWriter< OutputImageType > WriterType;

	WriterType::Pointer writer = WriterType::New();

	writer->SetInput( intensityRescaler->GetOutput() );

	writer->SetFileName( "output.nii" );

	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << " image file : " << argv[2] << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}
	// Software Guide : EndCodeSnippet

	std::cout << "Number of Iterations : ";
	std::cout << mrfFilter->GetNumberOfIterations() << std::endl;
	std::cout << "Stop condition: " << std::endl;
	std::cout << "  (1) Maximum number of iterations " << std::endl;
	std::cout << "  (2) Error tolerance:  "  << std::endl;
	std::cout << mrfFilter->GetStopCondition() << std::endl;
	
	system("pause");
	return EXIT_SUCCESS;
}


