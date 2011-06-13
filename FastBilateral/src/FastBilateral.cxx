#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <itkFastBilateralImageFilter.h>
#include <itkSubtractImageFilter.h>
 												
int main(int argc, char * argv[])				
{ 						
	if( argc < 6 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " input startSlice endSlice domainSigma rangeSigma" << std::endl;
		return EXIT_FAILURE;
	}

	ImageType::Pointer input = ReadDicom <ImageType> (argv[1],atoi(argv[2]),atoi(argv[3]));
	WriteITK <ImageType> (input,"input.nii");

	typedef itk::FastBilateralImageFilter<ImageType,ImageType> FastBilateralImageFitlerType;
	FastBilateralImageFitlerType::Pointer bilateralFilter = FastBilateralImageFitlerType::New();
	bilateralFilter->SetInput(input);
	
	bilateralFilter->SetDomainSigma( atoi( argv[4] ) );
	bilateralFilter->SetRangeSigma( atoi( argv[5] ) );
	bilateralFilter->Update();

	std::stringstream ss;
	ss << "out_" << argv[4] << "_" << argv[5] << ".nii";

	WriteITK <ImageType> (bilateralFilter->GetOutput(),ss.str());
	
	typedef itk::SubtractImageFilter<ImageType,ImageType> SubtracterType;
	SubtracterType::Pointer 

	return 0; 									
} 												
