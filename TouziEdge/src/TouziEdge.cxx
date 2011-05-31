#include "otbTouziEdgeDetectorImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "utils.h"

int main(int argc, char * argv[])
{
	// setup typedefs
	typedef  float         InternalPixelType;
	typedef  unsigned char OutputPixelType;

	typedef itk::Image<InternalPixelType,  2> InternalImageType;
	typedef itk::Image<OutputPixelType,  2>   OutputImageType;

	// load input
	InternalImageType::Pointer input = ReadDicom <InternalImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85);

	// rescale input
	typedef itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType> RescalerInputType;
	RescalerInputType::Pointer inputRescaler = RescalerInputType::New();
	inputRescaler->SetOutputMaximum(255);
	inputRescaler->SetOutputMinimum(0);
	inputRescaler->SetInput(input);
	inputRescaler->Update();
	input = inputRescaler->GetOutput();

	WriteITK <InternalImageType> (input,"input.nii");
	
	// define radius
	const unsigned int radius = 1;

	// setup filters
	typedef otb::TouziEdgeDetectorImageFilter<InternalImageType,InternalImageType> FilterType;

	FilterType::Pointer filter = FilterType::New();

	filter->SetInput(input);

	FilterType::SizeType Radius;
	Radius[0] = radius;
	Radius[1] = radius;

	filter->SetRadius(Radius);
	filter->Update();

	// rescale outputs
	typedef itk::RescaleIntensityImageFilter<InternalImageType, OutputImageType> RescalerOutputType;

	RescalerOutputType::Pointer rescaler = RescalerOutputType::New();

	rescaler->SetOutputMinimum(itk::NumericTraits<OutputPixelType>::min());
	rescaler->SetOutputMaximum(itk::NumericTraits<OutputPixelType>::max());

	rescaler->SetInput(filter->GetOutput());
	rescaler->Update();

	WriteITK <OutputImageType> (rescaler->GetOutput(),"edge.nii");

	rescaler->SetInput(filter->GetOutputDirection());
	rescaler->Update();

	WriteITK <OutputImageType> (rescaler->GetOutput(),"direction.nii");

	return EXIT_SUCCESS;
}
