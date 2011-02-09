#include "itkImage.h"
#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkGradientImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>

typedef itk::Image< float, 3 > ImageType;
typedef itk::CovariantVector< float,3> CovariantVectorType;
typedef itk::Image< CovariantVectorType, 3> ImageVectorType;

typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorTypeFloat4WithIndex;
typedef itk::ImageRegionIterator<ImageVectorType>    IteratorImageVectorType;

typedef itk::GradientImageFilter<ImageType> GradientFilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientMagnitudeFilterType;

void ReadITK(ImageType::Pointer &image, char * fileName);

int main()
{
	// Read input
	ImageType::Pointer input = ImageType::New();
	ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_100-105.hdr");

	GradientMagnitudeFilterType::Pointer gmFilter = GradientMagnitudeFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	ImageType::Pointer gm = gmFilter->GetOutput();
	gmFilter.~SmartPointer();
	IteratorTypeFloat4WithIndex gm_iter(gm, gm->GetLargestPossibleRegion() );

	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->Update();
	ImageVectorType::Pointer gradient = gradientFilter->GetOutput();
	IteratorImageVectorType gradient_iter(gradient, gm->GetLargestPossibleRegion() );

	std::ofstream myfile;
	myfile.open("gradient.txt");
	myfile << "gradient_magnitude\tgradient\n";

	for (gm_iter.GoToBegin(), gradient_iter.GoToBegin();
		 !gm_iter.IsAtEnd() && !gradient_iter.IsAtEnd();
		 ++gm_iter, ++gradient_iter)
	{
		CovariantVectorType grad = gradient_iter.Get();
		myfile << gm_iter.Get() << "\t" << grad.GetNorm() << "\n";
	}

	myfile.close();

	//system("pause");
	return 0;
}

void ReadITK(ImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	image = reader->GetOutput();
}