/*
Demonstrates bug with itkBSplineInterpolateImageFunction 
when image start index is nonzero
*/

#include "itkImage.h"
#include <itkBSplineInterpolateImageFunction.h>
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

typedef float PixelType;
typedef itk::Image<PixelType, 2> ImageType2D;
typedef itk::ImageRegionIteratorWithIndex< ImageType2D > IteratorType2D;

template <typename T>
void WriteITK(typename T::Pointer image, std::string name);

int main()
{
	// Create test image
	ImageType2D::Pointer image = ImageType2D::New();
	
	ImageType2D::IndexType index;
	index[0] = 0;
	index[1] = 0;

	ImageType2D::SizeType size;
	size[0] = 50;
	size[1] = 50;

	ImageType2D::RegionType region;
	region.SetIndex( index );
	region.SetSize( size );

	image->SetRegions( size );
	image->Allocate();
	image->FillBuffer(0);

	IteratorType2D image_iter(image,region);

	int count = 0;
	
	// Make a gradient from 0-255
	for (image_iter.GoToBegin(); !image_iter.IsAtEnd(); ++image_iter)
		image_iter.Set( count++ * 255 / ( size[0] * size[1] ) );

	WriteITK <ImageType2D> (image,"image.nii");

	// Extract arbitrary region
	ImageType2D::RegionType extractRegion;

	index[0] = 10;
	index[1] = 10;

	size[0] = 30;
	size[1] = 30;
	
	extractRegion.SetIndex( index );
	extractRegion.SetSize( size );

	typedef itk::ExtractImageFilter<ImageType2D,ImageType2D> ExtractImageFilterType;
	ExtractImageFilterType::Pointer extracter = ExtractImageFilterType::New();
	extracter->SetInput( image );
	extracter->SetExtractionRegion( extractRegion );
	extracter->Update();
	image = extracter->GetOutput();

	WriteITK <ImageType2D> (image, "extract.nii");

	// Interpolate
	typedef itk::BSplineInterpolateImageFunction<ImageType2D> BSplineInterpolateImageFunctionType;
	BSplineInterpolateImageFunctionType::Pointer interpolator = BSplineInterpolateImageFunctionType::New();
	interpolator->SetSplineOrder( 3 );
	interpolator->SetInputImage( image );
	
	BSplineInterpolateImageFunctionType::ContinuousIndexType startIndex = interpolator->GetStartIndex();
	BSplineInterpolateImageFunctionType::ContinuousIndexType endIndex = interpolator->GetEndIndex();

	std::cout << "Start Index: " << startIndex[0] << " " << startIndex[1] << std::endl;
	std::cout << "End Index: " << endIndex[0] << " " << endIndex[1] << std::endl;

	// Test interpolator at non-interpolated indices
	for (image_iter.GoToBegin(); !image_iter.IsAtEnd(); ++image_iter)
		std::cout << interpolator->EvaluateAtContinuousIndex( (BSplineInterpolateImageFunctionType::ContinuousIndexType) image_iter.GetIndex() ) << std::endl;

	system("pause");
	return 0;
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