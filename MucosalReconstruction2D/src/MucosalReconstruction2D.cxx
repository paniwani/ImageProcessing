#include "itkImage.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkShapeRelabelImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkGaussianBlurImageFunction.h>


typedef itk::Image< float, 2 > ImageType;
typedef itk::Image< unsigned char, 2 > ByteImageType;

typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageRegionIterator< ByteImageType >  ByteIteratorType;


ImageType::Pointer ReadITK(char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );

	try {
		reader->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}

	ImageType::Pointer image = reader->GetOutput();

	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;

	image->SetDirection( direction );
	
	return image;
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(ImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(IntImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< IntImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

int main()
{
	ImageType::Pointer input = ReadITK("C:/ImageData/Modified_mr10_092_13p.i0344_85_100_40_slice13_subtracted.hdr");

	WriteITK(input,"input.hdr");

	// Threshold to detect air
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing( input->GetSpacing() );
	air->Allocate();

	IteratorType input_iter(input,region);
	ByteIteratorType air_iter(air,region);

	for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter)
	{
		if (input_iter.Get() < -300)
		{
			air_iter.Set( 0 );
		} else {
			air_iter.Set( 1 );
		}
	}

	WriteITK(air,"air.hdr");

	// Remove non-air components that are not connected to the largest body of tissue
	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType > BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryFilter->SetInput( air );
	binaryFilter->SetBackgroundValue(0);
	binaryFilter->SetForegroundValue(1);
	binaryFilter->SetAttribute("Size");
	binaryFilter->SetNumberOfObjects(1);
	binaryFilter->Update();
	air = binaryFilter->GetOutput();

	WriteITK(air,"air_large_only.hdr");

	// Invert image to label air components
	air_iter = ByteIteratorType(air,region);
	for (air_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter)
	{
		air_iter.Set( !air_iter.Get() );
	}	

	// Isolate largest air component with largest size on border (background air)
	binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryFilter->SetInput( air );
	binaryFilter->SetBackgroundValue(0);
	binaryFilter->SetForegroundValue(1);
	binaryFilter->SetAttribute("SizeOnBorder");
	binaryFilter->SetNumberOfObjects(1);
	binaryFilter->Update();
	ByteImageType::Pointer bkg = binaryFilter->GetOutput();
	ByteIteratorType bkg_iter(bkg,region);
	
	for (air_iter.GoToBegin(), bkg_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++bkg_iter)
	{
		if (bkg_iter.Get() == 1)
			air_iter.Set( 0 );
	}

	WriteITK(air,"air_no_bkg.hdr");

	// Copy input to output and remove non-air
	ImageType::Pointer output = ImageType::New();
	output->SetRegions( region );
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation( input );
	output->Allocate();
	IteratorType output_iter(output,region);

	for (input_iter.GoToBegin(), output_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter, ++air_iter)
	{
		output_iter.Set( input_iter.Get() );

		if (air_iter.Get() == 1)
		{
			output_iter.Set( -1025 );
		}
	}

	/****
		Determine initial starting tissue value To by dilating air boundary, finding edges, and averaging intensity
	****/

	// Dilate by 1
	typedef itk::BinaryBallStructuringElement<unsigned char, 2> StructuringElementType;
	
	ByteImageType::SizeType radius;
	radius[0] = 1;
	radius[1] = 1;

	StructuringElementType structuringElement;
	structuringElement.SetRadius( radius );
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput( air );
	dilateFilter->SetForegroundValue(1);
	dilateFilter->SetBackgroundValue(0);
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->Update();
	ByteImageType::Pointer air_dilate = dilateFilter->GetOutput();

	WriteITK(air_dilate,"air_dilate_1px.hdr");

	// Find edges of air
	typedef itk::BinaryContourImageFilter< ByteImageType, ByteImageType > BinaryContourImageFilterType;
	BinaryContourImageFilterType::Pointer binaryEdgeFilter = BinaryContourImageFilterType::New();
	binaryEdgeFilter->SetInput( air_dilate );
	binaryEdgeFilter->SetForegroundValue(1);
	binaryEdgeFilter->SetBackgroundValue(0);
	binaryEdgeFilter->Update();

	ByteImageType::Pointer air_edge = binaryEdgeFilter->GetOutput();

	WriteITK(air_edge,"air_edge.hdr");

	// Calculate average starting tissue value (T0)
	ByteIteratorType air_edge_iter(air_edge,region);

	float T0 = 0;
	unsigned int count = 0;

	for (input_iter.GoToBegin(), air_edge_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_edge_iter)
	{
		if ( air_edge_iter.Get() == 1 )
		{
			if ( input_iter.Get() >= -250 && input_iter.Get() <= 150 )
			{
				T0 += input_iter.Get();
				count++;
			}
		}
	}

	T0 /= count;

	/***
		Erode air mask for each layer, find eges, and apply scaffolded tissue value
	***/

	// Make mask of set of scaffolding pixels
	ByteImageType::Pointer scaffold_mask = ByteImageType::New();
	scaffold_mask->SetRegions( region );
	scaffold_mask->SetSpacing( input->GetSpacing() );
	scaffold_mask->Allocate();
	scaffold_mask->FillBuffer(0);
	ByteIteratorType scaffold_mask_iter(scaffold_mask,region);

	float T = T0;

	unsigned int numOfLayers = 2;
	
	for (int i=0; i < numOfLayers; i++)
	{
		T += (-1000-T0)/(numOfLayers+1);

		if (i>0)
		{
			typedef itk::BinaryErodeImageFilter< ByteImageType, ByteImageType, StructuringElementType > BinaryErodeImageFilterType;
			BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
			erodeFilter->SetInput( air );
			erodeFilter->SetForegroundValue(1);
			erodeFilter->SetBackgroundValue(0);
			erodeFilter->SetKernel( structuringElement );
			erodeFilter->Update();
			air = erodeFilter->GetOutput();
		}

		binaryEdgeFilter = BinaryContourImageFilterType::New();
		binaryEdgeFilter->SetInput( air );
		binaryEdgeFilter->SetForegroundValue(1);
		binaryEdgeFilter->SetBackgroundValue(0);
		binaryEdgeFilter->SetFullyConnected(true);
		binaryEdgeFilter->Update();
		air_edge = binaryEdgeFilter->GetOutput();

		air_edge_iter = ByteIteratorType(air_edge,region);

		for (output_iter.GoToBegin(), air_edge_iter.GoToBegin(), scaffold_mask_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++air_edge_iter, ++scaffold_mask_iter)
		{
			if (air_edge_iter.Get() == 1)
			{
				output_iter.Set( T );

				scaffold_mask_iter.Set( 1 );
			}
		}
	}

	WriteITK(output,"output_unsmoothed.hdr");

	WriteITK(scaffold_mask,"scaffold_mask.hdr");

	// Dilate scaffold mask
	dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput( scaffold_mask );
	dilateFilter->SetForegroundValue(1);
	dilateFilter->SetBackgroundValue(0);
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->Update();
	scaffold_mask = dilateFilter->GetOutput();
	scaffold_mask_iter = ByteIteratorType(scaffold_mask,region);
	
	WriteITK(scaffold_mask,"scaffold_mask_dilated.hdr");

	// Smooth output within scaffold mask
	typedef itk::GaussianBlurImageFunction< ImageType > GaussianBlurImageFunctionType;
	GaussianBlurImageFunctionType::Pointer gaussianFunction = GaussianBlurImageFunctionType::New();
	gaussianFunction->SetInputImage( output );

	GaussianBlurImageFunctionType::ErrorArrayType setError;
	setError.Fill( 0.01 );
	gaussianFunction->SetMaximumError( setError );
	gaussianFunction->SetSigma( 0.7 );

	for (output_iter.GoToBegin(), scaffold_mask_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++scaffold_mask_iter)
	{
		if (scaffold_mask_iter.Get() == 1)
		{
			output_iter.Set( gaussianFunction->EvaluateAtIndex( output_iter.GetIndex() ) );
		}
	}

	WriteITK(output,"output_smoothed.hdr");

	system("pause");
	return 0;
}

