#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <iostream>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRegionalMaximaImageFilter.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>

typedef float													PixelType;
typedef int											ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D;
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >							ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >							ScalarImageType3D;	

typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;


template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );

int writeCount = 1;

int main() 
{ 
	// Setup IO
	typedef itk::ImageFileReader< ImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( "C:/ImageData/mr10_092_13p.i0344.hdr" );

	try
	{
		reader->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ImageType3D::Pointer input = reader->GetOutput();
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );

	WriteITK <ImageType3D> (input, "input.hdr");

	// Threshold to detect air lumen
	ScalarImageType3D::Pointer threshold = ScalarImageType3D::New();
	threshold->SetRegions( input->GetLargestPossibleRegion() );
	threshold->Allocate();
	ScalarIteratorType thresholdIt( threshold, input->GetLargestPossibleRegion() );

	for (	inputIt.GoToBegin(), thresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++inputIt, ++thresholdIt		) 
	{
		if ( inputIt.Get() <= -600 || inputIt.Get() >= 200 )	{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> ( threshold, "thresholdInput.hdr");

	// Create binary ball structuring element
	typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 4 );
    structuringElement.CreateStructuringElement();

	typedef itk::BinaryMorphologicalClosingImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryClosingFilterType;
	BinaryClosingFilterType::Pointer closingFilter = BinaryClosingFilterType::New();
	closingFilter->SetInput( threshold );
	closingFilter->SetKernel( structuringElement );
	closingFilter->SetForegroundValue( 1 );

	try {
		closingFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer closed = closingFilter->GetOutput();
	WriteITK <ScalarImageType3D> ( closed, "closedThreshold.hdr");

	/*
	std::cout << "Starting connected component" << std::endl;

	// Run connected component filter to remove background
	typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
	ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	connectedComponentFilter->SetInput( threshold );

	try
	{
		connectedComponentFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	WriteITK <ScalarImageType3D> ( connectedComponentFilter->GetOutput(), "connectedComponent.hdr");

	std::cout << "Starting relabel connecting component" << std::endl;

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
	//relabelFilter->SetMinimumObjectSize( 20 );

	try
	{
		relabelFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabel = relabelFilter->GetOutput();

	WriteITK <ScalarImageType3D> ( relabel, "relabel.hdr");
	*/

	system("pause"); 
	return 0; 
} 

template< class T >
void WriteITK (typename T::Pointer image , std::string ss )
{
	typedef itk::ImageFileWriter< T >	WriterType;
	WriterType::Pointer writer = WriterType::New();

	std::stringstream ss2;
	ss2 << writeCount++ << "_" << ss;
	writer->SetFileName(ss2.str().c_str());

	writer->SetInput(image);
	writer->GlobalWarningDisplayOff();
	writer->ReleaseDataFlagOn();
	std::cerr<<"Writing: "<< ss <<std::endl;
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return;
	}
}
