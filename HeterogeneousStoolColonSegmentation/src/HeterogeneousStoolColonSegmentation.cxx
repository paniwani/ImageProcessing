#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <iostream>
#include "itkBinaryThresholdImageFilter.h";

int main( int argc, char *argv[] )
{
/*
	Performs heterogeneous stool colon segmentation using the 
	ITK ConnectedTresholdImageFilter

	Parameters:
	1) Input image
	2) Output image
*/
	typedef itk::Image< float, 3 > ImageType;

	// Setup IO
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();

	reader->SetFileName( argv[1] );
	writer->SetFileName( argv[2] );

	try
	{
		reader->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

	ImageType::Pointer input = reader->GetOutput();
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );


	// Find min/max intensities
	typedef itk::MinimumMaximumImageCalculator< ImageType > MinimumMaximumType;
	MinimumMaximumType::Pointer minMaxFilter = MinimumMaximumType::New();
	minMaxFilter->SetImage( input );
	minMaxFilter->Compute();
	ImageType::PixelType min = minMaxFilter->GetMinimum();
	ImageType::PixelType max = minMaxFilter->GetMaximum();
	std::cout << "Minimum: " << min << std::endl; 
	std::cout << "Maximum: " << max << std::endl; 

	// Threshold to detect air lumen
	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType> ThresholdFilterType;
	ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
	thresholdFilter->SetInput( input );
	thresholdFilter->SetOutsideValue( 0 );
	thresholdFilter->SetInsideValue( 1 );
	thresholdFilter->SetLowerThreshold( min );
	thresholdFilter->SetUpperThreshold( -600 );

	try
	{
		thresholdFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ImageType::Pointer airImage = thresholdFilter->GetOutput();

	
	IteratorType airImageIt( airImage, airImage->GetLargestPossibleRegion() );






	writer->SetInput( thresholdFilter->GetOutput() );

	// Apply region growing filter
	typedef itk::ConnectedThresholdImageFilter< ImageType, ImageType > ConnectedFilterType;
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	connectedThreshold->SetInput( input );
	connectedThreshold->SetLower( 200 );
	connectedThreshold->SetUpper( max );
	connectedThreshold->SetReplaceValue( max*2 );

	for (	airImageIt.GoToBegin(), inputIt.GoToBegin();
			!airImageIt.IsAtEnd() && !inputIt.IsAtEnd();
			++airImageIt, ++inputIt		) 
	{
		if ( airImageIt.Get() == 1 ) {
			connectedThreshold->AddSeed( airImageIt.GetIndex() );
			inputIt.Set(200);
		}		
	}

	try
	{
		connectedThreshold->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	/*

	// Multiple Seeds
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType inputIt( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

	for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
	{
		if (inputIt.Get() <= -600)
		{
			connectedThreshold->AddSeed( inputIt.GetIndex() );
		}
	}

	*/


	writer->SetInput( connectedThreshold->GetOutput() );

	try
	{
		writer->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	system("PAUSE");
	return 0;
}