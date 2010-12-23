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

typedef float													PixelType;
typedef unsigned char											BytePixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D;
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< BytePixelType, 2 >							ByteImageType2D;
typedef itk::Image< BytePixelType, 3 >							ByteImageType3D;	

typedef itk::ImageLinearIteratorWithIndex< ImageType2D >		LinearIteratorType;
typedef itk::ImageLinearIteratorWithIndex< ByteImageType2D >	LinearByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType2D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType2D >	ByteIteratorType;

typedef itk::ImageSliceIteratorWithIndex< ImageType3D >			InputSliceIteratorType;
typedef itk::ImageSliceIteratorWithIndex< ByteImageType3D >		OutputSliceIteratorType;


template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );

int writeCount = 1;

int main() 
{ 
	// Setup IO
	typedef itk::ImageFileReader< ByteImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( "1_segmentedColonFinal3D_Closing1.hdr" );

	try
	{
		reader->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//ByteImageType3D::Pointer segmentedColon = reader->GetOutput();

	//WriteITK <ByteImageType3D> (segmentedColon, "input.hdr");

	std::cout << "Starting connected component" << std::endl;

	// Run connected component filter to remove background
	typedef itk::ConnectedComponentImageFilter< ByteImageType3D, ByteImageType3D > ConnectedComponentFilterType;
	ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	connectedComponentFilter->SetInput( reader->GetOutput() );

	try
	{
		connectedComponentFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	WriteITK <ByteImageType3D> ( connectedComponentFilter->GetOutput(), "connectedComponent.hdr");

	std::cout << "Starting relabel connecting component" << std::endl;

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ByteImageType3D, ByteImageType3D > RelabelFilterType;
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

	ByteImageType3D::Pointer relabel = relabelFilter->GetOutput();

	WriteITK <ByteImageType3D> ( relabel, "relabel.hdr");

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
