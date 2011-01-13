#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <iostream>
 
typedef short													PixelType;
typedef int														ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	

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

	// Threshold to detect bone
	ScalarImageType3D::Pointer threshold = ScalarImageType3D::New();
	threshold->SetRegions( input->GetLargestPossibleRegion() );
	threshold->Allocate();
	ScalarIteratorType thresholdIt( threshold, input->GetLargestPossibleRegion() );

	for (	inputIt.GoToBegin(), thresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++inputIt, ++thresholdIt		) 
	{
		if ( inputIt.Get() >= 200 )		{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> (threshold, "boneThreshold.hdr");

	// Run connected component filter to remove bone
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
	relabelFilter->SetMinimumObjectSize( 300000 );

	try
	{
		relabelFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabel = relabelFilter->GetOutput();

	std::cout << "Number of objects: " << relabelFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabel, "relabel.hdr");

	//// Threshold the relabel to only show top 10 components
	//ScalarIteratorType relabelIt( relabel, input->GetLargestPossibleRegion() );
	//for (	relabelIt.GoToBegin(); !relabelIt.IsAtEnd(); ++relabelIt ) 
	//{
	//	if ( relabelIt.Get() > 10 )		{ relabelIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> ( relabel, "relabelThreshold.hdr");

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
