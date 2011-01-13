#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <iostream>

#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
 
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
	// Read input
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
	ImageType3D::RegionType region = input->GetLargestPossibleRegion();
	ScalarImageType3D::SizeType size = region.GetSize();
	ImageType3D::SpacingType spacing = input->GetSpacing();
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );
	WriteITK <ImageType3D> (input, "input.hdr");

	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();
	reader2->SetFileName( "3_relabelAirThreshold.hdr" );
	reader2->Update();
	ScalarImageType3D::Pointer test = reader2->GetOutput();
	test->SetSpacing( spacing );

	/*
	// Threshold to detect air
	ScalarImageType3D::Pointer airThreshold = ScalarImageType3D::New();
	airThreshold->SetRegions( region );
	airThreshold->SetSpacing( spacing );
	airThreshold->Allocate();
	ScalarIteratorType airThresholdIt( airThreshold, region );

	for (	inputIt.GoToBegin(), airThresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++inputIt, ++airThresholdIt		) 
	{
		if ( inputIt.Get() <= -600 )	{ airThresholdIt.Set( 1 ); }
		else							{ airThresholdIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> (airThreshold, "airThreshold.hdr");

	// Run connected component filter to remove lungs from air threshold
	typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
	ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	connectedComponentFilter->SetInput( airThreshold );

	try
	{
		connectedComponentFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//WriteITK <ScalarImageType3D> ( connectedComponentFilter->GetOutput(), "connectedComponentStool.hdr");

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
	relabelFilter->SetMinimumObjectSize( 5000 );

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

	//WriteITK <ScalarImageType3D> ( relabel, "relabelAirThreshold.hdr");
	*/

	// Convert label image to label map
	typedef unsigned long LabelType;	
	typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
	typedef itk::LabelMap< LabelObjectType > LabelMapType;

	typedef itk::LabelImageToShapeLabelMapFilter< ScalarImageType3D, LabelMapType > 	ConverterType;
	ConverterType::Pointer converter = ConverterType::New();
	converter->SetInput( test );
	//converter->SetInput( relabel );
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();

	std::cout << "#\tPhysicalSize\tCentroid" << std::endl;

	for( unsigned int label=2; label<=labelMap->GetNumberOfLabelObjects(); label++ )	//ignore labels 0 and 1 (body and external air)
	{
		const LabelObjectType * lo = labelMap->GetLabelObject( label );
		
		float centroidX = lo->GetCentroid()[0] / spacing[0];
		float centroidY = size[1] - (lo->GetCentroid()[2] / spacing[1]);
		float ps = lo->GetPhysicalSize();
		std::cout << label << "\t" << lo->GetPhysicalSize() << "\t[" << centroidX << ", " << centroidY << "]\t";

		if ( centroidY < size[1]/2 && ps > 100000) {
			std::cout << "FOUND LUNG";
		}

		std::cout << std::endl;
	}


	//// Threshold to detect stool
	//ScalarImageType3D::Pointer stoolThreshold = ScalarImageType3D::New();
	//stoolThreshold->SetRegions( region );
	//stoolThreshold->Allocate();
	//ScalarIteratorType stoolThresholdIt( stoolThreshold, region );

	//for (	inputIt.GoToBegin(), stoolThresholdIt.GoToBegin();
	//		!inputIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
	//		++inputIt, ++stoolThresholdIt		) 
	//{
	//	if ( inputIt.Get() >= 200 )		{ stoolThresholdIt.Set( 1 ); }
	//	else							{ stoolThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> (stoolThreshold, "stoolThreshold.hdr");

	//// Run connected component filter to remove bone from stool threshold
	//typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
	//ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	//connectedComponentFilter->SetInput( stoolThreshold );

	//try
	//{
	//	connectedComponentFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	////WriteITK <ScalarImageType3D> ( connectedComponentFilter->GetOutput(), "connectedComponentStool.hdr");

	//// Relabel components
	//typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;
	//RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	//relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
	//relabelFilter->SetMinimumObjectSize( 300000 );

	//try
	//{
	//	relabelFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ScalarImageType3D::Pointer relabel = relabelFilter->GetOutput();

	//std::cout << "Number of objects: " << relabelFilter->GetNumberOfObjects() << std::endl;

	//WriteITK <ScalarImageType3D> ( relabel, "relabelStoolThreshold.hdr");

	//// Remove bone from stool threshold
	//ScalarIteratorType relabelIt( relabel, input->GetLargestPossibleRegion() );
	//for (	relabelIt.GoToBegin(), stoolThresholdIt.GoToBegin();
	//		!relabelIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
	//		++relabelIt, ++stoolThresholdIt		) 
	//{
	//	if ( relabelIt.Get() == 1 )		{ stoolThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> ( stoolThreshold, "stoolThresholdWithoutBone.hdr");

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
