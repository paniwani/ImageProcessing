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
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMedianImageFilter.h>
 
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



	/*// Test images
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();
	reader2->SetFileName( "trial1data/4_airThresholdWithoutLungs.hdr" );
	reader2->Update();
	ScalarImageType3D::Pointer airThreshold = reader2->GetOutput();
	airThreshold->SetSpacing( spacing );

	ReaderType2::Pointer reader3 = ReaderType2::New();
	reader3->SetFileName( "trial1data/7_stoolThresholdWithoutBone.hdr" );
	reader3->Update();
	ScalarImageType3D::Pointer stoolThreshold = reader3->GetOutput();
	stoolThreshold->SetSpacing( spacing );*/
	
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
	ConnectedComponentFilterType::Pointer connectedComponentAirFilter = ConnectedComponentFilterType::New();
	connectedComponentAirFilter->SetInput( airThreshold );

	try
	{
		connectedComponentAirFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//WriteITK <ScalarImageType3D> ( connectedComponentAirFilter->GetOutput(), "connectedComponentAir.hdr");

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;
	RelabelFilterType::Pointer relabelAirFilter = RelabelFilterType::New();
	relabelAirFilter->SetInput( connectedComponentAirFilter->GetOutput() );
	relabelAirFilter->SetMinimumObjectSize( 5000 );

	try
	{
		relabelAirFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabelAir = relabelAirFilter->GetOutput();

	std::cout << "Number of objects: " << relabelAirFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelAir, "relabelAirThreshold.hdr");
	
	// Convert label image to label map
	typedef unsigned long LabelType;	
	typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
	typedef itk::LabelMap< LabelObjectType > LabelMapType;

	typedef itk::LabelImageToShapeLabelMapFilter< ScalarImageType3D, LabelMapType > 	ConverterType;
	ConverterType::Pointer converter = ConverterType::New();
	//converter->SetInput( test );
	converter->SetInput( relabelAir );
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();

	std::cout << "#\tPhysicalSize\tCentroid\tModifiedCentroid" << std::endl;

	// Find labels for lungs and store
	std::vector< unsigned int > lungLabel;

	for( unsigned int label=2; label<=labelMap->GetNumberOfLabelObjects(); label++ )	//ignore labels 0 and 1 (body and external air)
	{
		const LabelObjectType * lo = labelMap->GetLabelObject( label );
		
		float centroidX = lo->GetCentroid()[0] / spacing[0];
		float centroidY = lo->GetCentroid()[1] / spacing[1];
		float ps = lo->GetPhysicalSize();
		std::cout << label << "\t" << lo->GetPhysicalSize() << "\t" << lo->GetCentroid() << "\t[" << centroidX << ", " << centroidY << "]\t";

		if ( centroidY > size[1]/2 && ps > 100000) {
			std::cout << "FOUND LUNG";
			lungLabel.push_back( label );
		}

		std::cout << std::endl;
	}

	// Remove external air and lungs from air threshold
	ScalarIteratorType relabelAirIt( relabelAir, region );

	for (	relabelAirIt.GoToBegin(), airThresholdIt.GoToBegin();
			!relabelAirIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++relabelAirIt, ++airThresholdIt		) 
	{
		if ( relabelAirIt.Get() == 1 )		{ airThresholdIt.Set( 0 ); }

		for (int i = 0; i < lungLabel.size() ; i++ ) {
			if ( relabelAirIt.Get() == lungLabel[i] ) {
				airThresholdIt.Set( 0 );
			}
		}
	}

	WriteITK <ScalarImageType3D> (airThreshold, "airThresholdWithoutLungs.hdr");


	// Threshold to detect stool
	ScalarImageType3D::Pointer stoolThreshold = ScalarImageType3D::New();
	stoolThreshold->SetRegions( region );
	stoolThreshold->Allocate();
	ScalarIteratorType stoolThresholdIt( stoolThreshold, region );

	for (	inputIt.GoToBegin(), stoolThresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
			++inputIt, ++stoolThresholdIt		) 
	{
		if ( inputIt.Get() >= 200 )		{ stoolThresholdIt.Set( 1 ); }
		else							{ stoolThresholdIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> (stoolThreshold, "stoolThreshold.hdr");

	// Run connected component filter to remove bone from stool threshold
	ConnectedComponentFilterType::Pointer connectedComponentStoolFilter = ConnectedComponentFilterType::New();
	connectedComponentStoolFilter->SetInput( stoolThreshold );

	try
	{
		connectedComponentStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//WriteITK <ScalarImageType3D> ( connectedComponentStoolFilter->GetOutput(), "connectedComponentStool.hdr");

	// Relabel components
	RelabelFilterType::Pointer relabelStoolFilter = RelabelFilterType::New();
	relabelStoolFilter->SetInput( connectedComponentStoolFilter->GetOutput() );
	relabelStoolFilter->SetMinimumObjectSize( 300000 );

	try
	{
		relabelStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabelStool = relabelStoolFilter->GetOutput();

	std::cout << "Number of objects: " << relabelStoolFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelStool, "relabelStoolThreshold.hdr");

	// Remove bone from stool threshold
	ScalarIteratorType relabelStoolIt( relabelStool, region );
	for (	relabelStoolIt.GoToBegin(), stoolThresholdIt.GoToBegin();
			!relabelStoolIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
			++relabelStoolIt, ++stoolThresholdIt		) 
	{
		if ( relabelStoolIt.Get() == 1 )		{ stoolThresholdIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> ( stoolThreshold, "stoolThresholdWithoutBone.hdr");

	// Remove noise from stool threshold image using mean filter
	typedef itk::BinaryMedianImageFilter< ScalarImageType3D, ScalarImageType3D> BinaryMedianFilterType;
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	
	ScalarImageType3D::SizeType indexRadius;
	indexRadius[0] = 1;
	indexRadius[1] = 1;
	indexRadius[2] = 1;
	
	medianFilter->SetRadius( indexRadius );
	medianFilter->SetInput( stoolThreshold );
	medianFilter->SetForegroundValue( 1 );
	medianFilter->Update();

	ScalarImageType3D::Pointer stoolThresholdSmooth = ScalarImageType3D::New();
	stoolThresholdSmooth->SetRegions( region );
	stoolThresholdSmooth->Allocate();
	stoolThresholdSmooth = medianFilter->GetOutput();
	ScalarIteratorType stoolThresholdSmoothIt( stoolThresholdSmooth, region );

	WriteITK <ScalarImageType3D> ( stoolThresholdSmooth, "stoolThresholdWithoutBoneSmoothed.hdr");

	// Add stool regions to air image
	for (	stoolThresholdSmoothIt.GoToBegin(), airThresholdIt.GoToBegin();
			!stoolThresholdSmoothIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++stoolThresholdSmoothIt, ++airThresholdIt		) 
	{
		if ( stoolThresholdSmoothIt.Get() == 1 ) { airThresholdIt.Set( 1 ); }
	}

	WriteITK <ScalarImageType3D> ( airThreshold, "taggedAll.hdr");

	// Create binary ball structuring element
	typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 1 );
    structuringElement.CreateStructuringElement();

	// Dilate
	typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( airThreshold );
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->SetDilateValue( 1 );

	try
	{
		dilateFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer dilateAir = dilateFilter->GetOutput();

	WriteITK <ScalarImageType3D> ( dilateAir , "airThresholdDilated_1.hdr");

	// Run connected component to remove erroneus regions
	ConnectedComponentFilterType::Pointer connectedComponentAirStoolFilter = ConnectedComponentFilterType::New();
	connectedComponentAirStoolFilter->SetInput( airThreshold );

	try
	{
		connectedComponentAirStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//WriteITK <ScalarImageType3D> ( connectedComponentAirStoolFilter->GetOutput(), "connectedComponentAirStool.hdr");

	// Relabel components
	RelabelFilterType::Pointer relabelAirStoolFilter = RelabelFilterType::New();
	relabelAirStoolFilter->SetInput( connectedComponentAirStoolFilter->GetOutput() );
	relabelAirStoolFilter->SetMinimumObjectSize( 5000 );

	try
	{
		relabelAirStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabelAirStool = relabelAirStoolFilter->GetOutput();

	std::cout << "Number of objects: " << relabelAirStoolFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelAirStool, "relabelAirStool.hdr");

	// Remove erroneous regions
	ScalarIteratorType relabelAirStoolIt( relabelAirStool, region );
	ScalarIteratorType dilateAirIt( dilateAir, region );
	for (	relabelAirStoolIt.GoToBegin(), dilateAirIt.GoToBegin();
			!relabelAirStoolIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
			++relabelAirStoolIt, ++dilateAirIt		) 
	{
		if ( relabelAirStoolIt.Get() != 1 ) { dilateAirIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> ( dilateAir , "segmentedColon.hdr");

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
