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
typedef itk::ImageRegionIteratorWithIndex< ImageType2D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType2D >	ByteIteratorType;

typedef itk::ImageSliceIteratorWithIndex< ImageType3D >			InputSliceIteratorType;
typedef itk::ImageSliceIteratorWithIndex< ByteImageType3D >		OutputSliceIteratorType;


template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );

void DoStuff( ImageType2D::Pointer image, LinearIteratorType imageIt);

int writeCount = 1;

int main( int argc, char *argv[] )
{
/*
	Performs heterogeneous stool colon segmentation using the 
	ITK ConnectedThresholdImageFilter
*/
	

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
	InputSliceIteratorType inputIt( input, input->GetLargestPossibleRegion() );

	// Setup 2D region from 3D input
	ImageType2D::RegionType region;
	ImageType2D::RegionType::SizeType size;
	ImageType2D::RegionType::IndexType index;

	ImageType3D::RegionType requestedRegion = input->GetRequestedRegion();

	index[ 0 ] = requestedRegion.GetIndex( 0 );
	index[ 1 ] = requestedRegion.GetIndex( 1 );
	size[ 0 ] = requestedRegion.GetSize()[ 0 ];
	size[ 1 ] = requestedRegion.GetSize()[ 1 ];

	region.SetSize( size );
	region.SetIndex( index );

	// Create 2D slice image and iterator
	ImageType2D::Pointer slice = ImageType2D::New();
	slice->SetRegions( region );
	slice->Allocate();

	LinearIteratorType sliceIt( slice, slice->GetLargestPossibleRegion() );

	// Setup 3D output image
	ImageType3D::Pointer output = ImageType3D::New();
	output->SetRegions( requestedRegion );
	output->Allocate();
	OutputSliceIteratorType outputIt( output, output->GetLargestPossibleRegion() );

	// Setup directions
	inputIt.SetFirstDirection( 0 );
	inputIt.SetSecondDirection( 1 );

	outputIt.SetFirstDirection( 0 );
	outputIt.SetSecondDirection( 1 );

	inputIt.GoToBegin();
	outputIt.GoToBegin();

	int sliceCounter = 1;
	
	while ( !inputIt.IsAtEnd() )	// Iterate by slice
	{
		sliceIt.GoToBegin();
		
		// Make 2D slice
		while ( !inputIt.IsAtEndOfSlice() )
		{
			while ( !inputIt.IsAtEndOfLine() )
			{
				sliceIt.Set( inputIt.Get() );

				++sliceIt;
				++inputIt;
			}

			inputIt.NextLine();
		}

		inputIt.NextSlice();

		// Perform operations on slice
		DoStuff(slice, sliceIt);

		// Place 2D slice in 3D output
		sliceIt.GoToBegin();

		while ( !outputIt.IsAtEndOfSlice() )
		{
			while ( !outputIt.IsAtEndOfLine() )
			{

				outputIt.Set( (byte) sliceIt.Get() );
				
				++sliceIt;
				++outputIt;
			}

			outputIt.NextLine();
		}
		
		outputIt.NextSlice();
		std::cout << "slice: " << sliceCounter++ << std::endl;
	}

	WriteITK <ImageType3D> (output, "output.hdr");

	system("PAUSE");
	return 0;
}

void FindColon( ImageType2D::Pointer input) {
	// Set regions
	ImageType2D::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType2D::IndexType endIndex = fullRegion.GetIndex();
	ImageType2D::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);

	// Threshold to detect air lumen
	ByteImageType2D::Pointer threshold = ByteImageType::New();
	threshold->SetRegions( fullRegion );
	threshold->Allocate();
	ByteIteratorType thresholdIt( threshold, fullRegion );

	for (	inputIt.GoToBegin(), thresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++inputIt, ++thresholdIt		) 
	{
		if ( inputIt.Get() <= -600 )	{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	//WriteITK( threshold, "thresholdInput.hdr");

	//// Run reigonal maximal filter
	//typedef itk::RegionalMaximaImageFilter< ImageType, ByteImageType> RegionalMaximaFilterType;
	//RegionalMaximaFilterType::Pointer regionalMaximaFilter = RegionalMaximaFilterType::New();
	//regionalMaximaFilter->SetInput( input );

	//try
	//{
	//	regionalMaximaFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//WriteITK( regionalMaximaFilter->GetOutput(), "regionalMaxima.hdr");

	//// Run connected component filter to remove background
	//typedef itk::ConnectedComponentImageFilter< ByteImageType, ByteImageType > ConnectedComponentFilterType;
	//ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	//connectedComponentFilter->SetInput( threshold );

	//try
	//{
	//	connectedComponentFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//WriteITK( connectedComponentFilter->GetOutput(), "connectedComponent.hdr");

	//// Relabel components
	//typedef itk::RelabelComponentImageFilter< ByteImageType, ByteImageType > RelabelFilterType;
	//RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	//relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
	//relabelFilter->SetMinimumObjectSize( 20 );

	//try
	//{
	//	relabelFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ByteImageType::Pointer relabel = relabelFilter->GetOutput();
	//ByteIteratorType relabelIt( relabel, fullRegion );

	//WriteITK( relabel, "relabelComponent.hdr");

	//// Remove outer air and colon labels by thresholding >= 2
	//for (	relabelIt.GoToBegin(), thresholdIt.GoToBegin();
	//		!relabelIt.IsAtEnd() && !thresholdIt.IsAtEnd();
	//		++relabelIt, ++thresholdIt		) 
	//{
	//	if ( relabelIt.Get() >= 2 )		{ thresholdIt.Set( 1 ); }
	//	else							{ thresholdIt.Set( 0 ); }
	//}

	//ByteImageType::Pointer airMask = threshold;
	//ByteIteratorType airMaskIt( airMask, airMask->GetLargestPossibleRegion() );

	//WriteITK( airMask, "relabelComponentTreshold.hdr");

	//// Apply region growing filter to detect tagged regions >= 200
	//typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType > ConnectedTresholdFilterType;
	//ConnectedTresholdFilterType::Pointer connectedThresholdFilter = ConnectedTresholdFilterType::New();
	//connectedThresholdFilter->SetLower( 200 );
	////connectedThresholdFilter->SetUpper( max );
	//connectedThresholdFilter->SetReplaceValue( 1 );

	//// Find edges of air
	//for (	airMaskIt.GoToBegin(), inputIt.GoToBegin();
	//		!airMaskIt.IsAtEnd() && !inputIt.IsAtEnd();
	//		++airMaskIt, ++inputIt		) 
	//{	
	//	if (airMaskIt.Get() == 1) {	// at air

	//		ImageType::IndexType index = airMaskIt.GetIndex();
	//		for(int i=-1;i<=1;i++) {
	//			if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
	//				for (int j=-1;j<=1;j++) {
	//					if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
	//						for (int k=-1;k<=1;k++) {
	//							if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
	//								
	//								ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};

	//								if ( airMask->GetPixel( neighborIndex ) == 0 )	// if neighbor is non-air
	//								{
	//									//airMask->SetPixel( index, 2); // mark as edge
	//									connectedThresholdFilter->AddSeed( index );
	//									inputIt.Set(200);
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//connectedThresholdFilter->SetInput( input );

	//try
	//{
	//	connectedThresholdFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ByteImageType::Pointer tagged = connectedThresholdFilter->GetOutput();
	//ByteIteratorType taggedIt( tagged, tagged->GetLargestPossibleRegion() );

	//WriteITK( tagged, "taggedRegionGrowing.hdr");

	//// Combine air and tagged regions into one image
	//for (	airMaskIt.GoToBegin(), taggedIt.GoToBegin();
	//		!airMaskIt.IsAtEnd() && !taggedIt.IsAtEnd();
	//		++airMaskIt, ++taggedIt		) 
	//{
	//	if ( airMaskIt.Get() == 1)
	//	{
	//		taggedIt.Set( 1 );
	//	}
	//}

	//WriteITK( tagged, "segmentedColon.hdr");

	//// Dilate
	//
	//// Create binary ball structuring element
	//typedef itk::BinaryBallStructuringElement< ByteImageType::PixelType, 3> StructuringElementType;
	//StructuringElementType structuringElement;
 //   structuringElement.SetRadius( 4 );
 //   structuringElement.CreateStructuringElement();

	//typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType > BinaryDilateFilterType;
	//BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	//dilateFilter->SetInput( tagged );
	//dilateFilter->SetKernel( structuringElement );
	//dilateFilter->SetDilateValue( 1 );

	//try
	//{
	//	dilateFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//WriteITK( dilateFilter->GetOutput() , "segmentedColonDilated.hdr");
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


/*

void WriteITK(ImageType2D::Pointer image, std::string ss) {
	std::cout << "Writing " << ss << std::endl;
	
	typedef itk::ImageFileWriter< ImageType2D > WriterType;
	WriterType::Pointer writer = WriterType::New();
	
	std::stringstream ss2;
	ss2 << writeCount++ << "_" << ss;
	writer->SetFileName(ss2.str().c_str());
	writer->SetInput(image);
	writer->GlobalWarningDisplayOff();

	try  
	{
		writer->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	}
}

void WriteITK(ByteImageType::Pointer image, std::string ss) {
	std::cout << "Writing " << ss << std::endl;
	
	typedef itk::ImageFileWriter< ByteImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	
	std::stringstream ss2;
	ss2 << writeCount++ << "_" << ss;
	writer->SetFileName(ss2.str().c_str());
	writer->SetInput(image);
	writer->GlobalWarningDisplayOff();

	try  
	{
		writer->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	} 
}

*/

void DoStuff( ImageType2D::Pointer image, LinearIteratorType imageIt) {

	//WriteITK < ImageType2D > (image, "pre_slice.hdr");

	for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); imageIt.NextLine() )
	{
		imageIt.GoToBeginOfLine();

		while( ! imageIt.IsAtEndOfLine() )
		{
			imageIt.Set( 100 );
			++imageIt;
		}
	}

	//WriteITK < ImageType2D > (image, "post_slice.hdr");
	
}