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

typedef itk::Image< float, 3 > ImageType;
void WriteITK(ImageType::Pointer image, std::string ss, int count);

int writeCount = 1;

int main( int argc, char *argv[] )
{
/*
	Performs heterogeneous stool colon segmentation using the 
	ITK ConnectedTresholdImageFilter

	Parameters:
	1) Input image
*/
	

	// Setup IO
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( argv[1] );

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

	// Set regions
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);

	// Run reigonal maximal filter
	typedef itk::RegionalMaximaImageFilter< ImageType, ImageType> RegionalMaximaFilterType;
	RegionalMaximaFilterType::Pointer regionalMaximaFilter = RegionalMaximaFilterType::New();
	regionalMaximaFilter->SetInput( input );

	try
	{
		regionalMaximaFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	WriteITK( regionalMaximaFilter->GetOutput(), "regionalMaxima.hdr", writeCount++);

	// Threshold to detect air lumen
	ImageType::Pointer threshold = ImageType::New();
	threshold->SetRegions( fullRegion );
	threshold->Allocate();
	IteratorType thresholdIt( threshold, fullRegion );

	for (	inputIt.GoToBegin(), thresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++inputIt, ++thresholdIt		) 
	{
		if ( inputIt.Get() <= -600 )	{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	WriteITK( threshold, "thresholdInput.hdr", writeCount++);

	// Run connected component filter to remove background
	typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConnectedComponentFilterType;
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

	WriteITK( connectedComponentFilter->GetOutput(), "connectedComponent.hdr", writeCount++);

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ImageType, ImageType > RelabelFilterType;
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
	relabelFilter->SetMinimumObjectSize( 20 );

	try
	{
		relabelFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ImageType::Pointer relabel = relabelFilter->GetOutput();
	IteratorType relabelIt( relabel, fullRegion );

	WriteITK( relabel, "relabelComponent.hdr", writeCount++);

	// Remove outer air and colon labels by thresholding >= 2
	for (	relabelIt.GoToBegin(), thresholdIt.GoToBegin();
			!relabelIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++relabelIt, ++thresholdIt		) 
	{
		if ( relabelIt.Get() >= 2 )		{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	ImageType::Pointer airMask = threshold;
	IteratorType airMaskIt( airMask, airMask->GetLargestPossibleRegion() );

	WriteITK( airMask, "relabelComponentTreshold.hdr", writeCount++);

	// Apply region growing filter to detect tagged regions >= 200
	typedef itk::ConnectedThresholdImageFilter< ImageType, ImageType > ConnectedTresholdFilterType;
	ConnectedTresholdFilterType::Pointer connectedThresholdFilter = ConnectedTresholdFilterType::New();
	connectedThresholdFilter->SetLower( 200 );
	//connectedThresholdFilter->SetUpper( max );
	connectedThresholdFilter->SetReplaceValue( 1 );

	// Find edges of air
	for (	airMaskIt.GoToBegin(), inputIt.GoToBegin();
			!airMaskIt.IsAtEnd() && !inputIt.IsAtEnd();
			++airMaskIt, ++inputIt		) 
	{	
		if (airMaskIt.Get() == 1) {	// at air

			ImageType::IndexType index = airMaskIt.GetIndex();
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							for (int k=-1;k<=1;k++) {
								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
									ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};

									if ( airMask->GetPixel( neighborIndex ) == 0 )	// if neighbor is non-air
									{
										//airMask->SetPixel( index, 2); // mark as edge
										connectedThresholdFilter->AddSeed( index );
										inputIt.Set(200);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	connectedThresholdFilter->SetInput( input );

	try
	{
		connectedThresholdFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ImageType::Pointer tagged = connectedThresholdFilter->GetOutput();
	IteratorType taggedIt( tagged, tagged->GetLargestPossibleRegion() );

	WriteITK( tagged, "taggedRegionGrowing.hdr", writeCount++);

	// Combine air and tagged regions into one image
	for (	airMaskIt.GoToBegin(), taggedIt.GoToBegin();
			!airMaskIt.IsAtEnd() && !taggedIt.IsAtEnd();
			++airMaskIt, ++taggedIt		) 
	{
		if ( airMaskIt.Get() == 1)
		{
			taggedIt.Set( 1 );
		}
	}

	WriteITK( tagged, "segmentedColon.hdr", writeCount++);

	// Dilate
	
	// Create binary ball structuring element
	typedef itk::BinaryBallStructuringElement< ImageType::PixelType, 3> StructuringElementType;
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 4 );
    structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( tagged );
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

	WriteITK( dilateFilter->GetOutput() , "segmentedColonDilated.hdr", writeCount++);

	system("PAUSE");
	return 0;
}

void WriteITK(ImageType::Pointer image, std::string ss, int count) {
	std::cout << "Writing " << ss << std::endl;
	
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	
	std::stringstream ss2;
	ss2 << count << "_" << ss;
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