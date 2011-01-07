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
//#include <itkRegionalMaximaImageFilter.h>
#include <itkImageSliceIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>

typedef short 													PixelType;
typedef short											ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D;
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >							ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >							ScalarImageType3D;	

typedef itk::ImageLinearIteratorWithIndex< ImageType2D >		LinearIteratorType;
typedef itk::ImageLinearIteratorWithIndex< ScalarImageType2D >	LinearByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType2D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType2D >	ByteIteratorType;

typedef itk::ImageSliceIteratorWithIndex< ImageType3D >			InputSliceIteratorType;
typedef itk::ImageSliceIteratorWithIndex< ScalarImageType3D >		OutputSliceIteratorType;


template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );

void FindColon( ImageType2D::Pointer input, ScalarImageType2D::Pointer &output );

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

	// Create input 2D slice image and iterator
	ImageType2D::Pointer inputSlice = ImageType2D::New();
	inputSlice->SetRegions( region );
	inputSlice->Allocate();

	LinearIteratorType inputSliceIt( inputSlice, inputSlice->GetLargestPossibleRegion() );

	// Create output 2D slice image and iterator
	ScalarImageType2D::Pointer outputSlice = ScalarImageType2D::New();
	outputSlice->SetRegions( region );
	outputSlice->Allocate();

	LinearByteIteratorType outputSliceIt( outputSlice, outputSlice->GetLargestPossibleRegion() );

	// Setup 3D output image
	ScalarImageType3D::Pointer output = ScalarImageType3D::New();
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
		inputSliceIt.GoToBegin();
		
		// Make 2D slice
		while ( !inputIt.IsAtEndOfSlice() )
		{
			while ( !inputIt.IsAtEndOfLine() )
			{
				inputSliceIt.Set( inputIt.Get() );

				++inputSliceIt;
				++inputIt;
			}

			inputIt.NextLine();
		}

		inputIt.NextSlice();

		//WriteITK <ImageType2D> ( inputSlice, "inputSlice.hdr" );

		// Perform operations on slice
		FindColon( inputSlice, outputSlice);

		//WriteITK <ScalarImageType2D> (outputSlice, "outputSlice.hdr");

		// Place 2D slice in 3D output
		outputSliceIt.GoToBegin();

		while ( !outputIt.IsAtEndOfSlice() )
		{
			while ( !outputIt.IsAtEndOfLine() )
			{

				outputIt.Set( outputSliceIt.Get() );
				
				++outputSliceIt;
				++outputIt;
			}

			outputIt.NextLine();
		}
		
		outputIt.NextSlice();
		std::cout << "slice: " << sliceCounter++ << std::endl;
	}

	WriteITK <ScalarImageType3D> (output, "segmentedColonFinal3D_Closing1.hdr");

	std::cout << "Starting connected component" << std::endl;

	// Run connected component filter in 3D to remove erroneous regions
	typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
	ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
	connectedComponentFilter->SetInput( output );

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

	system("PAUSE");
	return 0;
}

void FindColon( ImageType2D::Pointer input, ScalarImageType2D::Pointer &output ) {

	// Create input iterator
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );
	ByteIteratorType outputIt( output, output->GetLargestPossibleRegion() );

	// Set regions
	ImageType2D::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType2D::IndexType endIndex = fullRegion.GetIndex();
	ImageType2D::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);

	// Threshold to detect air lumen
	ScalarImageType2D::Pointer threshold = ScalarImageType2D::New();
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

	//WriteITK <ScalarImageType2D> ( threshold, "thresholdInput.hdr");

	// Run connected component filter to remove background
	typedef itk::ConnectedComponentImageFilter< ScalarImageType2D, ScalarImageType2D > ConnectedComponentFilterType;
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

	//WriteITK <ScalarImageType2D> ( connectedComponentFilter->GetOutput(), "connectedComponent.hdr");

	// Relabel components
	typedef itk::RelabelComponentImageFilter< ScalarImageType2D, ScalarImageType2D > RelabelFilterType;
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

	ScalarImageType2D::Pointer relabel = relabelFilter->GetOutput();
	ByteIteratorType relabelIt( relabel, fullRegion );

	//WriteITK <ScalarImageType2D> ( relabel, "relabelComponent.hdr");

	// Remove outer air and colon tissue labels by thresholding >= 2
	for (	relabelIt.GoToBegin(), thresholdIt.GoToBegin();
			!relabelIt.IsAtEnd() && !thresholdIt.IsAtEnd();
			++relabelIt, ++thresholdIt		) 
	{
		if ( relabelIt.Get() >= 2 )		{ thresholdIt.Set( 1 ); }
		else							{ thresholdIt.Set( 0 ); }
	}

	ScalarImageType2D::Pointer airMask = threshold;
	ByteIteratorType airMaskIt( airMask, airMask->GetLargestPossibleRegion() );

	//WriteITK <ScalarImageType2D> ( airMask, "relabelComponentTreshold.hdr");

	// Dilate air mask

	// Create binary ball structuring element
	typedef itk::BinaryBallStructuringElement< ScalarImageType2D::PixelType, 2> StructuringElementType;
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 4 );
    structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< ScalarImageType2D, ScalarImageType2D, StructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( airMask );
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

	ScalarImageType2D::Pointer dilateAir = dilateFilter->GetOutput();
	ByteIteratorType dilateAirIt( dilateAir, dilateAir->GetLargestPossibleRegion() );

	//WriteITK <ScalarImageType2D> ( dilateAir , "airMaskDilated.hdr");

	// Apply region growing filter to detect tagged regions >= 200
	typedef itk::ConnectedThresholdImageFilter< ImageType2D, ScalarImageType2D > ConnectedThresholdFilterType;
	ConnectedThresholdFilterType::Pointer connectedThresholdFilter = ConnectedThresholdFilterType::New();
	connectedThresholdFilter->SetLower( 200 );
	//connectedThresholdFilter->SetUpper( max );
	connectedThresholdFilter->SetReplaceValue( 1 );

	// Set edges of air as seeds
	for (	airMaskIt.GoToBegin(), inputIt.GoToBegin();
			!airMaskIt.IsAtEnd() && !inputIt.IsAtEnd();
			++airMaskIt, ++inputIt		) 
	{	
		if (airMaskIt.Get() == 1) {	// at air

			ImageType2D::IndexType index = airMaskIt.GetIndex();
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							
							//for (int k=-1;k<=1;k++) {
							//	if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
							ImageType2D::IndexType neighborIndex={index[0]+i,index[1]+j};//,index[2]+k};

									if ( airMask->GetPixel( neighborIndex ) == 0 )	// if neighbor is non-air
									{
										airMask->SetPixel( index, 2); // mark as edge
										connectedThresholdFilter->AddSeed( index );
										//inputIt.Set(200);
									}
								//}
							//}
						}
					}
				}
			}
		}
	}

	// Set dilated air region as seed in connected threshold

	for (	airMaskIt.GoToBegin(), dilateAirIt.GoToBegin();
			!airMaskIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
			++airMaskIt, ++dilateAirIt		) 
	{	
		if ( airMaskIt.Get() == 0 && dilateAirIt.Get() == 1 ) 
		{
			ImageType2D::IndexType index = airMaskIt.GetIndex();
			airMask->SetPixel( index , 2);
			connectedThresholdFilter->AddSeed( index );
		}
	}

	// Write updated air mask to show region growing seeds
	//WriteITK <ScalarImageType2D> ( airMask, "airMaskWithSeeds.hdr");


	connectedThresholdFilter->SetInput( input );

	try
	{
		connectedThresholdFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType2D::Pointer tagged = connectedThresholdFilter->GetOutput();
	ByteIteratorType taggedIt( tagged, tagged->GetLargestPossibleRegion() );

	//WriteITK <ScalarImageType2D> ( tagged, "taggedRegionGrowing.hdr");

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

	//WriteITK <ScalarImageType2D> ( tagged, "segmentedColon.hdr");

	// Dilate
	
	/*
	typedef itk::BinaryDilateImageFilter< ScalarImageType2D, ScalarImageType2D, StructuringElementType > BinaryDilateFilterType;
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

	ScalarImageType2D::Pointer dilate = dilateFilter->GetOutput();
	ByteIteratorType dilateIt( dilate, dilate->GetLargestPossibleRegion() );

	//WriteITK( dilateFilter->GetOutput() , "segmentedColonDilated.hdr");

	// Copy contents to output slice
	for (	outputIt.GoToBegin(), dilateIt.GoToBegin();
			!outputIt.IsAtEnd() && !dilateIt.IsAtEnd();
			++outputIt, ++dilateIt		) 
	{
		outputIt.Set( dilateIt.Get() );
	}
	*/

	// Create binary ball structuring element
	//typedef itk::BinaryBallStructuringElement< ScalarImageType2D::PixelType, 2> StructuringElementType;
	StructuringElementType structuringElement2;
    structuringElement2.SetRadius( 1 );
    structuringElement2.CreateStructuringElement();

	typedef itk::BinaryMorphologicalClosingImageFilter< ScalarImageType2D, ScalarImageType2D, StructuringElementType > BinaryClosingFilterType;
	BinaryClosingFilterType::Pointer closingFilter = BinaryClosingFilterType::New();
	closingFilter->SetInput( tagged );
	closingFilter->SetKernel( structuringElement2 );
	closingFilter->SetForegroundValue( 1 );

	try {
		closingFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType2D::Pointer closed = closingFilter->GetOutput();
	ByteIteratorType closedIt( closed, closed->GetLargestPossibleRegion() );

	//WriteITK <ScalarImageType2D> ( closed, "segmentedColonClosedBallRadius1.hdr");

	// Copy contents to output slice
	for (	outputIt.GoToBegin(), closedIt.GoToBegin();
			!outputIt.IsAtEnd() && !closedIt.IsAtEnd();
			++outputIt, ++closedIt		) 
	{
		outputIt.Set( closedIt.Get() );
	}
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

void WriteITK(ScalarImageType2D::Pointer image, std::string ss) {
	std::cout << "Writing " << ss << std::endl;
	
	typedef itk::ImageFileWriter< ScalarImageType2D > WriterType;
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