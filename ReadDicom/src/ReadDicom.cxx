#include <itkImage.h> 						
#include <iostream> 							
 								
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkDICOMSeriesFileNames.h"
 
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include <utils.h>

int main(int argc, char * argv[])				
{ 			
	const unsigned int InputDimension = 3;
	const unsigned int OutputDimension = 2;

	typedef signed short PixelType;

	typedef itk::Image< PixelType, InputDimension >
	InputImageType;
	typedef itk::Image< PixelType, OutputDimension >
	OutputImageType;
	typedef itk::ImageSeriesReader< InputImageType >
	ReaderType;
	typedef itk::GDCMImageIO
	ImageIOType;
	//typedef itk::GDCMSeriesFileNames
	typedef itk::DICOMSeriesFileNames
	InputNamesGeneratorType;
	typedef itk::NumericSeriesFileNames
	OutputNamesGeneratorType;
	typedef itk::ImageSeriesWriter< InputImageType, OutputImageType >
	SeriesWriterType;

	////////////////////////////////////////////////  
	// 1) Read the input series

	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	InputNamesGeneratorType::Pointer inputNames = InputNamesGeneratorType::New();
	//inputNames->SetGlobalWarningDisplay(false);
	inputNames->SetDirectory( "C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm" );
	//inputNames->SetSeriesFormat("%d");

	const ReaderType::FileNamesContainer & filenames = 
							inputNames->GetFileNames();

	ReaderType::Pointer reader = ReaderType::New();

	reader->SetImageIO( gdcmIO );
	reader->SetFileNames( filenames );
	try
	{
	reader->Update();
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while reading the series" << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	WriteITK <ImageType> (reader->GetOutput(),"out.nii");


	system("pause"); 							
	return 0; 									
} 												
