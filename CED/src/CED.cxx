#include "itkImage.h"
#include <iostream>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include "itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkAnisotropicHybridDiffusionImageFilter.h"
#include "itkAnisotropicEdgeEnhancementDiffusionImageFilter.h"

bool truncateOn = true;
float truncate_ar[2] = {85,100};

typedef itk::Image< float, 3 > ImageType;

typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorTypeFloat4WithIndex;

ImageType::Pointer ReadDicom( std::string path );
void WriteITK(ImageType::Pointer image, std::string name);
ImageType::Pointer ReadITK(char * fileName);

int main()
{
	
	ImageType::Pointer input = ReadDicom( "C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm" );
	//ImageType::Pointer input = ReadITK("Modified_mr10_092_13p.i0344_85_100_10_hessian_pre_process.hdr");

	WriteITK(input, "input.hdr");

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	typedef itk::AnisotropicEdgeEnhancementDiffusionImageFilter<ImageType,ImageType> DiffusionFilterType;
	DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
	diffusionFilter->SetInput( input );
	diffusionFilter->Update();
	ImageType::Pointer output = diffusionFilter->GetOutput();

	WriteITK(output, "output.hdr");
/*
	// Write input and output greater than 0
	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		if (input_iter.Get() > 0)
		{
			input_iter.Set( 1 );
		} else {
			input_iter.Set( 0 );
		}
	}

	for (output_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter)
	{
		if (output_iter.Get() > 0)
		{
			output_iter.Set( 1 );
		} else {
			output_iter.Set( 0 );
		}
	}

	WriteITK(input,"input_thres.hdr");
	WriteITK(output,"output_thres.hdr");
*/

	//system("pause");
	return 0;
}

ImageType::Pointer ReadDicom( std::string path )
{	
	// Create reader
	itk::ImageSeriesReader<ImageType>::Pointer reader = itk::ImageSeriesReader<ImageType>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	if (truncateOn)
	{
		names.erase( names.begin(), names.begin()+truncate_ar[0] );
		names.erase( names.begin()+truncate_ar[1]-truncate_ar[0], names.end() );
	}

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return 0;
	}

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput(reader->GetOutput());
    orienter->Update();

	// Set direction cosines to identity matrix
	ImageType::Pointer output = orienter->GetOutput();
	
	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	direction[2][2] = 1;
	
	output->SetDirection(direction);

    return output;

}

void WriteITK(ImageType::Pointer image, std::string name) 
{
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
	writer.~SmartPointer();
}

ImageType::Pointer ReadITK(char * fileName) 
{
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();

	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	direction[2][2] = 1;
	
	image->SetDirection(direction);

	return image;
}