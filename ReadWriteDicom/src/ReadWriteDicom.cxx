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
#include <itkNeighborhoodIterator.h>
#include <time.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
 
typedef short													PixelType;
typedef int														ScalarPixelType;
typedef itk::GDCMImageIO										ImageIOType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	

typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;

typedef itk::NeighborhoodIterator< ScalarImageType3D > NeighborhoodIteratorType;

typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;

typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
typedef itk::BinaryMedianImageFilter< ScalarImageType3D, ScalarImageType3D> BinaryMedianFilterType;

typedef itk::ImageSeriesReader< ImageType3D >        ImageSeriesReaderType;
typedef itk::RegularExpressionSeriesFileNames RegexFileNamesType;

template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );
int writeCount = 1;

ImageType3D::Pointer ReadDicom(char path[]);
 
int main( int argc, char* argv[] )
{ 
	// Verify the number of parameters in the command line
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory output\n";
		system("pause");
		return EXIT_FAILURE;
	}
	
	//std::string name = "";
/*
	char * pch;
	pch = strtok( argv[1], "/\\");
	int i = 0;
	while (pch != NULL) {
		printf("%d %s\n", i, pch);
		
		
		if (i == 4) {
			name = pch;
		}
		

		pch = strtok( NULL, "/\\");
		i++;
	}
*/

	//name += "_output.hdr";

	ImageType3D::Pointer input = ReadDicom( argv[1] );
	WriteITK < ImageType3D > ( input, argv[2] );
	
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

ImageType3D::Pointer ReadDicom(char path[])
{	
	// Create reader
	ImageSeriesReaderType::Pointer reader = ImageSeriesReaderType::New();
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	RegexFileNamesType::Pointer fit = RegexFileNamesType::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->NumericSortOn();
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();

    reader->SetFileNames( names );
    reader->Update();

    // Orient all input images into RAS orientation
    itk::OrientImageFilter<ImageType3D,ImageType3D>::Pointer orienter =     itk::OrientImageFilter<ImageType3D,ImageType3D>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI);
    orienter->SetInput(reader->GetOutput());
    orienter->Update();

    return orienter->GetOutput();

}

