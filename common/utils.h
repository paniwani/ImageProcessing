#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"

typedef short PixelType;
typedef itk::Image<PixelType, 3> ImageType;
typedef itk::Image<PixelType, 2> ImageType2D;
typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<unsigned char, 2> ByteImageType2D;

typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType2D > IteratorType2D;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType > ByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType2D > ByteIteratorType2D;

template <typename T>
typename T::Pointer ReadITK(char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< T > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );

	try {
		reader->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error reading image: " << err << std::endl;
		return NULL;
	}

	return reader->GetOutput();
}

template <typename T>
void WriteITK(typename T::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< T >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	
	try {
		writer->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error writing image: " << err << std::endl;
	}
}

template <typename T>
typename T::Pointer ReadDicom( std::string path, int slice1=0, int slice2=-1)
{	
	// Create reader
	itk::ImageSeriesReader<T>::Pointer reader = itk::ImageSeriesReader<T>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	/*
	
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	
	nameGenerator->SetDirectory( path );
	
	const std::vector< std::string > seriesUID = nameGenerator->GetSeriesUIDs();
	std::string seriesIdentifier = seriesUID.begin()->c_str();
	
	std::vector< std::string > names = nameGenerator->GetFileNames( seriesIdentifier );

	*/
	
	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	if (T::ImageDimension == 2 && slice2 == -1)
	{
		names.erase( names.begin(), names.begin()+slice1-1 );
		names.erase( names.begin()+1, names.end() );
	}
	
	if (slice2 > 0 && slice2 > slice1)
	{
		names.erase( names.begin(), names.begin()+slice1);
		names.erase( names.begin()+slice2-slice1, names.end() );
	}
	
    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error reading dicom: " << err << std::endl;
		return 0;
	}
	
	T::Pointer output = reader->GetOutput();
	
	/*
    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<T,T>::Pointer orienter = itk::OrientImageFilter<T,T>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput( output );
    orienter->Update();
	output = orienter->GetOutput();
	*/
	
	
    return output;
}