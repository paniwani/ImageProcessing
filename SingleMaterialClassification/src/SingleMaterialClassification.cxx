#include "itkImage.h"
#include <iostream>
#include <itkColonSegmentationFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"

enum VoxelType {
	Stool=1,
	Air=2,
	Tissue=3,
	Unclassified=4,
	StoolAir=5,
	TissueAir=6,
	TissueStool=7,
    ThinStool=8
};

typedef itk::Image< float, 3 > ImageType;
typedef itk::Image< unsigned char, 3 > ByteImageType;
typedef itk::Image< VoxelType, 3 > VoxelImageType;

typedef itk::ImageRegionIterator< ImageType > IteratorType;
typedef itk::ImageRegionIterator< ByteImageType > ByteIteratorType;
typedef itk::ImageRegionIterator< VoxelImageType> VoxelIteratorType;

std::string note="";


ImageType::Pointer ReadITK(char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );

	try {
		reader->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}

	ImageType::Pointer image = reader->GetOutput();

	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	
	if ( image->GetImageDimension() == 3 )
		direction[2][2] = 1;

	image->SetDirection( direction );
	
	return image;
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	name = note + "_" + name;
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(ImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	name = note + "_" + name;
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(VoxelImageType::Pointer vimage, std::string name) {
	// Get voxel type iter
	VoxelIteratorType vit( vimage, vimage->GetLargestPossibleRegion() );
	
	// Create helper image
	ByteImageType::Pointer temp = ByteImageType::New();
	temp->SetRegions( vimage->GetLargestPossibleRegion() );
	temp->Allocate();
	ByteIteratorType tit( temp, temp->GetLargestPossibleRegion() );

	for (vit.GoToBegin(), tit.GoToBegin(); !vit.IsAtEnd() && !tit.IsAtEnd(); ++vit, ++tit)
	{
		 switch(vit.Get()) {
            case Stool:
                tit.Set(1);
                break;
			case Air:
				tit.Set(2);
				break;
			case Tissue:
				tit.Set(3);
				break;
			case Unclassified:
				tit.Set(4);
				break;
			case StoolAir:
				tit.Set(5);
				break;
			case TissueAir:
				tit.Set(6);
				break;
			case TissueStool:
				tit.Set(7);
				break;
			case ThinStool:
				tit.Set(8);
				break;
			default:
                tit.Set(0);
                break;
        }

		tit.Set( static_cast<unsigned char> (tit.Get() * 255/8) );
	}

	WriteITK(temp, name);

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
	
	/*
	unsigned int firstSlice = 0;
	unsigned int lastSlice = 150;

	names.erase( names.begin(), names.begin()+firstSlice);
	names.erase( names.begin()+lastSlice-firstSlice, names.end() );
	*/

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

	ImageType::Pointer output = reader->GetOutput();

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput( output );
    orienter->Update();
	output = orienter->GetOutput();

	/*
	// Set direction cosines to identity matrix	
	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	direction[2][2] = 1;
	
	output->SetDirection(direction);
	*/

    return output;

}

std::vector<std::string> explode( const std::string &delimiter, const std::string &str)
{
	std::vector<std::string> arr;

    int strleng = str.length();
    int delleng = delimiter.length();
    if (delleng==0)
        return arr;//no change

    int i=0; 
    int k=0;
    while( i<strleng )
    {
        int j=0;
        while (i+j<strleng && j<delleng && str[i+j]==delimiter[j])
            j++;
        if (j==delleng)//found delimiter
        {
            arr.push_back(  str.substr(k, i-k) );
            i+=delleng;
            k=i;
        }
        else
        {
            i++;
        }
    }
    arr.push_back(  str.substr(k, i-k) );
    return arr;
}

int main(int argc, char * argv[])
{
	// Get input
	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}
	std::cerr<<"Started"<<std::endl;

	std::string dataset = argv[1];

	std::vector<std::string> datasetArr = explode( "/", dataset );
	std::string dsname = datasetArr[ datasetArr.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	// Set writer prefix
	note = dsname;
	
	// Load dicom
	ImageType::Pointer input = ReadDicom( dataset );

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	//WriteITK(input,"input.hdr");
	
	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilter;
	ColonSegmentationFilter::Pointer segmenter = ColonSegmentationFilter::New();
	segmenter->SetInput( input );
	segmenter->Update();
	ByteImageType::Pointer colon = segmenter->GetOutput();
	ByteIteratorType colon_iter(colon,region);

	WriteITK(colon,"colon_segmentation.hdr");
	
	// Compute gradient
	typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gradientMagFilter = GradientMagnitudeImageFilterType::New();
	gradientMagFilter->SetInput(input);
	gradientMagFilter->Update();
	ImageType::Pointer gmag = gradientMagFilter->GetOutput();
	IteratorType gmag_iter(gmag,region);
	
	// Apply single material thresholds
	VoxelImageType::Pointer vmap = VoxelImageType::New();
	vmap->SetRegions(region);
	vmap->SetSpacing(input->GetSpacing());
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);
	VoxelIteratorType vmap_iter(vmap,region);
	
	IteratorType input_iter(input,region);
	
	for (input_iter.GoToBegin(), gmag_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gmag_iter, ++vmap_iter, ++colon_iter)
	{
		if (colon_iter.Get() == 1)
		{
			VoxelType voxel = Unclassified;
			float input_pixel = input_iter.Get();
			float input_gradient_pixel = gmag_iter.Get();
		
			if ((input_pixel >=180 && input_gradient_pixel<=0.8*input_pixel)) {
				voxel = Stool;
			} else if (input_pixel<=-800 && input_gradient_pixel<=250) {
				voxel = Air;
			} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
				voxel = Tissue;
			}
			
			vmap_iter.Set(voxel);
		}
	}
	
	WriteITK(vmap,"vmap_carston.hdr");
	
	vmap.~SmartPointer();
	
	vmap = VoxelImageType::New();
	vmap->SetRegions(region);
	vmap->SetSpacing(input->GetSpacing());
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);
	vmap_iter = VoxelIteratorType(vmap,region);
	
	for (input_iter.GoToBegin(), gmag_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gmag_iter, ++vmap_iter, ++colon_iter)
	{
		if (colon_iter.Get() == 1)
		{
			VoxelType voxel = Unclassified;
			float input_pixel = input_iter.Get();
			float input_gradient_pixel = gmag_iter.Get();
		
			if ((input_pixel >=400 && input_gradient_pixel<=0.8*input_pixel) || input_pixel >= 1000) {
				voxel = Stool;
			} else if (input_pixel<=-700) {
				voxel = Air;
			} else if (input_pixel<=400  && input_pixel>=-250 && input_gradient_pixel<=400) {
				voxel = Tissue;
			}
			
			vmap_iter.Set(voxel);
		}
	}
	
	WriteITK(vmap,"vmap_modified.hdr");
	
	return 0;
}