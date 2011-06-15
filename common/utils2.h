#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include <itkMaskImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkConnectedComponentImageFilter.h>

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

struct point{
	int intensity;
	int size;
	int order;
	int border;
	float centroidY;
};

typedef struct point ptype;

typedef short PixelType;

typedef itk::Image<PixelType, Dimension> ImageType;
typedef itk::Image<unsigned char, Dimension> ByteImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Image<int, Dimension> IntImageType;
typedef itk::Image<unsigned int, Dimension> UIntImageType;
typedef itk::Image< VoxelType, Dimension> VoxelImageType;

typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType > ByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex< FloatImageType > FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex< UIntImageType > IntIteratorType;
typedef itk::ImageRegionIteratorWithIndex< VoxelImageType > VoxelIteratorType;
typedef itk::ImageRegionIteratorWithIndex< UIntImageType > UIntIteratorType;

template <typename T>
typename T::Pointer ReadITK(char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< T > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->SetGlobalWarningDisplay(false);

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
	
    return output;
}

template <typename T1, typename T2>
void Mask(typename T1::Pointer &im, typename T2::Pointer &m, typename T1::PixelType outsideValue)
{
	typedef itk::MaskImageFilter<T1,T2> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im);
	masker->SetInput2(m);
	masker->SetOutsideValue(outsideValue);
	masker->Update();
	im = masker->GetOutput();
}

template <typename T1, typename T2>
typename T1::Pointer AllocateImage(typename T2::Pointer &cim)
{
	T1::Pointer im = T1::New();
	im->SetRegions(cim->GetLargestPossibleRegion());
	im->CopyInformation(cim);
	im->Allocate();
	im->FillBuffer(0);
	return im;
}

template <typename T>
void BinaryDilate(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<T, T, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( im );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();
	im = dilater->GetOutput();
}

template <typename T>
void BinaryErode(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter<T, T, StructuringElementType> BinaryErodeImageFilterType;
	BinaryErodeImageFilterType::Pointer eroder = BinaryErodeImageFilterType::New();
	eroder->SetInput( im );
	eroder->SetKernel( se );
	eroder->SetForegroundValue(255);
	eroder->SetBackgroundValue(0);
	eroder->Update();
	im = eroder->GetOutput();
}

template <typename T>
typename T::RegionType BinaryCrop(typename T::Pointer &im, unsigned int pad=5, bool cropZ=false)
{
	T::RegionType region = im->GetLargestPossibleRegion();
	T::SizeType size = region.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0,minZ=0,maxZ=0;
	
	if (Dimension == 3 && cropZ)
	{
		minZ=size[2];
		maxZ=0;
	}

	typedef itk::ImageRegionIteratorWithIndex<T> IteratorType;
	IteratorType it(im,region);

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if ( it.Get() != 0 )
		{
			T::IndexType idx = it.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];

			if (Dimension == 3 && cropZ)
			{
				if (idx[2] < minZ)
					minZ = idx[2];
				if (idx[2] > maxZ)
					maxZ = idx[2];
			}
		}
	}

	T::IndexType edx;
	
	edx[0] = (minX-pad) > 0 ? minX-pad : 0;
	edx[1] = (minY-pad) > 0 ? minY-pad : 0;

	if (Dimension == 3)
	{
		if (cropZ)
		{
			edx[2] = (minZ-pad) > 0 ? minZ-pad : 0;
		} else {
			edx[2] = 0;
		}
	}

	T::SizeType esize;
	esize[0] = maxX-minX+2*pad+1 < size[0] ? maxX-minX+2*pad+1 : size[0];
	esize[1] = maxY-minY+2*pad+1 < size[1] ? maxY-minY+2*pad+1 : size[1];
	
	if (Dimension == 3)
	{
		if (cropZ)
		{
			esize[2] = maxZ-minZ+2*pad+1 < size[2] ? maxZ-minZ+2*pad+1 : size[2];
		} else {
			esize[2] = size[2];
		}
	}

	T::RegionType ror;
	ror.SetIndex( edx );
	ror.SetSize( esize );

	typedef itk::RegionOfInterestImageFilter<T,T> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( im );
	cropper->SetRegionOfInterest( ror );
	cropper->Update();
	im = cropper->GetOutput();

	return ror;
}

template <typename T>
void CropByRegion(typename T::Pointer &im, typename T::RegionType region)
{
	typedef itk::RegionOfInterestImageFilter<T,T> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( im );
	cropper->SetRegionOfInterest( region );
	cropper->Update();
	im = cropper->GetOutput();
}	

template <typename T>
typename T::Pointer Duplicate(typename T::Pointer &im)
{
	typedef itk::ImageDuplicator<T> DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(im);
	duplicator->Update();
	return duplicator->GetOutput();
}