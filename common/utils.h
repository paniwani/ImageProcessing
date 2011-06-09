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
typedef itk::Image<PixelType, 3> ImageType;
typedef itk::Image<PixelType, 2> ImageType2D;
typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<unsigned char, 2> ByteImageType2D;
typedef itk::Image<float, 3> FloatImageType;
typedef itk::Image<float, 2> FloatImageType2D;
typedef itk::Image<int, 3> IntImageType;
typedef itk::Image<int, 2> IntImageType2D;
typedef itk::Image<unsigned int, 3> UIntImageType;
typedef itk::Image<unsigned int, 2> UIntImageType2D;
typedef itk::Image< VoxelType, 3> VoxelImageType;

typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType > ByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex< FloatImageType > FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex< UIntImageType > IntIteratorType;
typedef itk::ImageRegionIteratorWithIndex< VoxelImageType > VoxelIteratorType;
typedef itk::ImageRegionIteratorWithIndex< UIntImageType > UIntIteratorType;


typedef itk::ImageRegionIteratorWithIndex< ImageType2D > IteratorType2D;
typedef itk::ImageRegionIteratorWithIndex< ByteImageType2D > ByteIteratorType2D;
typedef itk::ImageRegionIteratorWithIndex< FloatImageType2D > FloatIteratorType2D;
typedef itk::ImageRegionIteratorWithIndex< IntImageType2D > IntIteratorType2D;
typedef itk::ImageRegionIteratorWithIndex< UIntImageType2D > UIntIteratorType2D;



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

ByteImageType::Pointer AllocateByteImage(ImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(VoxelImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(FloatImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(UIntImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(ByteImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}


FloatImageType::Pointer AllocateFloatImage(ImageType::Pointer &in)
{
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

FloatImageType::Pointer AllocateFloatImage(FloatImageType::Pointer &in)
{
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer BinaryThreshold(ImageType::Pointer &in, PixelType t1, PixelType t2)
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	IteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() > t1 && it.Get() < t2)
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(ImageType::Pointer &in, PixelType t, bool reverse=false)
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	IteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() < t )
		{
			if (reverse)
			{
				oit.Set(0);
			} else {
				oit.Set(255);
			}
		} else {
			if (reverse)
			{
				oit.Set(255);
			} else {
				oit.Set(0);
			}
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(UIntImageType::Pointer &in, unsigned int t, bool reverse=false)
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	UIntIteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() < t )
		{
			if (reverse)
			{
				oit.Set(0);
			} else {
				oit.Set(255);
			}
		} else {
			if (reverse)
			{
				oit.Set(255);
			} else {
				oit.Set(0);
			}
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(FloatImageType::Pointer &in, float t1=itk::NumericTraits<float>::NonpositiveMin(), float t2=itk::NumericTraits<float>::max())
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	FloatIteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() > t1 && it.Get() < t2)
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(VoxelImageType::Pointer &vmap, VoxelType v)
{
	ByteImageType::Pointer out = AllocateByteImage(vmap);
	
	VoxelIteratorType it(vmap,vmap->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() == v )
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer Mask(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<ByteImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

FloatImageType::Pointer Mask(FloatImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<FloatImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

ImageType::Pointer Mask(ImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<ImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

//ByteImageType::Pointer BinaryKeeper(ByteImageType::Pointer &in, std::string attribute, unsigned int numOfObjects, bool reverse=false)
//{
//	typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType> KeeperType;
//	KeeperType::Pointer keeper = KeeperType::New();
//	keeper->SetInput(in);
//	keeper->SetAttribute(attribute);
//	keeper->SetNumberOfObjects(numOfObjects);
//	keeper->SetReverseOrdering(reverse);
//	keeper->SetForegroundValue(255);
//	keeper->SetBackgroundValue(0);
//	keeper->Update();
//	return keeper->GetOutput();
//}

ByteImageType::Pointer BinaryOr(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::OrImageFilter<ByteImageType> OrType;
	OrType::Pointer or = OrType::New();
	or->SetInput1(im1);
	or->SetInput2(im2);
	or->Update();
	return or->GetOutput();
}

ByteImageType::Pointer BinarySubtract(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{
	typedef itk::SubtractImageFilter<ByteImageType> SubtractType;
	SubtractType::Pointer subtracter = SubtractType::New();
	subtracter->SetInput1(im1);
	subtracter->SetInput2(im2);
	subtracter->Update();
	return subtracter->GetOutput();
}

template <typename T1, typename T2>
typename T2::Pointer Cast(typename T1::Pointer &in)
{
	typedef itk::CastImageFilter<T1,T2> CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput(in);
	caster->Update();
	return caster->GetOutput();
}

void Rescale(FloatImageType::Pointer &in, float min, float max)
{
	typedef itk::RescaleIntensityImageFilter<FloatImageType> RescaleType;
	RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput(in);
	rescaler->SetOutputMaximum(max);
	rescaler->SetOutputMinimum(min);
	rescaler->Update();
	in = rescaler->GetOutput();
}

ImageType::Pointer Crop(ImageType::Pointer &input, ImageType::RegionType requestedRegion)
{
	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> CropType;
	CropType::Pointer cropper = CropType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(requestedRegion);
	cropper->Update();
	return cropper->GetOutput();
}

void GetMinMax(ImageType::Pointer &im, PixelType &min, PixelType &max)
{
	typedef itk::MinimumMaximumImageCalculator<ImageType> CalcType;
	CalcType::Pointer calc = CalcType::New();
	calc->SetImage(im);
	calc->Compute();
	min = calc->GetMinimum();
	max = calc->GetMaximum();
}

void GetMinMax(UIntImageType::Pointer &im, unsigned int &min, unsigned int &max)
{
	typedef itk::MinimumMaximumImageCalculator<UIntImageType> CalcType;
	CalcType::Pointer calc = CalcType::New();
	calc->SetImage(im);
	calc->Compute();
	min = calc->GetMinimum();
	max = calc->GetMaximum();
}

void GetMinMax(ByteImageType::Pointer &im, unsigned char &min, unsigned char &max)
{
	typedef itk::MinimumMaximumImageCalculator<ByteImageType> CalcType;
	CalcType::Pointer calc = CalcType::New();
	calc->SetImage(im);
	calc->Compute();
	min = calc->GetMinimum();
	max = calc->GetMaximum();
}

UIntImageType::Pointer CC(ByteImageType::Pointer &im, bool fullConnect=false)
{
	typedef itk::ConnectedComponentImageFilter<ByteImageType,UIntImageType> CCFilterType;
	CCFilterType::Pointer ccer = CCFilterType::New();
	ccer->SetInput(im);
	ccer->SetFullyConnected(fullConnect);
	ccer->Update();
	return ccer->GetOutput();
}

bool compare(unsigned int i,unsigned int j) { return (i<j); }

bool compareSize(ptype a, ptype b)
{
	if(a.size > b.size)
	{
		return true;
	}
	else return false;
}

bool compareIntensity(ptype a, ptype b)
{
	if(a.intensity < b.intensity)
	{
		return true;
	}
	else return false;
}

bool compareBorder(ptype a, ptype b)
{
	if(a.border > b.border)
	{
		return true;
	}
	else return false;
}

bool compareCentroidY(ptype a, ptype b)
{
	if(a.centroidY > b.centroidY)
	{
		return true;
	}
	else return false;
}

void RelabelSize(UIntImageType::Pointer &image, unsigned int minSize = 0)
{
	UIntIteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		a.border = 0;
		a.centroidY = 0;
		count.push_back(a);
	}
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			++count[it.Get() - 1].size;
		}
	}
	

	sort(count.begin(), count.end(), compareSize);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareIntensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get() - 1].size > minSize)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}
}

void RelabelCentroidY(UIntImageType::Pointer &image, unsigned int minCentroidY = 0)
{
	UIntIteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		a.border = 0;
		a.centroidY = 0;
		count.push_back(a);
	}
	
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			ByteImageType::IndexType idx = it.GetIndex();
			count[it.Get() - 1].centroidY += idx[1];
			++count[it.Get() - 1].size;
		}
	}

	for(int i = 0; i < count.size(); i++)
	{
		count[i].centroidY /= (float) count[i].size;
	}

	sort(count.begin(), count.end(), compareCentroidY);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareIntensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get() - 1].centroidY > minCentroidY)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}
}

void RelabelCentroidY(ByteImageType::Pointer &image, unsigned int minCentroidY = 0)
{
	ByteIteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		a.border = 0;
		a.centroidY = 0;
		count.push_back(a);
	}
	
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			ByteImageType::IndexType idx = it.GetIndex();
			count[it.Get() - 1].centroidY += idx[1];
			++count[it.Get() - 1].size;
		}
	}

	for(int i = 0; i < count.size(); i++)
	{
		count[i].centroidY /= (float) count[i].size;
	}

	sort(count.begin(), count.end(), compareCentroidY);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareIntensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get() - 1].centroidY > minCentroidY)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}
}



void RelabelBorder(UIntImageType::Pointer &image, unsigned int minBorder = 0)
{
	UIntIteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		a.border = 0;
		a.centroidY = 0;
		count.push_back(a);
	}

	UIntImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			ByteImageType::IndexType idx = it.GetIndex();
			if (idx[0] == 0 || idx[0] == size[0] || idx[1] == 0 || idx[1] == size[1])
				++count[it.Get() - 1].border;
		}
	}
	

	sort(count.begin(), count.end(), compareBorder);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareIntensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get() - 1].border > minBorder)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}
}

void RelabelBorder(ByteImageType::Pointer &image, unsigned int minBorder = 0)
{
	ByteIteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		a.border = 0;
		a.centroidY = 0;
		count.push_back(a);
	}

	ByteImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			ByteImageType::IndexType idx = it.GetIndex();
			if (idx[0] == 0 || idx[0] == size[0] || idx[1] == 0 || idx[1] == size[1])
				++count[it.Get() - 1].border;
		}
	}
	

	sort(count.begin(), count.end(), compareBorder);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareIntensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get() - 1].border > minBorder)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}
}


void Threshold(UIntImageType::Pointer &im, unsigned int a, unsigned int b)
{
	UIntIteratorType it(im,im->GetLargestPossibleRegion());

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() < a || it.Get() > b)
		{
			it.Set(0);
		}
	}
}

ByteImageType::Pointer Binarize(UIntImageType::Pointer &im)
{
	ByteImageType::Pointer output = AllocateByteImage(im);
	ByteIteratorType oit(output,output->GetLargestPossibleRegion());

	UIntIteratorType it(im,im->GetLargestPossibleRegion());
	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if (it.Get() > 0)
		{
			oit.Set(255);
		} else {
			oit.Set(0);
		} 
	}

	return output;
}

ByteImageType::Pointer Binarize(ByteImageType::Pointer &im)
{
	ByteImageType::Pointer output = AllocateByteImage(im);
	ByteIteratorType oit(output,output->GetLargestPossibleRegion());

	ByteIteratorType it(im,im->GetLargestPossibleRegion());
	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if (it.Get() > 0)
		{
			oit.Set(255);
		} else {
			oit.Set(0);
		} 
	}

	return output;
}