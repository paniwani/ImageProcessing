#include <itkImage.h>
#include <itkFastDilateFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

typedef itk::Image<float, 2> ImageType;
typedef itk::Image<unsigned char, 2> ByteImageType;

typedef itk::ImageRegionIterator<ImageType> IteratorType;
typedef itk::ImageRegionIterator<ByteImageType> ByteIteratorType;

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

	image->SetDirection( direction );
	
	return image;
}

void WriteITK(ImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

int main()
{
	ImageType::Pointer input = ReadITK("C:/ImageData/Modified_mr10_092_13p.i0344_85_100_40_slice13.hdr");
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(input->GetSpacing());
	air->Allocate();
	ByteIteratorType air_iter(air,region);
	IteratorType input_iter(input,region);
	
	for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter)
	{
		if (input_iter.Get() < -600)
		{
			air_iter.Set(1);
		} else {
			air_iter.Set(0);
		}
	}
	
	WriteITK(air,"air.hdr");
	
	typedef itk::FastDilateFilter<ByteImageType> FastDilateFilterType;
	FastDilateFilterType::Pointer dilater = FastDilateFilterType::New();
	dilater->SetInput(air);
	dilater->Update();
	
	WriteITK(dilater->GetOutput(),"dilate.hdr");
	
	return 0;
}