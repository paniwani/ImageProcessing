#include "itkImage.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMucosalReconstructionFilter.h>

typedef itk::Image< float, 2 > ImageType;
typedef itk::Image< unsigned char, 2 > ByteImageType;

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
	ImageType::Pointer input = ReadITK("c:/imagedata/Modified_mr10_092_13p.i0344_85_100_40_slice13_subtracted.hdr");

	WriteITK(input,"input.hdr");

	// Make mask
	ByteImageType::Pointer mask = ByteImageType::New();
	mask->SetRegions( input->GetLargestPossibleRegion() );
	mask->SetSpacing( input->GetSpacing() );
	mask->SetOrigin( input->GetOrigin() );
	mask->Allocate();

	typedef itk::ImageRegionIteratorWithIndex< ByteImageType> ByteIteratorType;
	ByteIteratorType mask_iter(mask,input->GetLargestPossibleRegion());

	for (mask_iter.GoToBegin(); !mask_iter.IsAtEnd(); ++mask_iter)
	{	
		ByteImageType::IndexType idx = mask_iter.GetIndex();

		if (idx[0] > 250 && idx[0] < 512-50 && idx[1] > 50 && idx[1] < 512-50)
		{
			mask_iter.Set(1);
		} else {
			mask_iter.Set(0);
		}	
	}

	WriteITK(mask,"mask.hdr");

	typedef itk::MucosalReconstructionFilter<ImageType, ImageType, ByteImageType> MucosalReconstructionFilterType;
	MucosalReconstructionFilterType::Pointer reconstructionFilter = MucosalReconstructionFilterType::New();
	reconstructionFilter->SetInput( input );
	reconstructionFilter->SetMaskImage( mask );
	reconstructionFilter->Update();

	WriteITK(reconstructionFilter->GetOutput(),"reconstruction.hdr");


}