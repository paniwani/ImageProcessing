#include "itkImage.h"
#include <iostream>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkFastChamferDistanceImageFilter.h>


typedef float                                                               PixelType;
typedef itk::Image<PixelType, 3>                                            ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType>                        IteratorTypeFloat4WithIndex;

typedef itk::FastChamferDistanceImageFilter<ImageType,ImageType>			FastChamferFilterType;


void ReadITK(ImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	image = reader->GetOutput();
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

int main()
{
	ImageType::Pointer input = ImageType::New();
	ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_100-105.hdr");

	IteratorTypeFloat4WithIndex input_iter(input,input->GetLargestPossibleRegion());

	ImageType::Pointer th = ImageType::New();
	th->SetRegions(input->GetLargestPossibleRegion());
	th->Allocate();
	IteratorTypeFloat4WithIndex th_iter(th,input->GetLargestPossibleRegion());

	input_iter.GoToBegin();
	th_iter.GoToBegin();

	while (!input_iter.IsAtEnd())
	{

		if (input_iter.Get() > 200) {
			th_iter.Set(0);
		} else {
			th_iter.Set(1);
		}
	}
	
	FastChamferFilterType::Pointer fc = FastChamferFilterType::New();
	fc->SetInput(th);
	fc->Update();
	WriteITK(fc->GetOutput(), "fc.hdr");
	
	return 0;
}

