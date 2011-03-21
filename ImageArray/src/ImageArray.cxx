#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include <iostream>

typedef itk::Image< int, 3 > ImageType;

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
	ImageType::Pointer image = ImageType::New();

	ImageType::RegionType region;
	ImageType::IndexType start;
	ImageType::SizeType size;

	start[0] = start[1] = start[2] = 0;
	size[0]  = size[1]  = size[2]  = 200;

	region.SetSize( size );
	region.SetIndex( start );

	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(0);

	int * ptr = image->GetBufferPointer();
	ptr[40000] = 1;

	//ImageType::PixelContainerPointer imageArray = image->GetPixelContainer();
	//*imageArray[250] = 1;



	/*

	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	IteratorType it(image,image->GetLargestPossibleRegion());

	it.GoToBegin();

	int count=0;

	for (int i=0; i<400; i++)
	{
		++it;
		std::cout << count++ << std::endl;DAV
	}

	it.Set(1);
	*/

	WriteITK(image,"foo.hdr");

	//system("pause");
	return 0;
}


