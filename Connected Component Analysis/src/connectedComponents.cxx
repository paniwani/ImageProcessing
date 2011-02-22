#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

struct point{
	int intensity;
	int size;
	int order;
};

typedef struct point ptype;

bool compare(ptype a, ptype b)
{
	if(a.size > b.size)
	{
		return true;
	}
	else return false;
}

bool compareintensity(ptype a, ptype b)
{
	if(a.intensity < b.intensity)
	{
		return true;
	}
	else return false;
}

int main(int argc, char *argv[])
{
	// itk definitions
	typedef unsigned short PixelType;
	typedef itk::Image< short, 3 > ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
	
	// declare and set up input and output filters
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	ConnectedComponentType::Pointer ccFilter = ConnectedComponentType::New();
	ccFilter->SetInput(reader->GetOutput());

	try
	{
		ccFilter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}
	
	ImageType::Pointer image = ccFilter->GetOutput();

	IteratorType it(image, image->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		//short itvalue = it.Get();
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		std::cout << "working on component " << i << std::endl;
		ptype a;
		a.intensity = i;
		a.size = 0;
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if(it.Get() == i)
			{
				++(a.size);
			}
		}
		count.push_back(a);
	}
	

	sort(count.begin(), count.end(), compare);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
		//std::cout << count[i].intensity << std::endl;
	}

	sort(count.begin(), count.end(), compareintensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			it.Set(count[it.Get() - 1].order);
		}
	}




	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(argv[2]);
	


	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

system("pause");
return 0;
}
