#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "windows.h"
#include <iostream>
#include <itkCastImageFilter.h>

typedef itk::Image<unsigned short, 2> ImageType;
typedef itk::Image<unsigned int, 2> ImageType2;
typedef itk::ImageRegionIterator<ImageType> IteratorType;


template <class ObjectType>
void printRefCount(ObjectType a, char* message)
{
	std::cout << message << std::endl;
	std::cout << "Reference Count is: " << a->GetReferenceCount() - 1 << std::endl;
}

template <class rImageType>
typename rImageType::Pointer readImage(char *filename)
{
	typedef itk::ImageFileReader<rImageType> ReaderType;
	// read input labeled skeleton
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << err << std::endl;
		return (ImageType::Pointer)NULL;
	}
	printRefCount<ImageType::Pointer>(reader->GetOutput(), "Image in reader function");
	return reader->GetOutput();
}

template <class wImageType>
void writeImage(typename wImageType::Pointer image, char* filename)
{
	typedef itk::ImageFileWriter<wImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	printRefCount<wImageType::Pointer>(image, "Image writer before setting input");

	writer->SetInput(image);
	printRefCount<wImageType::Pointer>(image, "Image writer after setting input");
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << err << std::endl;
	}
}


int main(int argc, char* argv [] )
{
	// itk uses reference counting for memory allocation, when an object is created, it has an rc of 1, when another
	// object (i.e. a filter) needs it, it increments the rc to 2, when that filter is destroyed, it decrements the rc
	// back to one. When/if the reference count reaches 0 the image is deleted automatically


	// Testing of reading and writing functions
	ImageType::Pointer image = readImage<ImageType>(argv[1]);
	
	//printRefCount<ImageType::Pointer>(image, "After read function return");
	
	//writeImage<ImageType>(image, argv[2]);

	ImageType::IndexType idx;
	idx[0] = 50;
	idx[1] = 50;

	//std::cout << image->GetPixel( idx ) << std::endl;
	
	//printRefCount<ImageType::Pointer>(image, "After write function return");


	typedef itk::CastImageFilter<ImageType, ImageType2> CastFilterType;
	typedef itk::CastImageFilter<ImageType2, ImageType> CastFilterType2;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	CastFilterType2::Pointer castFilter2 = CastFilterType2::New();

	printRefCount<ImageType::Pointer>(image, "Before cast filter one");
	
	castFilter->SetInput( image );
	//castFilter->SetReleaseDataFlag(true);

	printRefCount<ImageType::Pointer>(image, "After cast filter one");

	//castFilter2->SetInput( castFilter->GetOutput() );
	//castFilter2->Update();

	//printRefCount<ImageType::Pointer>(image, "After cast filter two");

	//ImageType::Pointer foo = NULL;

	ImageType2::Pointer imageCast = castFilter->GetOutput();

	printRefCount<ImageType::Pointer>(image, "After creating image cast pointer");

	castFilter.~SmartPointer();

	printRefCount<ImageType::Pointer>(image, "After deleting filter");

	imageCast.~SmartPointer();

	castFilter = CastFilterType::New();







	

	/*

	ImageType2::IndexType idx2;
	idx2[0] = 50;
	idx2[1] = 50;

	std::cout << image->GetPixel( idx2 ) << std::endl;

	writeImage < ImageType2 > (imageCast, "imageCast.hdr" );


	ImageType2::IndexType idx;
	idx[0] = 4;
	idx[1] = 4;
	imageCast->GetPixel( idx );
	castFilter2->GetOutput();

	*/

	
	
	
	// Assignment Experiment
	// assignment to a smart pointer does increase the reference count, 
	//thus if you pass an ImageType::Pointer into a function, the reference count will increase
	// by one in that function, but return to normal after the function returns
	//ImageType::Pointer image2 = image;
	//printRefCount<ImageType::Pointer>(image, "After assignment to another smart pointer");

	//image->SetReferenceCount( 0 );
	//image = NULL;

	image.~SmartPointer();

	//image2->SetReferenceCount( 0 );
	//image2 = NULL;


	/* if you were to do the following, i.e. decrease the reference count to 0, the image data would be deleted
	//any further calls to member functions of image or image2 would cause errors as the image data was freed
	//In addition, the program will crash at the end of the main function since the pointers image and image2 will be
	//destructed, at which point they will attempt to further decrease the reference count of the image object, which has already been deleted
	Sleep(10000); // observe process memory consumption before image deletion
	image->UnRegister(); // decrements an objects reference count
	image->UnRegister();
	Sleep(10000); // confirm image data has been freed
	*/

	// iterator experiments (Questions of Interest: is the iterator object destructable? yes, but you must dynamically allocate it
	//, and does it affect an images reference count? no)
	
	/*
	IteratorType it(image, image->GetLargestPossibleRegion());
	printRefCount<ImageType::Pointer>(image, "After construction of iterator for image");
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
	}
	
	// dynamically allocating iterator
	IteratorType * ittest = new IteratorType(image, image->GetLargestPossibleRegion());
	printRefCount<ImageType::Pointer>(image, "After construction of dynamic iterator for image");
	for(ittest->GoToBegin(); !ittest->IsAtEnd(); ++(*ittest))
	{
	}
	delete ittest;

	// replacing ittest with new iterator
	ittest = new IteratorType(image2, image2->GetLargestPossibleRegion());
	for(ittest->GoToBegin(); !ittest->IsAtEnd(); ++(*ittest))
	{
	}
	delete ittest;

	*/

	
	system("pause");
	return 0;
	
	
}