/*
Compute all local textures from OTB and
average across offsets
*/

#include "itkImage.h"
#include <iostream>
#include <otbScalarImageToTexturesFilter.h>
#include <otbScalarImageToAdvancedTexturesFilter.h>
#include <otbScalarImageToHigherOrderTexturesFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNeighborhood.h>
#include <itkDivideByConstantImageFilter.h>
#include <itkAddImageFilter.h>

typedef float PixelType;
typedef itk::Image<PixelType,2> ImageType2D;
typedef itk::Image<unsigned char,2> ByteImageType2D;
typedef ImageType2D::OffsetType OffsetType;

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

int main()
{
	// load input
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("att.png");

	// set radius
	const unsigned int radius = 2;
	ImageType2D::SizeType rad;
	rad.Fill(radius);

	// allocate 28 outputs (8 haralick, 9 advanced, 11 higher order run length)
	const unsigned int numOutputs = 28;

	ImageType2D::Pointer outputArray[28];
	for (int i=0; i<numOutputs; i++)
	{
		outputArray[i] = ImageType2D::New();
		outputArray[i]->SetRegions(input->GetLargestPossibleRegion());
		outputArray[i]->CopyInformation(input);
		outputArray[i]->Allocate();
		outputArray[i]->FillBuffer(0);
	}

	// get all offsets
	typedef itk::Neighborhood<PixelType,2> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(1);

	std::vector<OffsetType> offsetVector;

	for (int i=0; i<(hood.Size()-1)/2; i++)
	{
		OffsetType of = hood.GetOffset(i);
		offsetVector.push_back(of);
	}

	unsigned int numOffsets = offsetVector.size();

	std::cout << "Number of outputs: " << numOutputs << std::endl;
	std::cout << "Number of offsets: " << numOffsets << std::endl;
	std::cout << "Radius: " << radius << std::endl;

	typedef otb::ScalarImageToTexturesFilter<ImageType2D,ImageType2D> ScalarImageToTexturesFilterType;
	typedef otb::ScalarImageToHigherOrderTexturesFilter<ImageType2D,ImageType2D> ScalarImageToHigherOrderTexturesFilterType;
	typedef otb::ScalarImageToAdvancedTexturesFilter<ImageType2D,ImageType2D> ScalarImageToAdvancedTexturesFilterType;
	
	for (int i=0; i<numOffsets; i++)
	{
		// compute haralick textures
		ScalarImageToTexturesFilterType::Pointer haralickFilter = ScalarImageToTexturesFilterType::New();
		haralickFilter->SetInput(input);
		haralickFilter->SetInputImageMinimum(0);
		haralickFilter->SetInputImageMaximum(255);
		haralickFilter->SetRadius(rad);
		haralickFilter->SetOffset(offsetVector[i]);
		haralickFilter->Update();

		// compute advanced textures
		ScalarImageToAdvancedTexturesFilterType::Pointer advancedFilter = ScalarImageToAdvancedTexturesFilterType::New();
		advancedFilter->SetInput(input);
		advancedFilter->SetInputImageMinimum(0);
		advancedFilter->SetInputImageMaximum(255);
		advancedFilter->SetRadius(rad);
		advancedFilter->SetOffset(offsetVector[i]);
		advancedFilter->Update();

		// compute run length textures
		ScalarImageToHigherOrderTexturesFilterType::Pointer runFilter = ScalarImageToHigherOrderTexturesFilterType::New();
		runFilter->SetInput(input);
		runFilter->SetInputImageMinimum(0);
		runFilter->SetInputImageMaximum(255);
		runFilter->SetRadius(rad);
		runFilter->SetOffset(offsetVector[i]);
		runFilter->Update();

		// sum texture outputs over all offsets
		for (int j=0; j<numOutputs; j++)
		{
			typedef itk::AddImageFilter<ImageType2D,ImageType2D> AddImageFilterType;
			AddImageFilterType::Pointer adder = AddImageFilterType::New();
			adder->SetInput1(outputArray[j]);
			
			if (j < 8)
			{
				adder->SetInput2(haralickFilter->GetOutput(j));
			} else if ( j >= 8 && j < 17) {
				adder->SetInput2(advancedFilter->GetOutput(j-8));
			} else {
				adder->SetInput2(runFilter->GetOutput(j-17));
			}

			adder->Update();
			outputArray[j] = adder->GetOutput();	
		}

		std::cout << "Completed offset number " << i+1 << std::endl;
	}

	// get average texture
	for (int j=0; j<numOutputs; j++)
	{
		typedef itk::DivideByConstantImageFilter<ImageType2D,PixelType,ImageType2D> DivideByConstantImageFilterType;
		DivideByConstantImageFilterType::Pointer divider = DivideByConstantImageFilterType::New();
		divider->SetInput(outputArray[j]);
		divider->SetConstant(numOffsets);
		divider->Update();
		outputArray[j] = divider->GetOutput();

		typedef itk::RescaleIntensityImageFilter<ImageType2D,ByteImageType2D> RescaleToByteType;
		RescaleToByteType::Pointer rbFilter = RescaleToByteType::New();
		rbFilter->SetInput(outputArray[j]);
		rbFilter->SetOutputMaximum(255);
		rbFilter->SetOutputMinimum(0);
		rbFilter->Update();

		std::stringstream ss;
		ss << "texture" << j << ".png";
		WriteITK <ByteImageType2D> (rbFilter->GetOutput(),ss.str());

	}

	system("pause");

	return 0;
}
