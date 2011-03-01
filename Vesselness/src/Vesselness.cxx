#include <string>
#include <fstream>

#include <iostream>
#include <sstream>

#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <stdlib.h>
//#include <stdio.h>



#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"

float vesselness(const float lambda1, const float lambda2, const float lambda3);
template <class wImageType>
void writeImage(typename wImageType::Pointer image, char* filename)
{
	typedef itk::ImageFileWriter<wImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << err << std::endl;
	}
}

std::vector<int> buildCumHistogramBottom(std::vector<int> histogramData)
{
	for(int i = histogramData.size() - 1; i >= 0; i--)
	{
		if(i < (histogramData.size() - 1))
		{
			histogramData[i] += histogramData[i + 1];
		}
	}
	return histogramData;
}

int main(int argc, char **argv)
{ 

	// itk pixel type defintions
	typedef short InputPixelType;
	typedef float OutputPixelType;
	typedef itk::SymmetricSecondRankTensor<OutputPixelType, 3> HessianPixelType;
	typedef itk::Vector<OutputPixelType, 3> EigenValuePixelType;
	typedef itk::Vector<EigenValuePixelType, 3> EigenVectorPixelType;

	// itk image type defintions
	typedef itk::Image<InputPixelType, 3> InputImageType;
	typedef itk::Image<OutputPixelType, 3> OutputImageType;
	typedef itk::Image<HessianPixelType, 3> HessianImageType;
	typedef itk::Image<short, 3> MaskImageType;
	typedef itk::ImageRegionIterator<OutputImageType> IteratorType;
	// itk filter defintions
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	typedef itk::CastImageFilter<InputImageType, OutputImageType> CastType;

	//read mask image
	MaskReaderType::Pointer maskReader = MaskReaderType::New();
	maskReader->SetFileName("C:/ImageData/Manual2A_Liver.hdr");
	try
	{
		maskReader->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}
	MaskImageType::Pointer maskImage = maskReader->GetOutput();

	// read input image and cast appropriately
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("C:/ImageData/2A_VenousContrastImage_Smooth.hdr");
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}
	std::cout << "Read Input Image" << std::endl;

	CastType::Pointer cast = CastType::New();
	cast->SetInput(reader->GetOutput());
	try
	{
		cast->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}
	std::cout << "Casted Input Image from (short) to (float)" << std::endl;


	// allocate output image
	OutputImageType::Pointer vesselnessImage = OutputImageType::New();
	vesselnessImage->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
	vesselnessImage->Allocate();
	vesselnessImage->SetSpacing(reader->GetOutput()->GetSpacing());
	vesselnessImage->FillBuffer(0.0);
	typedef itk::ImageRegionIteratorWithIndex<OutputImageType> OutputImageIteratorType;
	OutputImageType::Pointer lambda1;
	OutputImageType::Pointer lambda2;
	OutputImageType::Pointer lambda3;

	float sigmam = reader->GetOutput()->GetSpacing()[0];
	int n = 2;

	// begin multi scale loop
	for(int i = 0; i < n; i++)
	{
		float sigma = sigmam * pow(1.2f, (float)i);
		std::cout << "Computing Vesselness for sigma = " << sigma << std::endl;
		// Generate Hessian of Image
		typedef itk::HessianRecursiveGaussianImageFilter<OutputImageType, HessianImageType> HessianFilterType;
		HessianFilterType::Pointer hessianFilter = HessianFilterType::New();


		hessianFilter->SetSigma(sigma);
		hessianFilter->SetInput(cast->GetOutput());

		try
		{
			hessianFilter->Update();
		}
		catch(itk::ExceptionObject & excp)
		{
			std::cerr << excp << std::endl;
		}

		std::cout << "Generated Hessian Image" << std::endl;

		HessianImageType::Pointer hessianImage = hessianFilter->GetOutput();

		// Perform Eigen Analysis

		typedef itk::SymmetricEigenAnalysis<HessianPixelType, EigenValuePixelType, EigenVectorPixelType> EigenAnalysisType;
		EigenAnalysisType eigen;
		eigen.SetDimension(3);
		eigen.SetOrderEigenValues(true);

		typedef itk::ImageRegionIteratorWithIndex<HessianImageType> HessianIteratorType;
		HessianIteratorType hessianIt(hessianImage, hessianImage->GetLargestPossibleRegion());

		// set up images to hold eigen values
		lambda1 = OutputImageType::New();
		lambda1->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
		lambda1->Allocate();
		lambda1->SetSpacing(reader->GetOutput()->GetSpacing());
		lambda1->FillBuffer(0.0);

		lambda2 = OutputImageType::New();
		lambda2->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
		lambda2->Allocate();
		lambda2->SetSpacing(reader->GetOutput()->GetSpacing());
		lambda2->FillBuffer(0.0);

		lambda3 = OutputImageType::New();
		lambda3->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
		lambda3->Allocate();
		lambda3->SetSpacing(reader->GetOutput()->GetSpacing());
		lambda3->FillBuffer(0.0);



		
		OutputImageIteratorType lambda1It(lambda1, lambda1->GetLargestPossibleRegion());
		OutputImageIteratorType lambda2It(lambda2, lambda2->GetLargestPossibleRegion());
		OutputImageIteratorType lambda3It(lambda3, lambda3->GetLargestPossibleRegion());
		OutputImageIteratorType vesselnessIt(vesselnessImage, vesselnessImage->GetLargestPossibleRegion());

		typedef itk::ImageRegionIteratorWithIndex<MaskImageType> MaskImageIteratorType;
		MaskImageIteratorType maskIt(maskImage, maskImage->GetLargestPossibleRegion());


		// compute vesselness measure
		int counter = 0;
		float vesselnessMeasure;
		for(hessianIt.GoToBegin(), lambda1It.GoToBegin(), lambda2It.GoToBegin(), lambda3It.GoToBegin(), vesselnessIt.GoToBegin(), maskIt.GoToBegin();
			!hessianIt.IsAtEnd() && !lambda1It.IsAtEnd() && !lambda2It.IsAtEnd() && !lambda3It.IsAtEnd() && !vesselnessIt.IsAtEnd() && !maskIt.IsAtEnd();
			++hessianIt, ++ lambda1It, ++lambda2It, ++lambda3It, ++vesselnessIt, ++maskIt)
		{
			if(maskIt.Get() > (short)0)
			{
				HessianImageType::IndexType index;
				HessianImageType::PointType point;
				HessianPixelType hessianMatrix;
				EigenValuePixelType eigenValues;
				EigenVectorPixelType eigenVectors;

				index = hessianIt.GetIndex();
				hessianImage->TransformIndexToPhysicalPoint(index, point);
				hessianMatrix = hessianIt.Get();
				eigen.ComputeEigenValuesAndVectors(hessianMatrix, eigenValues, eigenVectors);

				lambda1It.Set(eigenValues[2]);
				lambda2It.Set(eigenValues[1]);
				lambda3It.Set(eigenValues[0]);
				vesselnessMeasure = /*pow(sigma, 2.0f) * */ vesselness(lambda1It.Get(), lambda2It.Get(), lambda3It.Get());
				if(vesselnessMeasure > vesselnessIt.Get())
					vesselnessIt.Set(vesselnessMeasure);
			}


			if(counter % 500000 == 0)
			{
				std::cout << "Processed " << counter << " voxels" << std::endl;
			}
			counter++;

		}
	
	}

	
	OutputImageIteratorType outputIt(vesselnessImage, vesselnessImage->GetLargestPossibleRegion());
	OutputImageIteratorType l1it(lambda1, lambda1->GetLargestPossibleRegion());
	OutputImageIteratorType l2it(lambda2, lambda2->GetLargestPossibleRegion());
	OutputImageIteratorType l3it(lambda3, lambda3->GetLargestPossibleRegion());
	double l1mean, l2mean, l3mean, l1max, l2max, l3max;
	l1mean = l2mean = l3mean = l1max = l2max = l3max = 0;
	long int counter = 0;
	for(outputIt.GoToBegin(), l1it.GoToBegin(), l2it.GoToBegin(), l3it.GoToBegin(); !outputIt.IsAtEnd(); ++l1it, ++l2it, ++l3it, ++outputIt)
	{
		if(l1it.Get() > l1max) l1max = l1it.Get();
		if(l2it.Get() > l2max) l2max = l2it.Get();
		if(l3it.Get() > l3max) l3max = l3it.Get();
	}
	std::cout << "l1max " << l1max << std::endl;
	std::cout << "l2max " << l2max << std::endl;
	std::cout << "l3max " << l3max << std::endl;


	/*
	// find max value
	int maxvalue = 0;
	for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
	{
		if((int)outputIt.Get() > maxvalue) maxvalue = (int)outputIt.Get();
	}
	std::cout << "Maximum Value = " << maxvalue << std::endl;

	// build histogram
	std::vector<int> histogram;
	histogram.resize(maxvalue + 1, 0);
	for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
	{
		if(outputIt.Get() > 0.0)
		{
			++histogram[(int)outputIt.Get()];
		}
	}
	histogram = buildCumHistogramBottom(histogram);
	int threshold;
	std::cout << histogram[0] << std::endl;
	for(int i = 0 ; i < histogram.size(); i++)
	{
		
		if((float)histogram[i] < (0.08 * (float)histogram[0]))
		{
			std::cout << "Threshold = " << i << std::endl;
			threshold = i;
			break;
		}
	}
	//writeImage<OutputImageType>(vesselnessImage, "test2.img");
	for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
	{
		if((int)outputIt.Get() < threshold) outputIt.Set(0.0);
	}

	

		



	std::cout << "Done Computing Vesselness Values" << std::endl;

	typedef itk::RescaleIntensityImageFilter<OutputImageType, OutputImageType> RescaleFilterType;
	RescaleFilterType::Pointer rescale = RescaleFilterType::New();
	rescale->SetInput(vesselnessImage);
	rescale->SetOutputMinimum(0.50);
	rescale->SetOutputMaximum(1.0);
	try
	{
		rescale->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}
	std::cout << "Done rescaling vesselness measure intensity" << std::endl;
	OutputImageType::Pointer outputImage = rescale->GetOutput();
	IteratorType it(outputImage, outputImage->GetLargestPossibleRegion());
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() == 0.50)
		{
			it.Set(0.0);
		}
	}
	*/
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(lambda1);
	writer->SetFileName("lambda1.img");
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}


	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(lambda2);
	writer2->SetFileName("lambda2.img");
	try
	{
		writer2->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}

	WriterType::Pointer writer3 = WriterType::New();
	writer3->SetInput(lambda3);
	writer3->SetFileName("lambda3.img");
	try
	{
		writer3->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}

	WriterType::Pointer vwriter = WriterType::New();
	vwriter->SetInput(vesselnessImage);
	vwriter->SetFileName("vesselness.hdr");
	try
	{
		vwriter->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
	}


	std::cout << "Done!" << std::endl;
	return 0;
}




float vesselness(float lambda1, float lambda2, float lambda3)
{
	float vesselnessMeasure;
	float vm1, vm2;
	/*
	if((0.0 >= lambda1) && (lambda1 > lambda2) && (lambda2 > lambda3))
	{
		vesselnessMeasure = (abs(lambda3) * pow(lambda2 / lambda3, gamma23)) * pow(1 + (lambda1 / abs(lambda2)), gamma12);
	}

	else if(((abs(lambda2) / alpha) > lambda1) && (lambda1 > 0.0) && (0.0 > lambda2) && (lambda2 > lambda3))
	{
		vesselnessMeasure = (abs(lambda3) * pow(lambda2 / lambda3, gamma23)) * pow(1 - (alpha * (lambda1 / abs(lambda2))), gamma12);
	}*/
	if((lambda1 <= 0) || (lambda2 <= 0) || (lambda3 <= 0))
	{
		vesselnessMeasure = 0.0;
	}
	else
	{
		
		if(lambda1 > 0.10)
		{
			//if((lambda2 / lambda1) > 0.05)
			if(true){
				vm1 = 1.0 - exp(-(lambda1/lambda3));
				vm2 = exp(-((lambda1/lambda2) - 1.0));
				if(vm1 > vm2) vesselnessMeasure = vm1;
				else vesselnessMeasure = vm2;/* exp(-((lambda1/lambda3) - 1.0)); //* exp(-((lambda2/lambda3) - 1.0));//*/ //1.0 - exp(-(lambda1/lambda2));
			}
			else
			{
				vesselnessMeasure = 0.0;
			}
		}
		else vesselnessMeasure = 0.0;
	}

	return vesselnessMeasure;
}