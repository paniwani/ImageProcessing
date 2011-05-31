#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <otbScalarImageToTexturesFilter.h>
#include <otbScalarImageToTexturesMaskFilter.h>
#include <itkColonSegmentationFilter.h>
#include <itkMaskImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkNeighborhoodIterator.h>
#include <itkZeroCrossingBasedEdgeDetectionImageFilter.h>
#include <itkFixedArray.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkCastImageFilter.h>
#include <time.h>
#include <itkDivideByConstantImageFilter.h>
#include <itkAddImageFilter.h>

typedef float PixelType;

const unsigned int Dimension = 2;

typedef itk::Image<PixelType,Dimension> ImageType;
typedef itk::Image<unsigned char,Dimension> ByteImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType> ByteIteratorType;

ImageType::Pointer StandardDeviation(ImageType::Pointer &input, unsigned int radius);
ImageType::Pointer Range(ImageType::Pointer &input, unsigned int radius);

int main(int argc, char * argv[])
{
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory radius";
		system("pause");
		return EXIT_FAILURE;
	}

	// start the clock
	clock_t init;
	init = clock();	

	// load image
	ImageType::Pointer input = ReadDicom <ImageType> (argv[1],85);
	WriteITK <ImageType> (input, "inputOriginal.nii");

	// set radius
	const unsigned int radius = atoi(argv[2]);

	ImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();

	ImageType::IndexType cropIndex;
	cropIndex[0] = 350;
	cropIndex[1] = 230;
	//cropIndex[2] = 0;

	ImageType::SizeType cropSize;
	cropSize[0] = 50;
	cropSize[1] = 50;
	//cropSize[2] = inputSize[2];

	ImageType::RegionType cropRegion;
	cropRegion.SetIndex(cropIndex);
	cropRegion.SetSize(cropSize);

	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(cropRegion);
	cropper->Update();
	input = cropper->GetOutput();
	WriteITK <ImageType> (input,"inputCropped.nii");

	// text output
	std::stringstream ss;

	//// segment colon
	//typedef itk::ColonSegmentationFilter<ImageType, ByteImageType> ColonSegmentationFilterType;
	//ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	//colonSegmenter->SetInput(input);
	//colonSegmenter->SetOutputForegroundValue(255);
	//colonSegmenter->SetOutputBackgroundValue(0);
	//colonSegmenter->SetRemoveBoneLung(false);
	//colonSegmenter->Update();
	//ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	//// crop to colon region
	//Crop(input,colon);

	//// mask input with colon
	//typedef itk::MaskImageFilter<ImageType,ByteImageType,ImageType> MaskImageFilterType;
	//MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	//masker->SetInput1(input);
	//masker->SetInput2(colon);
	//masker->SetOutsideValue(-1024);
	//masker->Update();

	//input = masker->GetOutput();

	//WriteITK <ImageType> (input,"input.nii");
	//WriteITK <ByteImageType> (colon,"colon.nii");
	
	/*
	// median filter on input
	typedef itk::MedianImageFilter<ImageType,ImageType> MedianImageFilterType;
	MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
	medianFilter->SetInput(input);

	ImageType::SizeType medianRadius;
	medianRadius.Fill(radius);

	medianFilter->SetRadius(medianRadius);
	medianFilter->Update();
	input = medianFilter->GetOutput();

	ss.str("");
	ss << "inputMedian" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType> (input,ss.str());
	*/

	// get sobel
	typedef itk::SobelEdgeDetectionImageFilter<ImageType,ImageType> SobelEdgeDetectionImageFilterType;
	SobelEdgeDetectionImageFilterType::Pointer sobelFilter = SobelEdgeDetectionImageFilterType::New();
	sobelFilter->SetInput(input);
	sobelFilter->Update();
	WriteITK <ImageType> (sobelFilter->GetOutput(),"sobel.nii");

	// get gradient
	typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	WriteITK <ImageType> (gmFilter->GetOutput(),"gm.nii");

	// get gradient smoothed
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gsmFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gsmFilter->SetInput(input);

	float gsigma[6] = {0.1,0.25,0.5,1,2,5};

	for (int i=0; i<6; i++)
	{
		gsmFilter->SetSigma(gsigma[i]);
		gsmFilter->Update();

		ss.str("");
		ss << "gsm_" << gsigma[i] << ".nii";
		WriteITK <ImageType> (gsmFilter->GetOutput(),ss.str());
	}

	// get laplacian zero crossings
	typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter<ImageType,ImageType> ZeroCrossingBasedEdgeDetectionImageFilterType;
	ZeroCrossingBasedEdgeDetectionImageFilterType::Pointer zeroCrossingFilter = ZeroCrossingBasedEdgeDetectionImageFilterType::New();
	zeroCrossingFilter->SetInput(input);
	zeroCrossingFilter->SetBackgroundValue(0);
	zeroCrossingFilter->SetForegroundValue(255);

	// set variance
	typedef itk::FixedArray<double,2> ArrayType;
	ArrayType variance;

	for (int i=5; i<=20; i+=5)
	{
		variance[0] = i;
		variance[1] = i;

		zeroCrossingFilter->SetVariance(variance);
		zeroCrossingFilter->Update();

		ss.str("");
		ss << "zeroCrossing_" << i << ".nii";
		WriteITK <ImageType> (zeroCrossingFilter->GetOutput(),ss.str());
	}

	// get standard deviation image
	ImageType::Pointer std = StandardDeviation(input,radius);
	
	ss.str("");
	ss << "std" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType> (std,ss.str());

	// get range image
	ImageType::Pointer range = Range(input,radius);
	
	ss.str("");
	ss << "range" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType> (range,ss.str());	

	// get textures

	//// compress image to 0-255 unsigned char
	//typedef itk::RescaleIntensityImageFilter<ImageType,ByteImageType> RescaleIntensityImageFilterType;
	//RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
	//rescaler->SetInput(input);
	//rescaler->SetOutputMaximum(255);
	//rescaler->SetOutputMinimum(0);
	//rescaler->Update();
	//WriteITK <ByteImageType> (rescaler->GetOutput(),"inputCompressed.nii");

	// rescale input to 0-255
	typedef itk::RescaleIntensityImageFilter<ImageType> RescaleIntensityImageFilterType;
	RescaleIntensityImageFilterType::Pointer rescaler = RescaleIntensityImageFilterType::New();
	rescaler->SetInput(input);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	WriteITK <ImageType> (rescaler->GetOutput(),"inputRescaled.nii");

	// get image min and max
	/*typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage(input);
	minMaxCalc->Compute();*/

	// set all offets to average across
	// half of all possible directions from center pixel
	std::vector<ImageType::OffsetType> offsetVector;

	typedef itk::Neighborhood<PixelType,ImageType::ImageDimension> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(1);

	for (unsigned int i=0; i<(hood.Size()-1)/2; i++)
	{
		ImageType::OffsetType offset = hood.GetOffset(i);
		offsetVector.push_back(offset);
	}

	// allocate all 8 outputs
	ImageType::Pointer outputArray[8];
	for (int i=0; i<8; i++)
	{
		outputArray[i] = ImageType::New();
		outputArray[i]->SetRegions(input->GetLargestPossibleRegion());
		outputArray[i]->SetSpacing(input->GetSpacing());
		outputArray[i]->Allocate();
		outputArray[i]->FillBuffer(0);
	}

	// get texture for all
	unsigned int numOffsets = offsetVector.size();
	
	std::cout << "Number of offsets: " << numOffsets << std::endl;

	for (int i=0; i<numOffsets; i++)
	{
		typedef otb::ScalarImageToTexturesFilter<ImageType,ImageType> ScalarImageToTexturesFilterType;
		ScalarImageToTexturesFilterType::Pointer textureFilter = ScalarImageToTexturesFilterType::New();
		textureFilter->SetInput(rescaler->GetOutput());
		//textureFilter->SetMaskImage(colon);
		textureFilter->SetInputImageMinimum(0);
		textureFilter->SetInputImageMaximum(255);
		textureFilter->SetNumberOfThreads(1);
		
		std::cout << "Number of threads: " << textureFilter->GetNumberOfThreads() << std::endl;

		ImageType::SizeType rad;
		rad.Fill(radius);

		textureFilter->SetRadius(rad);

		/*
		ImageType::OffsetType offset;
		offset[0] = 1;
		offset[1] = 1;
		*/

		textureFilter->SetOffset(offsetVector[i]);
		
		textureFilter->Update();
		
		// get all outputs
		for (unsigned int j=0; j<8; j++)
		{
			// scale output to 0-255
			rescaler = RescaleIntensityImageFilterType::New();
			rescaler->SetInput(textureFilter->GetOutput(j));
			rescaler->SetOutputMaximum(255);
			rescaler->SetOutputMinimum(0);
			rescaler->Update();

			// cast to byte
			typedef itk::CastImageFilter<ImageType,ByteImageType> CastImageFilterType;
			CastImageFilterType::Pointer caster = CastImageFilterType::New();
			caster->SetInput(rescaler->GetOutput());
			caster->Update();

			ss.str("");
			ss << "texture" << j << "_offset_" << offsetVector[i][0] << "_" << offsetVector[i][1] << "_" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
			WriteITK <ByteImageType> (caster->GetOutput(),ss.str());

			// sum textures from each offet
			typedef itk::AddImageFilter<ImageType,ByteImageType,ImageType> AddImageFilterType;
			AddImageFilterType::Pointer adder = AddImageFilterType::New();
			adder->SetInput1(outputArray[j]);
			adder->SetInput2(caster->GetOutput());
			adder->Update();
			outputArray[j] = adder->GetOutput();			
		}

		std::cout << "Completed offset number " << i+1 << std::endl;
		std::cout << "Time elapsed: " << ((double) clock()-init) / (double) CLOCKS_PER_SEC << std::endl;
		init = clock();
	}

	// average all 8 textures across all offsets
	// output as ImageType...
	for (int i=0; i<8; i++)
	{
		typedef itk::DivideByConstantImageFilter<ImageType,PixelType,ImageType> DivideByConstantImageFilterType;
		DivideByConstantImageFilterType::Pointer divider = DivideByConstantImageFilterType::New();
		divider->SetInput(outputArray[i]);
		divider->SetConstant( numOffsets );
		divider->Update();

		typedef itk::CastImageFilter<ImageType,ByteImageType> CastImageFilterType;
		CastImageFilterType::Pointer caster = CastImageFilterType::New();
		caster->SetInput(divider->GetOutput());
		caster->Update();

		ss.str("");
		ss << "texture" << i << "_" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
		WriteITK <ByteImageType> (caster->GetOutput(),ss.str());
	}

	std::cout << "Time elapsed: " << ((double) clock()-init) / (double) CLOCKS_PER_SEC << std::endl;
	system("pause");
	return 0;
}

//void Crop(ImageType::Pointer &input, ByteImageType::Pointer &colon)
//{
//	ImageType::RegionType region = input->GetLargestPossibleRegion();
//	ImageType::SizeType size = region.GetSize();
//
//	long minX=size[0],minY=size[1],maxX=0,maxY=0;
//	short paddingXY = 5;
//	
//	IteratorType2D inputIt(input,region);
//	ByteIteratorType2D colonIt(colon,region);
//
//	for (inputIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt)
//	{
//		if ( colonIt.Get() != 0 )
//		{
//			ImageType::IndexType idx = inputIt.GetIndex();
//
//			if (idx[0] < minX)
//				minX = idx[0];
//			if (idx[0] > maxX)
//				maxX = idx[0];
//			if (idx[1] < minY)
//				minY = idx[1];
//			if (idx[1] > maxY)
//				maxY = idx[1];
//		}
//	}
//
//	ImageType::IndexType edx;
//	
//	edx[0] = (minX-paddingXY) > 0 ? minX-paddingXY : 0;
//	edx[1] = (minY-paddingXY) > 0 ? minY-paddingXY : 0;
//
//	ImageType::SizeType esize;
//	esize[0] = maxX-minX+2*paddingXY+1 < size[0] ? maxX-minX+2*paddingXY+1 : size[0];
//	esize[1] = maxY-minY+2*paddingXY+1 < size[1] ? maxY-minY+2*paddingXY+1 : size[1];
//
//	ImageType::RegionType extractRegion;
//	extractRegion.SetIndex( edx );
//	extractRegion.SetSize( esize );
//
//	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> RegionOfInterestImageFilterType;
//	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
//	cropper->SetInput( input );
//	cropper->SetRegionOfInterest( extractRegion );
//	cropper->Update();
//	input = cropper->GetOutput();
//
//	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> RegionOfInterestImageFilterByteType;
//	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
//	cropperByte->SetInput( colon );
//	cropperByte->SetRegionOfInterest( extractRegion );
//	cropperByte->Update();
//	colon = cropperByte->GetOutput();
//}

// NOTE IMAGE TYPE MUST BE FLOAT
ImageType::Pointer StandardDeviation(ImageType::Pointer &input, unsigned int radius)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	ImageType::Pointer out = ImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	IteratorType outIt(out,region);
	
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt)
	{
		float mean = 0;
		float std = 0;

		// get mean
		for (int i=0; i<nIt.Size(); i++)
		{
			mean += (float) nIt.GetPixel(i);
		}

		mean /= nIt.Size();

		//std::cout << mean << std::endl;

		// get standard deviation
		for (int i=0; i<nIt.Size(); i++)
		{
			std += ( (float) nIt.GetPixel(i) - mean ) * ( (float) nIt.GetPixel(i) - mean );
		}

		std /= (nIt.Size()-1);
		std = sqrt(std);

		outIt.Set(std);
	}

	return out;
}

ImageType::Pointer Range(ImageType::Pointer &input, unsigned int radius)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	ImageType::Pointer out = ImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	IteratorType outIt(out,region);
	
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt)
	{
		
		PixelType max = itk::NumericTraits<PixelType>::NonpositiveMin();
		PixelType min = itk::NumericTraits<PixelType>::max();
		
		for (int i=0; i<nIt.Size(); i++)
		{
			PixelType val = nIt.GetPixel(i);

			if (val > max)
				max = val;

			if (val < min)
				min = val;
		}

		outIt.Set(max-min);
	}

	return out;
}