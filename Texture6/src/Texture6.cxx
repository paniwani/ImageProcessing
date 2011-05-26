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

void Crop(ImageType2D::Pointer &input, ByteImageType2D::Pointer &colon);
ImageType2D::Pointer StandardDeviation(ImageType2D::Pointer &input, ByteImageType2D::Pointer &mask, unsigned int radius);
ImageType2D::Pointer Range(ImageType2D::Pointer &input, ByteImageType2D::Pointer &mask, unsigned int radius);

int main()
{
	// load image
	ImageType2D::Pointer input = ReadDicom <ImageType2D> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",86);
	WriteITK <ImageType2D> (input, "inputOriginal.nii");

	// set radius
	const unsigned int radius = 1;

	// text output
	std::stringstream ss;

	// segment colon
	typedef itk::ColonSegmentationFilter<ImageType2D, ByteImageType2D> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetOutputForegroundValue(255);
	colonSegmenter->SetOutputBackgroundValue(0);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType2D::Pointer colon = colonSegmenter->GetOutput();

	// crop to colon region
	Crop(input,colon);

	// mask input with colon
	typedef itk::MaskImageFilter<ImageType2D,ByteImageType2D,ImageType2D> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(input);
	masker->SetInput2(colon);
	masker->SetOutsideValue(-1024);
	masker->Update();

	input = masker->GetOutput();

	WriteITK <ImageType2D> (input,"input.nii");
	WriteITK <ByteImageType2D> (colon,"colon.nii");
	
	/*
	// median filter on input
	typedef itk::MedianImageFilter<ImageType2D,ImageType2D> MedianImageFilterType;
	MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
	medianFilter->SetInput(input);

	ImageType2D::SizeType medianRadius;
	medianRadius.Fill(radius);

	medianFilter->SetRadius(medianRadius);
	medianFilter->Update();
	input = medianFilter->GetOutput();

	ss.str("");
	ss << "inputMedian" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType2D> (input,ss.str());
	*/

	// get sobel
	typedef itk::SobelEdgeDetectionImageFilter<ImageType2D,ImageType2D> SobelEdgeDetectionImageFilterType;
	SobelEdgeDetectionImageFilterType::Pointer sobelFilter = SobelEdgeDetectionImageFilterType::New();
	sobelFilter->SetInput(input);
	sobelFilter->Update();
	WriteITK <ImageType2D> (sobelFilter->GetOutput(),"sobel.nii");

	// get gradient
	typedef itk::GradientMagnitudeImageFilter<ImageType2D,ImageType2D> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	WriteITK <ImageType2D> (gmFilter->GetOutput(),"gm.nii");

	// get gradient smoothed
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType2D> GradientMagnitudeRecursiveGaussianImageFilterType;
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gsmFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gsmFilter->SetInput(input);

	float gsigma[6] = {0.1,0.25,0.5,1,2,5};

	for (int i=0; i<6; i++)
	{
		gsmFilter->SetSigma(gsigma[i]);
		gsmFilter->Update();

		ss.str("");
		ss << "gsm_" << gsigma[i] << ".nii";
		WriteITK <ImageType2D> (gsmFilter->GetOutput(),ss.str());
	}

	// get laplacian zero crossings
	typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter<ImageType2D,ImageType2D> ZeroCrossingBasedEdgeDetectionImageFilterType;
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
		WriteITK <ImageType2D> (zeroCrossingFilter->GetOutput(),ss.str());
	}

	// get standard deviation image
	ImageType2D::Pointer std = StandardDeviation(input,colon,radius);
	
	ss.str("");
	ss << "std" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType2D> (std,ss.str());

	// get range image
	ImageType2D::Pointer range = Range(input,colon,radius);
	
	ss.str("");
	ss << "range" << 2*radius+1 << "x" << 2*radius+1 << ".nii";
	WriteITK <ImageType2D> (range,ss.str());

	// get image min and max
	typedef itk::MinimumMaximumImageCalculator<ImageType2D> MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage(input);
	minMaxCalc->Compute();

	// get textures
	typedef otb::ScalarImageToTexturesMaskFilter<ImageType2D,ImageType2D,ByteImageType2D> ScalarImageToTexturesFilterType;
	ScalarImageToTexturesFilterType::Pointer textureFilter = ScalarImageToTexturesFilterType::New();
	textureFilter->SetInput(input);
	textureFilter->SetMaskImage(colon);
	textureFilter->SetInputImageMinimum(minMaxCalc->GetMinimum());
	textureFilter->SetInputImageMaximum(minMaxCalc->GetMaximum());

	ImageType2D::SizeType rad;
	rad.Fill(radius);

	textureFilter->SetRadius(rad);

	ImageType2D::OffsetType offset;
	offset[0] = 1;
	offset[1] = 1; 

	textureFilter->SetOffset(offset);
	
	textureFilter->Update();
	
	// get all outputs
	for (unsigned int i=0; i<8; i++)
	{
		std::stringstream ss;
		ss << "texture" << i << ".nii";
		WriteITK <ImageType2D> (textureFilter->GetOutput(i),ss.str());
	}

	system("pause");
	return 0;
}

void Crop(ImageType2D::Pointer &input, ByteImageType2D::Pointer &colon)
{
	ImageType2D::RegionType region = input->GetLargestPossibleRegion();
	ImageType2D::SizeType size = region.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0;
	short paddingXY = 5;
	
	IteratorType2D inputIt(input,region);
	ByteIteratorType2D colonIt(colon,region);

	for (inputIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt)
	{
		if ( colonIt.Get() != 0 )
		{
			ImageType2D::IndexType idx = inputIt.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];
		}
	}

	ImageType2D::IndexType edx;
	
	edx[0] = (minX-paddingXY) > 0 ? minX-paddingXY : 0;
	edx[1] = (minY-paddingXY) > 0 ? minY-paddingXY : 0;

	ImageType2D::SizeType esize;
	esize[0] = maxX-minX+2*paddingXY+1 < size[0] ? maxX-minX+2*paddingXY+1 : size[0];
	esize[1] = maxY-minY+2*paddingXY+1 < size[1] ? maxY-minY+2*paddingXY+1 : size[1];

	ImageType2D::RegionType extractRegion;
	extractRegion.SetIndex( edx );
	extractRegion.SetSize( esize );

	typedef itk::RegionOfInterestImageFilter<ImageType2D,ImageType2D> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( input );
	cropper->SetRegionOfInterest( extractRegion );
	cropper->Update();
	input = cropper->GetOutput();

	typedef itk::RegionOfInterestImageFilter<ByteImageType2D,ByteImageType2D> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();
}

// NOTE IMAGE TYPE MUST BE FLOAT
ImageType2D::Pointer StandardDeviation(ImageType2D::Pointer &input, ByteImageType2D::Pointer &mask, unsigned int radius)
{
	// get region
	ImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType2D::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	ImageType2D::Pointer out = ImageType2D::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	IteratorType2D outIt(out,region);
	ByteIteratorType2D maskIt(mask,region);
	
	typedef itk::NeighborhoodIterator<ImageType2D> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(), maskIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt, ++maskIt)
	{
		// inside mask
		if (maskIt.Get() != 0)
		{
			float mean = 0;
			float std = 0;

			// get mean
			for (int i=0; i<nIt.Size(); i++)
			{
				mean += nIt.GetPixel(i);
			}

			mean /= nIt.Size();

			//std::cout << mean << std::endl;

			// get standard deviation
			for (int i=0; i<nIt.Size(); i++)
			{
				std += ( nIt.GetPixel(i) - mean ) * ( nIt.GetPixel(i) - mean );
			}

			std /= (nIt.Size()-1);
			std = sqrt(std);

			

			outIt.Set(std);
		}
	}

	return out;
}

ImageType2D::Pointer Range(ImageType2D::Pointer &input, ByteImageType2D::Pointer &mask, unsigned int radius)
{
	// get region
	ImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType2D::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	ImageType2D::Pointer out = ImageType2D::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	IteratorType2D outIt(out,region);
	ByteIteratorType2D maskIt(mask,region);
	
	typedef itk::NeighborhoodIterator<ImageType2D> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(), maskIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt, ++maskIt)
	{
		// inside mask
		if (maskIt.Get() != 0)
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
	}

	return out;
}