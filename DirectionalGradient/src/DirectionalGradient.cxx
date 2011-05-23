#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkColonSegmentationFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>
#include <itkMorphologicalDistanceTransformImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkDirectionalGradientImageFilter.h>
#include <itkDirectionalGradientImageFilter2.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkSubtractImageFilter.h>

int main()
{
	// load image
	ImageType2D::Pointer input = ReadDicom <ImageType2D> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",86);
	WriteITK <ImageType2D> (input,"inputFull.nii");

	ImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// segment colon
	typedef itk::ColonSegmentationFilter< ImageType2D, ByteImageType2D > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( input );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	colonSegmenter->SetRemoveBoneLung(false);
	//colonSegmenter->SetPrintImages(true);
	colonSegmenter->Update();
	ByteImageType2D::Pointer colon = colonSegmenter->GetOutput();
	WriteITK <ByteImageType2D> (colon,"colonFull.nii");

	// crop images
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
	WriteITK <ImageType2D> (input,"input.nii");

	typedef itk::RegionOfInterestImageFilter<ByteImageType2D,ByteImageType2D> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();
	WriteITK <ByteImageType2D> (colon,"colon.nii");

	// reset to new region
	region = input->GetLargestPossibleRegion();

	// mask input with colon
	typedef itk::MaskImageFilter<ImageType2D,ByteImageType2D,ImageType2D> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(input);
	masker->SetInput2(colon);
	masker->SetOutsideValue(-1024);
	masker->Update();
	input = masker->GetOutput();
	WriteITK <ImageType2D> (input,"inputMasked.nii");

	// get gradient
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType2D,ImageType2D> GradientMagnitudeRecursiveGaussianImageFilterType;
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->SetSigma(input->GetSpacing()[0]);
	gradientFilter->Update();
	ImageType2D::Pointer gradient = gradientFilter->GetOutput();
	WriteITK <ImageType2D> (gradient,"gradient.nii");

	// separate tissue and stool using otsu threshold	
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType2D > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1500);
	otsuCalculator->Compute();

	PixelType tissueStoolThreshold = otsuCalculator->GetThreshold();
	std::cout << "Tissue Stool Otsu Threshold: " << tissueStoolThreshold << std::endl;

	// mark air and tissue
	ByteImageType2D::Pointer vmap = ByteImageType2D::New();
	vmap->SetRegions(region);
	vmap->SetSpacing(input->GetSpacing());
	vmap->Allocate();
	vmap->FillBuffer(0);

	ByteIteratorType2D vmapIt(vmap,region);
	IteratorType2D gradientIt(gradient,region);

	inputIt = IteratorType2D(input,region);
	colonIt = ByteIteratorType2D(colon,region);

	for (inputIt.GoToBegin(), vmapIt.GoToBegin(), gradientIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd();
		 ++inputIt, ++vmapIt, ++gradientIt, ++colonIt)
	{
		if (colonIt.Get() == 255)
		{
			unsigned char v = 0;

			PixelType I = inputIt.Get();
			PixelType G = gradientIt.Get();

			if ( ( I >= tissueStoolThreshold && G < 0.8*I ) || I > 1000 )
			{
				v = 1;
			} else if ( I <= -800 ) {
				v = 2;
			} else if ( I < tissueStoolThreshold && I > -300 && G <= 400 ) {
				v = 3;
			}

			vmapIt.Set(v);
		}
	}

	typedef itk::MultiplyByConstantImageFilter<ByteImageType2D, unsigned char, ByteImageType2D> MultiplyByConstantImageFilterType;
	MultiplyByConstantImageFilterType::Pointer multiplier = MultiplyByConstantImageFilterType::New();
	multiplier->SetInput(vmap);
	multiplier->SetConstant( 255/3 );
	multiplier->Update();
	WriteITK <ByteImageType2D> (multiplier->GetOutput(),"vmap.nii");

	// air mask
	ByteImageType2D::Pointer airMask = ByteImageType2D::New();
	airMask->SetRegions(region);
	airMask->SetSpacing(input->GetSpacing());
	airMask->Allocate();
	airMask->FillBuffer(0);
	ByteIteratorType2D airMaskIt(airMask,region);

	for (vmapIt.GoToBegin(), airMaskIt.GoToBegin(); !vmapIt.IsAtEnd(); ++vmapIt, ++airMaskIt)
	{
		if (vmapIt.Get() == 2)
			airMaskIt.Set(255);
	}
	WriteITK <ByteImageType2D> (airMask,"airMask.nii");

	// distance from air
	typedef itk::MorphologicalDistanceTransformImageFilter<ByteImageType2D, ImageType2D> MorphologicalDistanceTransformImageFilterType;
	MorphologicalDistanceTransformImageFilterType::Pointer distanceFilter = MorphologicalDistanceTransformImageFilterType::New();
	distanceFilter->SetInput(airMask);
	distanceFilter->SetOutsideValue(255);
	distanceFilter->Update();

	// set max dist
	ImageType2D::Pointer airDist = distanceFilter->GetOutput();
	IteratorType2D airDistIt(airDist,region);

	for (airDistIt.GoToBegin(); !airDistIt.IsAtEnd(); ++airDistIt)
	{
		if (airDistIt.Get() > 20)
			airDistIt.Set(20);
	}

	WriteITK <ImageType2D> (airDist,"airDistance.nii");

	// dilate air mask
	typedef itk::BinaryBallStructuringElement<PixelType, 2>	StructuringElementType;StructuringElementType se;
		
	unsigned int radius = 3;

	se.SetRadius( radius );
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<ByteImageType2D, ByteImageType2D, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( airMask );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();
	
	typedef itk::SubtractImageFilter<ByteImageType2D> SubtractImageFilterType;
	SubtractImageFilterType::Pointer subtracter = SubtractImageFilterType::New();
	subtracter->SetInput1(dilater->GetOutput());
	subtracter->SetInput2(airMask);
	subtracter->Update();
	ByteImageType2D::Pointer airBorder = subtracter->GetOutput();
	ByteIteratorType2D airBorderIt(airBorder,region);

	for (airBorderIt.GoToBegin(), vmapIt.GoToBegin(); !airBorderIt.IsAtEnd(); ++airBorderIt, ++vmapIt)
	{
		if (vmapIt.Get() != 0)
			airBorderIt.Set(0);
	}

	WriteITK <ByteImageType2D> (airBorder,"airBorder.nii");

	// gradient away from air
	typedef itk::DirectionalGradientImageFilter<ImageType2D,ByteImageType2D,ImageType2D> DirectionalGradientImageFilterType;
	DirectionalGradientImageFilterType::Pointer dgFilter = DirectionalGradientImageFilterType::New();
	dgFilter->SetInput(input);
	dgFilter->SetMaskImage(airMask);
	dgFilter->SetOutsideValue(255);
	
	double sigma[6] = {.1,.25,.5,.75,1,2};

	for (int i=0; i<6; i++)
	{
		dgFilter->SetSigma(sigma[i]);
		dgFilter->Update();

		MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
		masker->SetInput1(dgFilter->GetOutput());
		masker->SetInput2(airBorder);
		masker->Update();
		
		std::stringstream ss;
		ss << "dg_" << sigma[i] << ".nii";
		WriteITK <ImageType2D> (masker->GetOutput(),ss.str());

	}

	system("pause");
	return 0;
}