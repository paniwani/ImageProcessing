const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkConfidenceConnectedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkNoiseImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBlackTopHatImageFilter.h>
#include <itkConvolutionImageFilter.h>
 												
int main(int argc, char * argv[])				
{ 					
	// Load image
	ImageType::Pointer input = ReadITK <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm/mr10_092_13p_i0091.dcm");
	ImageType::RegionType region = input->GetLargestPossibleRegion();
		
	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ImageType> (input,colonRegion);
	//Mask <ImageType,ByteImageType> (input,colon,-1025);

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();	

	// Get gradient magnitude
	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
	gmFilter->SetInput(input);
	gmFilter->Update();
	FloatImageType::Pointer gm = gmFilter->GetOutput();
	WriteITK <FloatImageType> (gm,"gm.nii");

	// Allocate output classified image
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// Show classification with hard threshold at 200
	FloatIteratorType gmIt(gm,region);
	IteratorType inputIt(input,region);
	ByteIteratorType colonIt(colon,region);
	ByteIteratorType outIt(out,region);

	PixelType T = 200;

	for (gmIt.GoToBegin(), inputIt.GoToBegin(), colonIt.GoToBegin(), outIt.GoToBegin(); !gmIt.IsAtEnd(); ++gmIt, ++inputIt, ++colonIt, ++outIt)
	{
		if (colonIt.Get() != 0)
		{
			PixelType I = inputIt.Get();
			float G = gmIt.Get();

			if ( ( I >= T && G < 0.8*I ) || I > 1000 )
			{
				outIt.Set(1);
			} else if ( I <= -700 ) {
				outIt.Set(2);
			} else if ( I < T && I > -300 && G <= 300 ) {
				outIt.Set(3);
			}
		}
	}

	WriteITK <ByteImageType> (out,"smc.nii");

	// Show standard deviation image
	typedef itk::NoiseImageFilter<ImageType,FloatImageType> NoiseFilterType;
	NoiseFilterType::Pointer noiser = NoiseFilterType::New();
	noiser->SetInput(input);
	
	ImageType::SizeType radius;
	radius.Fill(1);
	
	noiser->SetRadius(radius);
	noiser->Update();
	WriteITK <FloatImageType> (noiser->GetOutput(),"noise.nii");

	// Find local minima
	typedef itk::BinaryBallStructuringElement<BytePixelType,Dimension> StructuringElementType;  
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(3);
	se.SetRadius( rad );
	se.CreateStructuringElement();

	typedef itk::BlackTopHatImageFilter<ImageType,ImageType,StructuringElementType> BlackTopHatFilterType;
	BlackTopHatFilterType::Pointer blackTopHat = BlackTopHatFilterType::New();
	blackTopHat->SetInput(input);
	blackTopHat->SetKernel(se);
	blackTopHat->Update();

	ImageType::Pointer bth = blackTopHat->GetOutput();
	WriteITK <ImageType> (bth,"bth.nii");

	// Enhance edges with sharpening Laplacian kernel
	ImageType::Pointer laplacian = ImageType::New();
	
	ImageType::RegionType lr;
	ImageType::IndexType li;
	ImageType::SizeType ls;
	ls[0] = 3;
	ls[1] = 3;

	li[0] = 0;
	li[1] = 0;

	lr.SetIndex(li);
	lr.SetSize(ls);

	laplacian->SetRegions(lr);
	laplacian->Allocate();
	laplacian->FillBuffer(-1);
	
	li[0] = 1; li[1] = 1;

	laplacian->SetPixel(li,9);

	// Convolve
	typedef itk::ConvolutionImageFilter<ImageType> ConvolutionFilterType;
	ConvolutionFilterType::Pointer convolver = ConvolutionFilterType::New();
	convolver->SetInput(input);
	convolver->SetImageKernelInput(laplacian);
	convolver->Update();

	WriteITK <ImageType> (convolver->GetOutput(),"laplacian.nii");




	return 0;

}