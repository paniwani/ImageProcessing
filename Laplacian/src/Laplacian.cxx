#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryContourImageFilter.h>

int main()
{
	// load image
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85,90);
	WriteITK <ImageType> (input,"input.nii");

	// get laplacian
	typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType,ImageType> LaplacianRecursiveGaussianImageFilterType;
	LaplacianRecursiveGaussianImageFilterType::Pointer laplacianFilter = LaplacianRecursiveGaussianImageFilterType::New();
	laplacianFilter->SetInput(input);
	laplacianFilter->SetSigma( 2*input->GetSpacing()[0] );
	laplacianFilter->Update();
	ImageType::Pointer laplacian = laplacianFilter->GetOutput();
	WriteITK <ImageType> (laplacian,"laplacian.nii");

	// binarize laplacian
	typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> BinaryThresholdImageFilterType;
	BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
	thresholdFilter->SetInput(laplacian);
	thresholdFilter->SetInsideValue(255);
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetLowerThreshold(0);
	thresholdFilter->Update();
	WriteITK <ByteImageType> (thresholdFilter->GetOutput(),"laplacianBinary.nii");

	// get zero crossings
	typedef itk::BinaryContourImageFilter<ByteImageType,ByteImageType> BinaryContourImageFilterType;
	BinaryContourImageFilterType::Pointer edgeFinder = BinaryContourImageFilterType::New();
	edgeFinder->SetInput(thresholdFilter->GetOutput());
	edgeFinder->SetBackgroundValue(0);
	edgeFinder->SetForegroundValue(255);
	edgeFinder->Update();
	WriteITK <ByteImageType> (edgeFinder->GetOutput(),"laplacianZeroCrossings.nii");

	system("pause");
	return 0;
}