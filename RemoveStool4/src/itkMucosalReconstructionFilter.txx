#ifndef __itkMucosalReconstructionFilter_txx
#define __itkMucosalReconstructionFilter_txx

#include "itkMucosalReconstructionFilter.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage >
MucosalReconstructionFilter<TInputImage,TOutputImage >
::MucosalReconstructionFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_NumOfLayers = 2;
}


/**
 * GenerateData Performs the reconstruction
 */
template <class TInputImage, class TOutputImage >
void
MucosalReconstructionFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
	// Setup
	typename Superclass::InputImageConstPointer  input = this->GetInput();
	typename Superclass::OutputImagePointer output = this->GetOutput(0);

	InputImageType::RegionType region = input->GetLargestPossibleRegion();

	output->SetRegions( region );
	output->SetOrigin( input->GetOrigin() );
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation( input );
	output->Allocate();
  
	typedef itk::Image<unsigned char, 2> ByteImageType;

	typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > 	InputIteratorType;	
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > 		IteratorType;
	typedef itk::ImageRegionIterator< ByteImageType >  					ByteIteratorType;
  
	// Threshold to detect air
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing( input->GetSpacing() );
	air->Allocate();

	InputIteratorType input_iter(input,region);
	ByteIteratorType air_iter(air,region);

	for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter)
	{
		if (input_iter.Get() < -300) // Chosen to reconstruct entire mucosa (including regions with no stool) for uniformity
		{
			air_iter.Set( 0 );
		} else {
			air_iter.Set( 1 );
		}
	}

	
	// Remove non-air components that are not connected to the largest body of tissue
	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType > BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryFilter->SetInput( air );
	binaryFilter->SetBackgroundValue(0);
	binaryFilter->SetForegroundValue(1);
	binaryFilter->SetAttribute("Size");
	binaryFilter->SetNumberOfObjects(1);
	binaryFilter->Update();
	air = binaryFilter->GetOutput();

	// Invert image to label air components
	air_iter = ByteIteratorType(air,region);
	for (air_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter)
	{
		air_iter.Set( !air_iter.Get() );
	}	

	// Isolate largest air component with largest size on border (background air)
	binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryFilter->SetInput( air );
	binaryFilter->SetBackgroundValue(0);
	binaryFilter->SetForegroundValue(1);
	binaryFilter->SetAttribute("SizeOnBorder");
	binaryFilter->SetNumberOfObjects(1);
	binaryFilter->Update();
	ByteImageType::Pointer bkg = binaryFilter->GetOutput();
	ByteIteratorType bkg_iter(bkg,region);

	for (air_iter.GoToBegin(), bkg_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++bkg_iter)
	{
		if (bkg_iter.Get() == 1)
			air_iter.Set( 0 );
	}

	IteratorType output_iter(output,region);

	for (input_iter.GoToBegin(), output_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter, ++air_iter)
	{
		output_iter.Set( input_iter.Get() );

		if (air_iter.Get() == 1)
		{
			output_iter.Set( -1025 );
		}
	}

	/****
		Determine initial starting tissue value To by dilating air boundary, finding edges, and averaging intensity
	****/

	// Dilate by 1
	typedef itk::BinaryBallStructuringElement<unsigned char, 2> StructuringElementType;

	ByteImageType::SizeType radius;
	radius[0] = 1;
	radius[1] = 1;

	StructuringElementType structuringElement;
	structuringElement.SetRadius( radius );
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput( air );
	dilateFilter->SetForegroundValue(1);
	dilateFilter->SetBackgroundValue(0);
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->Update();
	ByteImageType::Pointer air_dilate = dilateFilter->GetOutput();

	// Find edges of air
	typedef itk::BinaryContourImageFilter< ByteImageType, ByteImageType > BinaryContourImageFilterType;
	BinaryContourImageFilterType::Pointer binaryEdgeFilter = BinaryContourImageFilterType::New();
	binaryEdgeFilter->SetInput( air_dilate );
	binaryEdgeFilter->SetForegroundValue(1);
	binaryEdgeFilter->SetBackgroundValue(0);
	binaryEdgeFilter->Update();

	ByteImageType::Pointer air_edge = binaryEdgeFilter->GetOutput();

	// Calculate average starting tissue value (T0)
	ByteIteratorType air_edge_iter(air_edge,region);

	float T0 = 0;
	unsigned int count = 0;

	for (input_iter.GoToBegin(), air_edge_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_edge_iter)
	{
		if ( air_edge_iter.Get() == 1 )
		{
			if ( input_iter.Get() >= -250 && input_iter.Get() <= 150 )
			{
				T0 += input_iter.Get();
				count++;
			}
		}
	}

	T0 /= count;

	/***
		Erode air mask for each layer, find eges, and apply scaffolded tissue value
	***/

	// Make mask of set of scaffolding pixels
	ByteImageType::Pointer scaffold_mask = ByteImageType::New();
	scaffold_mask->SetRegions( region );
	scaffold_mask->SetSpacing( input->GetSpacing() );
	scaffold_mask->Allocate();
	scaffold_mask->FillBuffer(0);
	ByteIteratorType scaffold_mask_iter(scaffold_mask,region);
	
	typedef itk::BinaryErodeImageFilter< ByteImageType, ByteImageType, StructuringElementType > BinaryErodeImageFilterType;

	float T = T0;

	for (int i=0; i < m_NumOfLayers; i++)
	{
		T += (-1000-T0)/(m_NumOfLayers+1);

		if (i>0)
		{
			BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
			erodeFilter->SetInput( air );
			erodeFilter->SetForegroundValue(1);
			erodeFilter->SetBackgroundValue(0);
			erodeFilter->SetKernel( structuringElement );
			erodeFilter->Update();
			air = erodeFilter->GetOutput();
		}

		binaryEdgeFilter = BinaryContourImageFilterType::New();
		binaryEdgeFilter->SetInput( air );
		binaryEdgeFilter->SetForegroundValue(1);
		binaryEdgeFilter->SetBackgroundValue(0);
		binaryEdgeFilter->SetFullyConnected(true);
		binaryEdgeFilter->Update();
		air_edge = binaryEdgeFilter->GetOutput();

		air_edge_iter = ByteIteratorType(air_edge,region);

		for (output_iter.GoToBegin(), air_edge_iter.GoToBegin(), scaffold_mask_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++air_edge_iter, ++scaffold_mask_iter)
		{
			if (air_edge_iter.Get() == 1)
			{
				output_iter.Set( T );

				scaffold_mask_iter.Set( 1 );
			}
		}
	}

	// Dilate scaffold mask
	dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput( scaffold_mask );
	dilateFilter->SetForegroundValue(1);
	dilateFilter->SetBackgroundValue(0);
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->Update();
	scaffold_mask = dilateFilter->GetOutput();
	scaffold_mask_iter = ByteIteratorType(scaffold_mask,region);

	// Smooth output within scaffold mask
	typedef itk::GaussianBlurImageFunction< OutputImageType > GaussianBlurImageFunctionType;
	GaussianBlurImageFunctionType::Pointer gaussianFunction = GaussianBlurImageFunctionType::New();
	gaussianFunction->SetInputImage( output );

	GaussianBlurImageFunctionType::ErrorArrayType setError;
	setError.Fill( 0.01 );
	gaussianFunction->SetMaximumError( setError );
	gaussianFunction->SetSigma( 0.7 );

	for (output_iter.GoToBegin(), scaffold_mask_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++scaffold_mask_iter)
	{
		if (scaffold_mask_iter.Get() == 1)
		{
			output_iter.Set( gaussianFunction->EvaluateAtIndex( output_iter.GetIndex() ) );
		}
	}
	
}

template <class TInputImage, class TOutputImage >
void
MucosalReconstructionFilter<TInputImage,TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "NumOfLayers: " << m_NumOfLayers << std::endl;
}

} // end namespace itk

#endif
