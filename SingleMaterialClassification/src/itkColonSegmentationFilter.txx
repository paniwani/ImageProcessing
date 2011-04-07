#ifndef ColonSegmentationFilter_txx
#define ColonSegmentationFilter_txx

#include "itkColonSegmentationFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage >
ColonSegmentationFilter<TInputImage,TOutputImage >
::ColonSegmentationFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_TaggedValue= 200;
  m_ForegroundValue=1;
  m_BackgroundValue=0;
}


/**
 * GenerateData Performs the reconstruction
 */
template <class TInputImage, class TOutputImage >
void
ColonSegmentationFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
	// Setup
	typename Superclass::InputImageConstPointer  input = this->GetInput();
	typename Superclass::OutputImagePointer output = this->GetOutput(0);

	InputImageType::RegionType region = input->GetLargestPossibleRegion();
	InputImageType::SpacingType spacing = input->GetSpacing();

	output->SetRegions( region );
	output->SetOrigin( input->GetOrigin() );
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation( input );
	output->Allocate();
  
	typedef itk::Image<unsigned char, InputImageDimension> ByteImageType;

	typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > 	InputIteratorType;	
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > 		IteratorType;
	typedef itk::ImageRegionIterator< ByteImageType >  					ByteIteratorType;
	
	
	// Find air components
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(spacing);
	air->Allocate();
	
	InputIteratorType input_iter(input,region);
	ByteIteratorType air_iter(air,region);

	for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter)
	{
		if (input_iter.Get() < -600)
		{
			air_iter.Set(1);
		} else {
			air_iter.Set(0);
		}
	}

	//WriteITK(air, "air_lumen.hdr" );

	// Apply median filter to remove noise
	typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType > BinaryMedianFilterType;
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( air );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 1 );

	ByteImageType::SizeType radius;
	radius.Fill(0);
	radius[0] = 4;
	radius[1] = 4;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	air = medianFilter->GetOutput();
	air_iter = ByteIteratorType(air,region);

	//WriteITK(air, "air_median.hdr");

	// Remove background air component
	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType > BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
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
	
	binaryFilter.~SmartPointer();
	bkg.~SmartPointer();
	
	// Dilate air
	typedef itk::BinaryBallStructuringElement<unsigned char, InputImageDimension> StructuringElementType;

	radius.Fill(0);
	radius[0] = 1;
	radius[1] = 1;

	StructuringElementType ball;
	ball.SetRadius( radius );
	ball.CreateStructuringElement();
	
	typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType> DilateFilterType;

	DilateFilterType::Pointer dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(1);
	dilater->SetInput( air );
	dilater->SetKernel( ball );
	dilater->Update();
	ByteImageType::Pointer air_dilated = dilater->GetOutput();
	ByteIteratorType air_dilated_iter(air_dilated,region);

	//WriteITK(air_dilated, "air_dilated.hdr");

	// Subtract original air area from dilation
	for (air_iter.GoToBegin(), air_dilated_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++air_dilated_iter)
	{
		air_dilated_iter.Set( air_dilated_iter.Get() - air_iter.Get() );
	}

	//WriteITK(air_dilated,"air_dilated_only.hdr");

	// Setup region growing filter
	typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType> ConnectedThresholdImageFilterType;
	ConnectedThresholdImageFilterType::Pointer grower = ConnectedThresholdImageFilterType::New();
	grower->SetInput( input );
	grower->SetLower(m_TaggedValue);
	grower->SetUpper(1200);
	grower->SetReplaceValue(1);

	// Set tagged regions within dilated area as seeds of region growing
	for (input_iter.GoToBegin(), air_dilated_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_dilated_iter)
	{
		if (air_dilated_iter.Get() == 1 && input_iter.Get() >= m_TaggedValue)
		{
			grower->AddSeed( input_iter.GetIndex() );
		}
	}

	grower->Update();
	ByteImageType::Pointer tagged = grower->GetOutput();
	ByteIteratorType tagged_iter(tagged,region);

	//WriteITK(tagged,"tagged.hdr");

	// Add air lumen to tagged region
	for (tagged_iter.GoToBegin(), air_iter.GoToBegin(); !tagged_iter.IsAtEnd(); ++tagged_iter, ++air_iter)
	{
		if ( air_iter.Get() == 1 )
			tagged_iter.Set(1);
	}

	//WriteITK(tagged,"tagged_air.hdr");

	// Dilate
	dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(1);
	dilater->SetInput( tagged );
	
	radius.Fill(0);
	radius[0] = 8;
	radius[1] = 8;

	ball.SetRadius( radius );
	ball.CreateStructuringElement();

	dilater->SetKernel( ball );
	dilater->Update();

	//WriteITK(dilater->GetOutput(),"colon.hdr");
	
	ByteImageType::Pointer colon = dilater->GetOutput();
	
	// Cast to output image type
	IteratorType output_iter(output,region);
	ByteIteratorType colon_iter(colon,region);
	
	for (output_iter.GoToBegin(), colon_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++colon_iter)
	{
		if ( colon_iter.Get() == 1 )
		{
			output_iter.Set( m_ForegroundValue );
		} else {
			output_iter.Set( m_BackgroundValue );
		}
	}	
}

template <class TInputImage, class TOutputImage >
void
ColonSegmentationFilter<TInputImage,TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "TaggedValue: " << m_TaggedValue << std::endl;
  os << indent << "ForegroundValue: " << m_ForegroundValue << std::endl;
  os << indent << "BackgroundValue: " << m_BackgroundValue << std::endl;
}

} // end namespace itk

#endif
