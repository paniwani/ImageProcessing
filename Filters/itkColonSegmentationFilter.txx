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
	InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
	
	// Orient input into LAI orientation (spine is at top of image, lungs at z=0)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<InputImageType,InputImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI);
    orienter->SetInput(input);
    orienter->Update();
	input = orienter->GetOutput();	
	
	InputImageType::RegionType region = input->GetLargestPossibleRegion();
	InputImageType::SpacingType spacing = input->GetSpacing();
	
	typename Superclass::OutputImagePointer output = this->GetOutput(0);
	output->SetRegions( region );
	output->SetOrigin( input->GetOrigin() );
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation( input );
	output->Allocate();
	
	// Detect bone using region growing with seeds from first slice
	// Assumption: only high intensity object at z=0 is bone
	typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType> ConnectedThresholdImageFilterType;
	ConnectedThresholdImageFilterType::Pointer grower = ConnectedThresholdImageFilterType::New();
	grower->SetInput( input );
	grower->SetLower(200);
	grower->SetUpper(1200);
	grower->SetReplaceValue(1);
	
	InputIteratorType input_iter(input,region);

	// Set tagged regions within dilated area as seeds of region growing
	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		InputImageType::IndexType idx = input_iter.GetIndex();
	
		if (idx[2] == 0)
		{
			if (input_iter.Get() > 200)
			{
				grower->AddSeed(idx);
			}
		
		}
	}

	grower->Update();
	ByteImageType::Pointer bone = grower->GetOutput();
	
	WriteITK(bone,"bone.hdr");	
	
	// Find air components
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(spacing);
	air->Allocate();
	
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

	WriteITK(air, "air_lumen.hdr" );

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

	WriteITK(air, "air_median.hdr");

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
	
	WriteITK(air,"air_no_bkg.hdr");	
	
	// Detect lungs using region growing with seeds from first slice
	// Assumption: only low intensity object at z=0 is lung
	typedef itk::ConnectedThresholdImageFilter<ByteImageType, ByteImageType> ConnectedThresholdImageFilterByteType;
	ConnectedThresholdImageFilterByteType::Pointer lungGrower = ConnectedThresholdImageFilterByteType::New();
	lungGrower->SetInput( air );
	lungGrower->SetLower(1);
	lungGrower->SetUpper(1);
	lungGrower->SetReplaceValue(1);

	// Set air regions within dilated area as seeds of region growing
	for (air_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter)
	{
		ByteImageType::IndexType idx = air_iter.GetIndex();
	
		if (idx[2] == 0)
		{
			if (air_iter.Get() == 1)
			{
				lungGrower->AddSeed(idx);
			}
		
		}
	}

	lungGrower->Update();
	ByteImageType::Pointer lung = lungGrower->GetOutput();
	ByteIteratorType lung_iter(lung,region);

	WriteITK(lung,"lung.hdr");
	
	// Remove lung from air mask
	for (air_iter.GoToBegin(), lung_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++lung_iter)
	{
		if (lung_iter.Get() == 1)
			air_iter.Set( 0 );
	}

	WriteITK(air,"air_no_lung.hdr");
	lung.~SmartPointer();
	lungGrower.~SmartPointer();
	
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

	WriteITK(air_dilated, "air_dilated.hdr");

	// Subtract original air area from dilation
	for (air_iter.GoToBegin(), air_dilated_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++air_dilated_iter)
	{
		air_dilated_iter.Set( air_dilated_iter.Get() - air_iter.Get() );
	}

	WriteITK(air_dilated,"air_dilated_only.hdr");

	// Setup region growing filter
	grower = ConnectedThresholdImageFilterType::New();
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
	
	// Remove bone from tagged regions
	ByteIteratorType bone_iter(bone,region);
	
	for (bone_iter.GoToBegin(), tagged_iter.GoToBegin(); !bone_iter.IsAtEnd(); ++bone_iter, ++tagged_iter)
	{
		if (bone_iter.Get() == 1)
		{
			tagged_iter.Set(0);
		}
	}
	

	WriteITK(tagged,"tagged.hdr");

	// Add air lumen to tagged region
	for (tagged_iter.GoToBegin(), air_iter.GoToBegin(); !tagged_iter.IsAtEnd(); ++tagged_iter, ++air_iter)
	{
		if ( air_iter.Get() == 1 )
			tagged_iter.Set(1);
	}

	// Smooth tagged air via median filter
	medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( tagged );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 1 );

	radius.Fill(0);
	radius[0] = 1;
	radius[1] = 1;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	tagged = medianFilter->GetOutput();
	tagged_iter = ByteIteratorType(tagged,region);

	WriteITK(tagged,"tagged_air.hdr");

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

	WriteITK(dilater->GetOutput(),"colon.hdr");
	
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

template <class TInputImage, class TOutputImage >
void
ColonSegmentationFilter<TInputImage,TOutputImage>::
WriteITK(typename ByteImageType::Pointer image, std::string name)
{
	/*
	// Find max of image
	unsigned char max = 0;
	ByteIteratorType image_iter(image,image->GetLargestPossibleRegion());
	for (image_iter.GoToBegin(); !image_iter.IsAtEnd(); ++image_iter)
	{
		if (image_iter.Get() > max)
			max = image_iter.Get();
	}

	// Scale to 255
	if (max > 0)
	{
		for (image_iter.GoToBegin(); !image_iter.IsAtEnd(); ++image_iter)
		{
			image_iter.Set( image_iter.Get() * 255 / max );
		}
	}
	*/

	/*
	typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
	*/
}








} // end namespace itk

#endif
