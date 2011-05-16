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
  m_TaggedValue=200;
  m_OutputForegroundValue=255;
  m_OutputBackgroundValue=0;
  m_PrintImages=false;
  m_RemoveBoneLung=true;
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
	
	InputImageType::RegionType region = input->GetLargestPossibleRegion();
	InputImageType::SpacingType spacing = input->GetSpacing();
	InputIteratorType input_iter(input,region);
	
	typename Superclass::OutputImagePointer output = this->GetOutput(0);
	output->SetRegions( region );
	output->SetOrigin( input->GetOrigin() );
	output->SetSpacing( input->GetSpacing() );
	output->SetDirection(input->GetDirection());
	output->CopyInformation( input );
	output->Allocate();

	// Typedefs
	typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType > BinaryMedianFilterType;

	typedef Image< unsigned char, 2 > ByteImageType2D;

	typedef itk::BinaryShapeOpeningImageFilter< ByteImageType2D > BinaryShapeOpeningImageFilter2D;

	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, BinaryShapeOpeningImageFilter2D > SliceBySliceImageFilterBackgroundType;

	typedef itk::ConnectedThresholdImageFilter<ByteImageType, ByteImageType> ConnectedThresholdImageFilterByteType;

	typedef itk::BinaryBallStructuringElement<unsigned char, InputImageDimension> StructuringElementType;

	typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType> DilateFilterType;

	typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType> ConnectedThresholdImageFilterType;

	typedef itk::IsolatedConnectedImageFilter< ImageType, ByteImageType > IsolatedConnectedImageFilterType;

	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType > BinaryShapeKeepNObjectsImageFilterType;
		
	
	// Find air components
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(spacing);
	air->SetDirection(input->GetDirection());
	air->Allocate();
	
	ByteIteratorType air_iter(air,region);

	for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter)
	{
		if (input_iter.Get() < -600)
		{
			air_iter.Set(255);
		} else {
			air_iter.Set(0);
		}
	}

	WriteITK(air, "air_lumen.nii" );

	// Apply median filter to remove noise
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( air );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 255 );

	ByteImageType::SizeType radius;
	radius.Fill(0);
	radius[0] = 4;
	radius[1] = 4;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	air = medianFilter->GetOutput();
	air_iter = ByteIteratorType(air,region);

	WriteITK(air, "air_median.nii");

	// Remove background air component
	// Remove components in 2D that are touching XY slice border [REQUIRES 3D INPUT]
	
	BinaryShapeOpeningImageFilter2D::Pointer bkgFilter2D = BinaryShapeOpeningImageFilter2D::New();
	bkgFilter2D->SetAttribute("SizeOnBorder");
	bkgFilter2D->SetBackgroundValue(0);
	bkgFilter2D->SetForegroundValue(255);
	bkgFilter2D->SetLambda(0);
	bkgFilter2D->SetReverseOrdering(true);
	
	SliceBySliceImageFilterBackgroundType::Pointer bkgRemover = SliceBySliceImageFilterBackgroundType::New();
	bkgRemover->SetInput( air );
	bkgRemover->SetFilter( bkgFilter2D );
	bkgRemover->Update();
	air = bkgRemover->GetOutput();
	air_iter = ByteIteratorType(air,region);

	WriteITK(air,"air_no_bkg.nii");	
	
	if ( m_RemoveBoneLung )
	{
		// Detect lungs using region growing with seeds from first slice
		// Assumption: only low intensity object at z=0 is lung
	
		ConnectedThresholdImageFilterByteType::Pointer lungGrower = ConnectedThresholdImageFilterByteType::New();
		lungGrower->SetInput( air );
		lungGrower->SetLower(255);
		lungGrower->SetUpper(255);
		lungGrower->SetReplaceValue(255);

		// Set air regions within dilated area as seeds of region growing
		for (air_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter)
		{
			ByteImageType::IndexType idx = air_iter.GetIndex();
		
			if (idx[2] == 0)
			{
				if (air_iter.Get() == 255)
				{
					lungGrower->AddSeed(idx);
				}
			
			}
		}

		lungGrower->Update();
		ByteImageType::Pointer lung = lungGrower->GetOutput();
		ByteIteratorType lung_iter(lung,region);

		WriteITK(lung,"lung.nii");
		
		// Remove lung from air mask
		for (air_iter.GoToBegin(), lung_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++lung_iter)
		{
			if (lung_iter.Get() == 255)
				air_iter.Set( 0 );
		}

		WriteITK(air,"air_no_lung.nii");
		lung.~SmartPointer();
		lungGrower.~SmartPointer();
	
	}
	
	// Dilate air

	radius.Fill(0);
	radius[0] = 3;
	radius[1] = 3;

	StructuringElementType ball;
	ball.SetRadius( radius );
	ball.CreateStructuringElement();	

	DilateFilterType::Pointer dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(255);
	dilater->SetInput( air );
	dilater->SetKernel( ball );
	dilater->Update();
	ByteImageType::Pointer air_dilated = dilater->GetOutput();
	ByteIteratorType air_dilated_iter(air_dilated,region);

	WriteITK(air_dilated, "air_dilated.nii");

	// Subtract original air area from dilation
	for (air_iter.GoToBegin(), air_dilated_iter.GoToBegin(); !air_iter.IsAtEnd(); ++air_iter, ++air_dilated_iter)
	{
		air_dilated_iter.Set( air_dilated_iter.Get() - air_iter.Get() );
	}

	WriteITK(air_dilated,"air_dilated_only.nii");

	// Find tagged regions, with or without bone

	ByteImageType::Pointer tagged;

	if ( m_RemoveBoneLung )
	{
		// Detect bone using region growing with seeds from first slice
		// Assumption: only high intensity object at z=0 is bone
		InputImagePixelType bone_threshold = 250;

		ConnectedThresholdImageFilterType::Pointer grower = ConnectedThresholdImageFilterType::New();
		grower->SetInput( input );
		grower->SetLower( bone_threshold );
		grower->SetUpper(1500); //1200
		grower->SetReplaceValue(255);

		// Set tagged regions within dilated area as seeds of region growing
		for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
		{
			InputImageType::IndexType idx = input_iter.GetIndex();
		
			if (idx[2] == 0)
			{
				if (input_iter.Get() > bone_threshold )
				{
					grower->AddSeed(idx);
				}
			}
		}

		grower->Update();
		ByteImageType::Pointer bone = grower->GetOutput();
		ByteIteratorType bone_iter(bone,region);

		WriteITK(bone,"bone.nii");
		
		IsolatedConnectedImageFilterType::Pointer isolatedGrower = IsolatedConnectedImageFilterType::New();
		isolatedGrower->SetInput( input );
		isolatedGrower->SetLower(m_TaggedValue);
		isolatedGrower->SetUpper( 1500 ); // 1200
		isolatedGrower->SetReplaceValue( 255 );
		isolatedGrower->FindUpperThresholdOff();

		// Set tagged regions within dilated area as seeds of region growing
		for (input_iter.GoToBegin(), air_dilated_iter.GoToBegin(), bone_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_dilated_iter, ++bone_iter)
		{

			if ( ! ( air_dilated_iter.Get() == 255 && bone_iter.Get() == 255 ) )
			{
				if (air_dilated_iter.Get() == 255 && input_iter.Get() >= m_TaggedValue)
				{
					isolatedGrower->AddSeed1( input_iter.GetIndex() );
				}

				if (bone_iter.Get() == 255)
				{
					isolatedGrower->AddSeed2( input_iter.GetIndex() );
				}
			}
		}

		isolatedGrower->Update();

		std::cout << "Bone / Stool Isolating Threshold: " << isolatedGrower->GetIsolatedValue() << std::endl;

		tagged = isolatedGrower->GetOutput();

		WriteITK(tagged,"tagged_without_bone.nii");	
	
	} else {

		// Setup region growing filter
		ConnectedThresholdImageFilterType::Pointer grower2 = ConnectedThresholdImageFilterType::New();
		grower2->SetInput( input );
		grower2->SetLower(m_TaggedValue);
		grower2->SetUpper(1500); //1200
		grower2->SetReplaceValue(255);

		// Set tagged regions within dilated area as seeds of region growing
		for (input_iter.GoToBegin(), air_dilated_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_dilated_iter)
		{
			if (air_dilated_iter.Get() == 255 && input_iter.Get() >= m_TaggedValue)
			{
				grower2->AddSeed( input_iter.GetIndex() );
			}
		}

		grower2->Update();
		tagged = grower2->GetOutput();
		
		
		WriteITK(tagged,"tagged.nii");

	}

	ByteIteratorType tagged_iter(tagged,region);
	

	// Add air lumen to tagged region
	for (tagged_iter.GoToBegin(), air_iter.GoToBegin(); !tagged_iter.IsAtEnd(); ++tagged_iter, ++air_iter)
	{
		if ( air_iter.Get() == 255 )
			tagged_iter.Set(255);
	}

	// Smooth tagged air via median filter
	medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( tagged );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 255 );

	radius.Fill(0);
	radius[0] = 1;
	radius[1] = 1;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	tagged = medianFilter->GetOutput();
	tagged_iter = ByteIteratorType(tagged,region);

	WriteITK(tagged,"tagged_air.nii");

	// Dilate
	dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(255);
	dilater->SetInput( tagged );
	
	radius.Fill(0);
	radius[0] = 8;
	radius[1] = 8;

	ball.SetRadius( radius );
	ball.CreateStructuringElement();

	dilater->SetKernel( ball );
	dilater->Update();

	WriteITK(dilater->GetOutput(),"colon.nii");
	
	ByteImageType::Pointer colon = dilater->GetOutput();
	
	// Cast to output image type
	IteratorType output_iter(output,region);
	ByteIteratorType colon_iter(colon,region);
	
	for (output_iter.GoToBegin(), colon_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++colon_iter)
	{
		if ( colon_iter.Get() == 255 )
		{
			output_iter.Set( m_OutputForegroundValue );
		} else {
			output_iter.Set( m_OutputBackgroundValue );
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
  os << indent << "OutputForegroundValue: " << m_OutputForegroundValue << std::endl;
  os << indent << "OutputBackgroundValue: " << m_OutputBackgroundValue << std::endl;
  os << indent << "PrintImages: " << m_PrintImages << std::endl;
  os << indent << "RemoveBoneLung: " << m_RemoveBoneLung << std::endl;
}

template <class TInputImage, class TOutputImage >
void
ColonSegmentationFilter<TInputImage,TOutputImage>::
WriteITK(typename ByteImageType::Pointer image, std::string name)
{
	if (m_PrintImages)
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

		
		
		typedef itk::ImageFileWriter< ByteImageType >  WriterType;
		WriterType::Pointer writer = WriterType::New();
		WriterType::SetGlobalWarningDisplay(false);
		writer->SetFileName(name.c_str());
		writer->SetInput(image);
		std::cout<<"Writing: "<<name<<std::endl;
		writer->Update();
	}
}








} // end namespace itk

#endif
