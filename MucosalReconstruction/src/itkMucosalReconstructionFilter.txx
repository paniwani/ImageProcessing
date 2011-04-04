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
  m_NumOfLayers = 3;
}


/**
 * GenerateData Performs the reconstruction
 */
template <class TInputImage, class TOutputImage >
void
MucosalReconstructionFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
  typename Superclass::InputImageConstPointer  input = this->GetInput();
  typename Superclass::OutputImagePointer output = this->GetOutput(0);
  
  InputImageType::RegionType region = input->GetLargestPossibleRegion();

  output->SetRegions( region );
  output->Allocate();
  
  // Find air components
  
  typedef itk::Image< unsigned char, 2 >		ByteImageType;
  typedef itk::Image< int, 2>					IntImageType;
  
  ByteImageType::Pointer air = ByteImageType::New();
  air->SetRegions( region );
  air->Allocate();
  
  typedef ImageRegionConstIterator<TInputImage> InputIteratorType;
  typedef ImageRegionIterator<TOutputImage>     OutputIteratorType;
  typedef ImageRegionIterator<ByteImageType>    ByteIteratorType;
  
  InputIteratorType input_iter(input,region);
  ByteIteratorType air_iter(air,region);
  
  for (input_iter.GoToBegin(), air_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
  {
	if ( input_iter.Get() < -600 )
	{
		air_iter.Set( 1 );
	} else {
		air_iter.Set( 0 );
	}
  }
  
  // Median filter to air mask
 typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType > BinaryMedianFilterType;
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( air );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 1 );

	ByteImageType::SizeType radius;
	radius[0] = 4;
	radius[1] = 4;
	radius[2] = 0;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	air = medianFilter->GetOutput();
	air = ByteIteratorType( air, region);
  
  typedef ConnectedComponentImageFilter< ByteImageType, IntImageType > ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer ccFilter = ConnectedComponentImageFilterType::New();
  ccFilter->SetInput( air );
    
  // Remove external air by finding largest air component touching the image border
  
  typedef itk::ShapeRelabelImageFilter< IntImageType > ShapeRelabelImageFilterType;
  
	ShapeRelabelImageFilterType::Pointer relabeler = ShapeRelabelImageFilterType::New();
	
	relabeler->SetInput( ccFilter->GetOutput() );
	relabeler->SetAttribute("SizeOnBorder");
	relabeler->SetBackgroundValue(0);
	//relabeler->SetReverseOrdering(true);
	relabeler->Update();
	
	IntImageType::Pointer cc = relabeler->GetOutput();
	
	
  
  
  
  
  
  // Dilate
  // Reconstruct in direction of gradient
  
  
  
  
  
  
  

  typedef ImageLinearConstIteratorWithIndex<TInputImage>  InputIterator;
  typedef ImageLinearIteratorWithIndex<TOutputImage>      OutputIterator;

  InputIterator  inputIt(  inputPtr,  inputPtr->GetRequestedRegion() );
  OutputIterator outputIt( outputPtr, outputPtr->GetRequestedRegion() );

       
  inputIt.SetDirection( m_Direction );
  outputIt.SetDirection( m_Direction );

  inputIt.GoToBegin();
  outputIt.GoToBegin();

  while( !inputIt.IsAtEnd() ) 
    {

    outputIt.GoToEndOfLine();
    --outputIt;
    while( !inputIt.IsAtEndOfLine() ) 
      {
      outputIt.Set( inputIt.Get() );
      ++inputIt;
      --outputIt;
      progress.CompletedPixel();
      }

    inputIt.NextLine();
    outputIt.GoToEndOfLine(); // NextLine() assumes that the 
    outputIt.NextLine();      // iterator is at the end of line.
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
