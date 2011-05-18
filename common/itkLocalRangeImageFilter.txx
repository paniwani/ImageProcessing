/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalRangeImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 16:30:50 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLocalRangeImageFilter_txx
#define __itkLocalRangeImageFilter_txx

#include "itkLocalRangeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkProgressReporter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage >
LocalRangeImageFilter<TInputImage,TOutputImage >
::LocalRangeImageFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_Size.Fill(0);
}


/**
 * GenerateData Performs the reflection
 */
template <class TInputImage, class TOutputImage >
void
LocalRangeImageFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
  typename Superclass::InputImageConstPointer  inputPtr = this->GetInput();
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput(0);

  outputPtr->SetRequestedRegion( inputPtr->GetRequestedRegion() );
  outputPtr->SetBufferedRegion( inputPtr->GetBufferedRegion() );
  outputPtr->SetLargestPossibleRegion( inputPtr->GetLargestPossibleRegion() );
  outputPtr->Allocate();

  typedef ConstNeighborhoodIterator<InputImageType> InputIteratorType;
  typedef ImageRegionIteratorWithIndex<TOutputImage> OutputIteratorType;
  
  InputIteratorType  it(  m_Size, inputPtr, inputPtr->GetRequestedRegion() );
  OutputIteratorType oit( outputPtr, outputPtr->GetRequestedRegion() );

  ProgressReporter progress(this, 0,  inputPtr->GetRequestedRegion().GetNumberOfPixels() );
  
  for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
  {
	InputImagePixelType min = NumericTraits<InputImagePixelType>::max();
	InputImagePixelType max = NumericTraits<InputImagePixelType>::NonpositiveMin();
	
	for (int i=0; i < it.Size(); i++)
	{
		InputImagePixelType val = it.GetPixel(i);
		
		if (val > max)
			max = val;
			
		if (val < min)
			min = val;
	}
	
	oit.Set( max - min);
	
	progress.CompletedPixel();
  }
}

template <class TInputImage, class TOutputImage >
void
LocalRangeImageFilter<TInputImage,TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "Size: " << m_Size << std::endl;
}

} // end namespace itk

#endif
