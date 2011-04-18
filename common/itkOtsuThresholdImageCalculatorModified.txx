/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkOtsuThresholdImageCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-26 21:45:54 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOtsuThresholdImageCalculator_txx
#define __itkOtsuThresholdImageCalculator_txx

#include "itkOtsuThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
OtsuThresholdImageCalculatorModified<TInputImage>
::OtsuThresholdImageCalculatorModified()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 0;
  m_RegionSetByUser = false;
  m_HistogramMin = 0;
  m_HistogramMax = 0;
  m_PrintHistogram="";
  m_MinMax = false;
}


/*
 * Compute the Otsu's threshold
 */
template<class TInputImage>
void
OtsuThresholdImageCalculatorModified<TInputImage>
::Compute(void)
{

  unsigned int j;

  if ( !m_Image ) { return; }
  if( !m_RegionSetByUser )
    {
    m_Region = m_Image->GetRequestedRegion();
    }

	double totalPixels = 0;
	
	// get total number of pixels within min and max
	typedef ImageRegionConstIteratorWithIndex<TInputImage> Iterator;
	Iterator iter( m_Image, m_Region );
	iter.GoToBegin();
	while ( !iter.IsAtEnd() )
	{
		if ( iter.Get() >= m_HistogramMin && iter.Get() <= m_HistogramMax )
		{
			totalPixels++;
		}
		++iter;
	}

	// default histogram values to min max of image
	// compute image max and min
	if ( !m_MinMax )
	{
		typedef MinimumMaximumImageCalculator<TInputImage> RangeCalculator;
		typename RangeCalculator::Pointer rangeCalculator = RangeCalculator::New();
		rangeCalculator->SetImage( m_Image );
		rangeCalculator->Compute();

		m_HistogramMin = rangeCalculator->GetMinimum();
		m_HistogramMax = rangeCalculator->GetMaximum();

		if ( m_HistogramMin >= m_HistogramMax )
		{
			m_Threshold = m_HistogramMin;
			return;
		}
	}
	
	//////////////////////////////////////////////////////////////////////
	// validate inputs and image
	if ( totalPixels == 0 ) { return; }
	
	// by default, use every pixel within min and max
	if (m_NumberOfHistogramBins == 0)
	{
		m_NumberOfHistogramBins = m_HistogramMax - m_HistogramMin;
	}

	std::cout << "Total number of pixels: " << totalPixels << std::endl;
	std::cout << "HistogramMin: " << m_HistogramMin << std::endl;
	std::cout << "HistogramMax: " << m_HistogramMax << std::endl;
	std::cout << "PrintHistogram: " << m_PrintHistogram << std::endl;
	
	/////////////////////////////////////////////////////////////////////
	
	
  PixelType imageMin = m_HistogramMin;
  PixelType imageMax = m_HistogramMax;


  // create a histogram
  std::vector<double> relativeFrequency;
  relativeFrequency.resize( m_NumberOfHistogramBins );
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    relativeFrequency[j] = 0.0;
    }

  double binMultiplier = (double) m_NumberOfHistogramBins /
    (double) ( imageMax - imageMin );

  iter.GoToBegin();

  while ( !iter.IsAtEnd() )
    {
		if ( iter.Get() >= m_HistogramMin && iter.Get() <= m_HistogramMax )
		{
			unsigned int binNumber;
			PixelType value = iter.Get();

			if ( value == imageMin ) 
			  {
			  binNumber = 0;
			  }
			else
			  {
			  binNumber = (unsigned int) vcl_ceil((value - imageMin) * binMultiplier ) - 1;
			  if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
				{
				binNumber -= 1;
				}
			  }

			relativeFrequency[binNumber] += 1.0;
		}
		
		++iter;

    }
 
  // normalize the frequencies
  double totalMean = 0.0;
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    relativeFrequency[j] /= totalPixels;
    totalMean += (j+1) * relativeFrequency[j];
    }
	
	
	// optional: print histogram
	if ( m_PrintHistogram != "" )
	{
		std::ofstream file;
		file.open( m_PrintHistogram.c_str() );
		file << "Bin,Frequency\n";
	
		for ( j = 0; j < m_NumberOfHistogramBins; j++ )
		{
			file << j+m_HistogramMin << "," << relativeFrequency[j] << "\n";
		}
		
		file.close();
	}


  // compute Otsu's threshold by maximizing the between-class
  // variance
  double freqLeft = relativeFrequency[0];
  double meanLeft = 1.0;
  double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

  double maxVarBetween = freqLeft * ( 1.0 - freqLeft ) *
    vnl_math_sqr( meanLeft - meanRight );
  int maxBinNumber = 0;

  double freqLeftOld = freqLeft;
  double meanLeftOld = meanLeft;

  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    freqLeft += relativeFrequency[j];
    meanLeft = ( meanLeftOld * freqLeftOld + 
                 (j+1) * relativeFrequency[j] ) / freqLeft;
    if (freqLeft == 1.0)
      {
      meanRight = 0.0;
      }
    else
      {
      meanRight = ( totalMean - meanLeft * freqLeft ) / 
        ( 1.0 - freqLeft );
      }
    double varBetween = freqLeft * ( 1.0 - freqLeft ) *
      vnl_math_sqr( meanLeft - meanRight );
   
    if ( varBetween > maxVarBetween )
      {
      maxVarBetween = varBetween;
      maxBinNumber = j;
      }

    // cache old values
    freqLeftOld = freqLeft;
    meanLeftOld = meanLeft; 

    } 

  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( maxBinNumber + 1 ) / binMultiplier );
}

template<class TInputImage>
void
OtsuThresholdImageCalculatorModified<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
OtsuThresholdImageCalculatorModified<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
  os << indent << "HistogramMin: " << m_HistogramMin << std::endl;
  os << indent << "HistogramMax: " << m_HistogramMax << std::endl;
  os << indent << "PrintHistogram: " << m_PrintHistogram << std::endl;
  os << indent << "MinMax: " << m_MinMax << std::endl;
}

} // end namespace itk

#endif
