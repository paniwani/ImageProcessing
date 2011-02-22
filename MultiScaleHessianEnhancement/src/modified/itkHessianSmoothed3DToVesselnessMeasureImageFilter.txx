/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHessianSmoothed3DToVesselnessMeasureImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 20:59:44 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHessianSmoothed3DToVesselnessMeasureImageFilter_txx
#define __itkHessianSmoothed3DToVesselnessMeasureImageFilter_txx

#include "itkHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

#define EPSILON  1e-03

namespace itk
{

/**
 * Constructor
 */
template < typename TPixel >
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::HessianSmoothed3DToVesselnessMeasureImageFilter()
{
  m_Alpha = 0.5;
  m_Beta  = 0.3;
  m_Gamma = 0.3;
  m_Eta   = 0.2;

  m_C = 10e-6;

  m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  m_SymmetricEigenValueFilter->SetDimension( ImageDimension );
  // m_SymmetricEigenValueFilter->OrderEigenValuesBy( EigenAnalysisFilterType::FunctorType::OrderByValue );
  m_SymmetricEigenValueFilter->OrderEigenValuesBy( EigenAnalysisFilterType::FunctorType::OrderByMagnitude );
  
  // By default, scale the vesselness measure by the largest eigenvalue
  m_ScaleVesselnessMeasure  = false;
  
}


template < typename TPixel >
void 
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::GenerateData()
{
  itkDebugMacro(
      << "HessianSmoothed3DToVesselnessMeasureImageFilter generating data ");

  m_SymmetricEigenValueFilter->SetInput( this->GetInput() );
  
  typename OutputImageType::Pointer output = this->GetOutput();

  typedef typename EigenAnalysisFilterType::OutputImageType
                                            EigenValueImageType;

  m_SymmetricEigenValueFilter->Update();
  
  const typename EigenValueImageType::ConstPointer eigenImage = 
                    m_SymmetricEigenValueFilter->GetOutput();
  
  // walk the region of eigen values and get the vesselness measure
  EigenValueArrayType eigenValue;
  ImageRegionConstIterator<EigenValueImageType> it;
  it = ImageRegionConstIterator<EigenValueImageType>(
      eigenImage, eigenImage->GetRequestedRegion());
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());
  oit.GoToBegin();
  it.GoToBegin();
  while (!it.IsAtEnd())
    {
    // Get the eigen value
    eigenValue = it.Get();


    // Find the smallest eigenvalue
    double smallest = vnl_math_abs( eigenValue[0] );
    double Lambda1 = eigenValue[0];
 
    for ( unsigned int i=1; i <=2; i++ )
      {
      if ( vnl_math_abs( eigenValue[i] ) < smallest )
        {
        Lambda1 = eigenValue[i];
        smallest = vnl_math_abs( eigenValue[i] );
        }
      }

    // Find the largest eigenvalue
    double largest = vnl_math_abs( eigenValue[0] );
    double Lambda3 = eigenValue[0];
 
    for ( unsigned int i=1; i <=2; i++ )
      {
      if (  vnl_math_abs( eigenValue[i] ) > largest  )
        {
        Lambda3 = eigenValue[i];
        largest = vnl_math_abs( eigenValue[i] );
        }
      }


    //  find Lambda2 so that |Lambda1| < |Lambda2| < |Lambda3|
    double Lambda2 = eigenValue[0];

    for ( unsigned int i=0; i <=2; i++ )
      {
      if ( eigenValue[i] != Lambda1 && eigenValue[i] != Lambda3 )
        {
        Lambda2 = eigenValue[i];
        break;
        }
      }

    if ( Lambda2 >= 0.0 ||  Lambda3 >= 0.0 || 
         vnl_math_abs( Lambda2) < EPSILON  || 
         vnl_math_abs( Lambda3 ) < EPSILON )
      {
      oit.Set( NumericTraits< OutputPixelType >::Zero );
      } 
    else
      {
   
      double Lambda1Abs = vnl_math_abs( Lambda1 );
      double Lambda2Abs = vnl_math_abs( Lambda2 );
      double Lambda3Abs = vnl_math_abs( Lambda3 );

	  double Ra = Lambda1Abs / vcl_sqrt( Lambda2Abs * Lambda3Abs );
	  double Fa = vcl_exp(  - vnl_math_sqr( Ra ) / (2 * vnl_math_sqr( m_Alpha ) ) );

	  double Rb = Lambda2Abs / Lambda3Abs;
	  double Fb = vcl_exp( - vnl_math_sqr( Rb - m_Gamma ) / (2 * vnl_math_sqr( m_Beta ) ) );

	  double Frut = Fa * Fb;

	  double Fcup = (1.0 - vcl_exp( - vnl_math_sqr( Lambda1Abs / Lambda2Abs ) / (2 * vnl_math_sqr( m_Eta ) ) ) ) * (1.0 - vcl_exp( - vnl_math_sqr( Lambda2Abs / Lambda3Abs ) / (2 * vnl_math_sqr( m_Eta ) ) ) );

	  double H = vnl_math_max( Frut, Fcup );

      if(  m_ScaleVesselnessMeasure ) 
        {
        oit.Set( static_cast< OutputPixelType >(
                                     Lambda3Abs*H ) );
        }
      else
        {
        oit.Set( static_cast< OutputPixelType >( H ) );
        }
      }
    ++it;
    ++oit;
    }
    
}

template < typename TPixel >
void
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta:  " << m_Beta  << std::endl;
  os << indent << "Gamma: " << m_Gamma << std::endl;
  os << indent << "Eta: "   << m_Eta   << std::endl;

  os << indent << "C: " << m_C << std::endl;
}


} // end namespace itk
  
#endif
