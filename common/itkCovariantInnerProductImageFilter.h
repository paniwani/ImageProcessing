/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCovariantInnerProductImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:19 $
  Version:   $Revision: 1.23 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCovariantInnerProductImageFilter_h
#define __itkCovariantInnerProductImageFilter_h

#include "itkBinaryFunctorImageFilter.h"

namespace itk
{
  
/** \class CovariantInnerProductImageFilter
 *
 * \brief Take a pair of images of vectors as input and produce an image with the
 *  dot product of those vectors.
 *
 * The filter expects the input image pixel type to be a vector and 
 * the output image pixel type to be a scalar.
 *
 * This filter assumes that the PixelType of the input image is a VectorType
 * that provides a * method.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */

namespace Functor {  
  
template< class TInput1, class TInput2, class TOutput>
class InnerProduct
{
public:
  InnerProduct() {}
  ~InnerProduct() {}

  bool operator!=( const InnerProduct & ) const
    {
    return false;
    }
  bool operator==( const InnerProduct & other ) const
    {
    return !(*this != other);
    }

  void SetScale(double i){m_Scale = i;}

  inline TOutput operator()( const TInput1 & A, const TInput2 &B ) const
    {
      return static_cast<TOutput>(m_Scale*( A * B ));
    }
private:
  double m_Scale;
}; 
}
 
template <class TInputImage1, class TInputImage2, class TOutputImage>
class ITK_EXPORT CovariantInnerProductImageFilter :
    public
BinaryFunctorImageFilter<TInputImage1, TInputImage2, TOutputImage, 
                        Functor::InnerProduct< typename TInputImage1::PixelType, 
					       typename TInputImage2::PixelType,
					       typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef CovariantInnerProductImageFilter Self;
  typedef BinaryFunctorImageFilter<
  TInputImage1, TInputImage2, TOutputImage, 
      Functor::InnerProduct< typename TInputImage1::PixelType, typename TInputImage2::PixelType,
    typename TOutputImage::PixelType> >
                                         Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;
    
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(CovariantInnerProductImageFilter, 
               BinaryFunctorImageFilter);

  // mostly for making the gradient negative, so we don't need a
  // separate filter to do it.
  void SetScale(double scale)
  {
    this->GetFunctor().SetScale(scale);
    this->Modified();
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  /** End concept checking */
#endif

protected:
  CovariantInnerProductImageFilter() {this->SetScale(1.0);}
  virtual ~CovariantInnerProductImageFilter() {}
    
private:
  CovariantInnerProductImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};
 
} // end namespace itk


#endif
