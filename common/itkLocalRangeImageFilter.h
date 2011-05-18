/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalRangeImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:28:02 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLocalRangeImageFilter_h
#define __itkLocalRangeImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
  
/** \class LocalRangeImageFilter
 * \brief Implements a Reflection of an image along a selected direction.
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.  
 * 
 * \ingroup   IntensityImageFilters     Singlethreaded
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT LocalRangeImageFilter : public ImageToImageFilter<TInputImage,TOutputImage> 
{
public:
  /** Standard class typedefs. */
  typedef LocalRangeImageFilter                            Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(LocalRangeImageFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType; 
  typedef typename    InputImageType::PixelType  InputImagePixelType; 

  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;

  /** Set size used for local neighborhood */  
  itkGetConstMacro( Size, typename InputImageType::SizeType);
  itkSetMacro( Size, typename InputImageType::SizeType);

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<InputImagePixelType, OutputImagePixelType>));
  /** End concept checking */
#endif

protected:
  LocalRangeImageFilter();
  virtual ~LocalRangeImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method implements the actual reflection of the image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void GenerateData(void);

private:
  LocalRangeImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename InputImageType::SizeType m_Size;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLocalRangeImageFilter.txx"
#endif

#endif
