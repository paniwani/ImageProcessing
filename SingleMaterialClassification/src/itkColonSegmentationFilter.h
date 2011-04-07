#ifndef __itkColonSegmentationFilter_h
#define __itkColonSegmentationFilter_h

#include "itkImageToImageFilter.h"
#include <itkBinaryContourImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>

namespace itk
{
  
/** \class ColonSegmentationFilter
 * \brief Implements uncleansed colon segmentation
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.  
 * 
 */
template <	class TInputImage, 
			class TOutputImage=Image< unsigned char, ::itk::GetImageDimension<TInputImage>::ImageDimension >  >
			
class ITK_EXPORT ColonSegmentationFilter : public ImageToImageFilter<TInputImage,TOutputImage> 
{
public:
  /** Standard class typedefs. */
  typedef ColonSegmentationFilter                            Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ColonSegmentationFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType; 
  typedef typename    InputImageType::PixelType  InputImagePixelType; 

  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;  
  

  /** Set tagged value */
  itkGetConstMacro( TaggedValue, unsigned int );
  itkSetMacro( TaggedValue, unsigned int );
  
  itkSetMacro(ForegroundValue, OutputImagePixelType);
  itkGetConstMacro(ForegroundValue, OutputImagePixelType);

  itkSetMacro(BackgroundValue, OutputImagePixelType);
  itkGetConstMacro(BackgroundValue, OutputImagePixelType);

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
  ColonSegmentationFilter();
  virtual ~ColonSegmentationFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method implements the actual reflection of the image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void GenerateData(void);

private:
  ColonSegmentationFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_TaggedValue;
  OutputImagePixelType m_ForegroundValue;
  OutputImagePixelType m_BackgroundValue;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkColonSegmentationFilter.txx"
#endif

#endif
