#ifndef __itkFastDilateFilter_h
#define __itkFastDilateFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
  
/** \class FastDilateFilter
 * \brief Implements uncleansed colon segmentation
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.  
 * 
 */
template <	class TInputImage, class TOutputImage=TInputImage>
			
class ITK_EXPORT FastDilateFilter : public ImageToImageFilter<TInputImage,TOutputImage> 
{
public:
  /** Standard class typedefs. */
  typedef FastDilateFilter                            Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(FastDilateFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType; 
  typedef typename    InputImageType::PixelType  InputImagePixelType; 

  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;  
  
  itkSetMacro(NumOfIterations, unsigned int);
  itkGetConstMacro(NumOfIterations, unsigned int);
  
  itkSetMacro(FullyConnected, bool);
  itkGetConstMacro(FullyConnected, bool);

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
  FastDilateFilter();
  virtual ~FastDilateFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method implements the actual reflection of the image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void GenerateData(void);

private:
  FastDilateFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_NumOfIterations;
  unsigned int m_FullyConnected;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastDilateFilter.txx"
#endif

#endif
