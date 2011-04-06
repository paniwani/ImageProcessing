#ifndef __itkMucosalReconstructionFilter_h
#define __itkMucosalReconstructionFilter_h

#include "itkImageToImageFilter.h"
#include <itkBinaryContourImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <itkGaussianBlurImageFunction.h>

namespace itk
{
  
/** \class MucosalReconstructionFilter
 * \brief Implements a mucosal reconstruction post electronic subtraction
 *
 * This class is parameterized over the type of the input image and
 * the type of the output image.  
 * 
 */
template <class TInputImage, class TOutputImage, class TMaskImage>
class ITK_EXPORT MucosalReconstructionFilter : public ImageToImageFilter<TInputImage,TOutputImage> 
{
public:
  /** Standard class typedefs. */
  typedef MucosalReconstructionFilter                            Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MucosalReconstructionFilter, ImageToImageFilter);

  /** Some convenient typedefs. */
  typedef TInputImage                            InputImageType;
  typedef typename    InputImageType::Pointer    InputImagePointer;
  typedef typename    InputImageType::RegionType InputImageRegionType; 
  typedef typename    InputImageType::PixelType  InputImagePixelType; 

  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;  
  
  typedef TMaskImage							MaskImageType;
  typedef typename MaskImageType::Pointer   	MaskImagePointer;
  

  /** Set number of layers to reconstruct */
  itkGetConstMacro( NumOfLayers, unsigned int );
  itkSetMacro( NumOfLayers, unsigned int );

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
					  
	void SetMaskImage(TMaskImage* mask)
	{
	this->SetNthInput(1, const_cast<TMaskImage *>( mask ));
	}
					  
	const TMaskImage* GetMaskImage() const
	{
	return (static_cast<const TMaskImage*>(this->ProcessObject::GetInput(1)));
	}
					  

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<InputImagePixelType, OutputImagePixelType>));
  /** End concept checking */
#endif

protected:
  MucosalReconstructionFilter();
  virtual ~MucosalReconstructionFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method implements the actual reflection of the image.
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */
  void GenerateData(void);

private:
  MucosalReconstructionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_NumOfLayers;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMucosalReconstructionFilter.txx"
#endif

#endif
