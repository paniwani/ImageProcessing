#ifndef __itkDirectionalGradientImageFilter_h
#define __itkDirectionalGradientImageFilter_h

#include <itkCovariantVector.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include "itkCovariantInnerProductImageFilter.h"
#include <itkGradientImageFilter.h>
#include "itkMorphologicalDistanceTransformImageFilter.h"

namespace itk {
/**
 * \class DirectionalGradientImageFilter
 * \brief Compute the component of gradient orthogonal to a mask.
 *
 * The direction gradient is a way of computing components of image
 * gradient orthogonal to a mask contour. This is done by computing
 * the gradient of the distance transform of the mask and taking the
 * dot product with the gradient of the greyscale image. 
 * This is useful in situations where we are looking for the sign of
 * change of brightness as we move towards or away from some
 * object. For example, a cell wall might be bright on a dark
 * background, and we are interested in the outer edge. The outer edge
 * will be indicated by a decreasing brightness as we move away from
 * the cell centre, while the inner edge will be an increasing
 * edge. This filter allows the edges to be distinguised easily.
 *
 * This is a mini-pipeline filter that uses the recursive gradient
 * filters, morphological distance transform filters and gradient
 * filters.
 *
 * Warning: The mask must have values of 0 and 1 for the distance transform to
 * work correctly (0 and 255, for example, won't work). This is
 * because the morphological distrance transform is used
 * internally. It is very fast, but hasn't been set up to work on
 * arbitary mask values.
 *
 * \author Richard Beare, Monash University, Department of Medicine
*/


template<class TInputImage, class TMaskImage, class TOutputImage>
class ITK_EXPORT DirectionalGradientImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DirectionalGradientImageFilter         Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef TInputImage  InputImageType;
  typedef TMaskImage   MaskImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::PixelType InputPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DirectionalGradientImageFilter, ImageToImageFilter);

  void SetMaskImage(TMaskImage *input)
  {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<TMaskImage *>(input) );
  }
  MaskImageType * GetMaskImage()
  {
    return static_cast<MaskImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(1)));

  }

  /**
 * Set/Get the Outside value used in the distance transform. Adjust
 * this depending on whether you want the gradient to be computed
 * inside or outside your mask.
**/
  itkSetMacro(OutsideValue, InputPixelType);
  itkGetConstReferenceMacro(OutsideValue, InputPixelType);

  /** 
* Set/Get the smoothing parameter in the gradient computation
**/

  itkSetMacro(Sigma, double);
  itkGetConstReferenceMacro(Sigma, double);

/** 
* A scaling factor applied during the dot product phase. Useful for
*  flipping the sign of the output 
**/

  itkSetMacro(Scale, double);
  itkGetConstReferenceMacro(Scale, double);

protected:
  DirectionalGradientImageFilter();
  ~DirectionalGradientImageFilter() {};

  void GenerateData();
private:
  double m_Sigma;
  double m_Scale;
  InputPixelType m_OutsideValue;
  typedef typename itk::Image<float, TInputImage::ImageDimension> FImType;
  typedef typename itk::Image<itk::CovariantVector<float, TInputImage::ImageDimension>, TInputImage::ImageDimension  > GradImType;
  typedef typename itk::GradientRecursiveGaussianImageFilter<InputImageType , GradImType >  GaussGradFiltType;
  typedef typename itk::CovariantInnerProductImageFilter< GradImType, GradImType, OutputImageType  > InnerProductType;

  typedef typename itk::MorphologicalDistanceTransformImageFilter<MaskImageType, FImType > DistTransType;
  // don't need a fancy gradient for the distance transform
  typedef typename itk::GradientImageFilter<FImType> GradFiltType;


  typename DistTransType::Pointer m_DT;
  typename GradFiltType::Pointer m_GradDT;
  typename InnerProductType::Pointer m_Innerprod;
  typename GaussGradFiltType::Pointer m_RawGrad;


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDirectionalGradientImageFilter.txx"
#endif


#endif
