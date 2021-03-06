#ifndef __itkDirectionalGradientImageFilter2_txx
#define __itkDirectionalGradientImageFilter2_txx

#include "itkDirectionalGradientImageFilter2.h"
#include "itkProgressAccumulator.h"
namespace itk
{
template<class TInputImage, class TMaskImage, class TOutputImage>
DirectionalGradientImageFilter2<TInputImage, TMaskImage, TOutputImage>
::DirectionalGradientImageFilter2()
{
  m_DT = DistTransType::New();
  m_GradDT = GradFiltType::New();
  m_Innerprod = InnerProductType::New();
  m_RawGrad = GaussGradFiltType::New();
  m_Scale=1.0;
  m_OutsideValue = 0;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
void
DirectionalGradientImageFilter2<TInputImage, TMaskImage, TOutputImage>
::GenerateData()
{

  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typename InputImageType::ConstPointer input  = this->GetInput();
  typename MaskImageType::ConstPointer mask = this->GetMaskImage();

  // construct mini pipeline
  m_DT->SetInput(mask);
  m_DT->SetOutsideValue(m_OutsideValue);
  m_GradDT->SetInput(m_DT->GetOutput());

  m_RawGrad->SetInput(input);

  m_Innerprod->SetInput(m_GradDT->GetOutput());
  m_Innerprod->SetInput2(m_RawGrad->GetOutput());
  m_Innerprod->SetScale(m_Scale);

  progress->RegisterInternalFilter(m_DT, 0.2f);
  progress->RegisterInternalFilter(m_GradDT, 0.2f);
  progress->RegisterInternalFilter(m_RawGrad, 0.4f);
  progress->RegisterInternalFilter(m_Innerprod, 0.2f);
  
  m_Innerprod->GraftOutput(this->GetOutput());
  m_Innerprod->Update();
  this->GraftOutput(m_Innerprod->GetOutput());
}


}

#endif
