/*=========================================================================

  Program:   ORFEO Toolbox
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See OTBCopyright.txt for details.


     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __otbScalarImageToTexturesMaskFilter_txx
#define __otbScalarImageToTexturesMaskFilter_txx

#include "otbScalarImageToTexturesMaskFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace otb
{
template <class TInputImage, class TOutputImage, class TMaskImage>
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::ScalarImageToTexturesMaskFilter() : m_Radius(),
  m_Offset(),
  m_NumberOfBinsPerAxis(8),
  m_InputImageMinimum(0),
  m_InputImageMaximum(256)
{
  // There are 8 outputs corresponding to the 8 textures indices
  this->SetNumberOfOutputs(8);

  // Create the 8 outputs
  this->SetNthOutput(0, OutputImageType::New());
  this->SetNthOutput(1, OutputImageType::New());
  this->SetNthOutput(2, OutputImageType::New());
  this->SetNthOutput(3, OutputImageType::New());
  this->SetNthOutput(4, OutputImageType::New());
  this->SetNthOutput(5, OutputImageType::New());
  this->SetNthOutput(6, OutputImageType::New());
  this->SetNthOutput(7, OutputImageType::New());
}

template <class TInputImage, class TOutputImage, class TMaskImage>
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::~ScalarImageToTexturesMaskFilter()
{}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetEnergyOutput()
{
  if (this->GetNumberOfOutputs() < 1)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(0));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetEntropyOutput()
{
  if (this->GetNumberOfOutputs() < 2)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(1));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetCorrelationOutput()
{
  if (this->GetNumberOfOutputs() < 3)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(2));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetInverseDifferenceMomentOutput()
{
  if (this->GetNumberOfOutputs() < 4)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(3));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetInertiaOutput()
{
  if (this->GetNumberOfOutputs() < 5)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(4));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetClusterShadeOutput()
{
  if (this->GetNumberOfOutputs() < 6)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(5));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetClusterProminenceOutput()
{
  if (this->GetNumberOfOutputs() < 7)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(6));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
typename ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::OutputImageType *
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GetHaralickCorrelationOutput()
{
  if (this->GetNumberOfOutputs() < 8)
    {
    return 0;
    }
  return static_cast<OutputImageType *>(this->GetOutput(7));
}

template <class TInputImage, class TOutputImage, class TMaskImage>
void
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::GenerateInputRequestedRegion()
{
  // First, call superclass implementation
  Superclass::GenerateInputRequestedRegion();

  // Retrieve the input and output pointers
  InputImagePointerType  inputPtr = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointerType outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
    {
    return;
    }

  // Retrieve the output requested region
  // We use only the first output since requested regions for all outputs are enforced to be equal
  // by the default GenerateOutputRequestedRegiont() implementation
  OutputRegionType outputRequestedRegion = outputPtr->GetRequestedRegion();

  typename OutputRegionType::IndexType outputIndex = outputRequestedRegion.GetIndex();
  typename OutputRegionType::SizeType  outputSize   = outputRequestedRegion.GetSize();
  typename InputRegionType::IndexType  inputIndex;
  typename InputRegionType::SizeType   inputSize;

  // First, apply offset
  for (unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim)
    {
    inputIndex[dim] = std::min(outputIndex[dim], outputIndex[dim] + m_Offset[dim]);
    inputSize[dim] =
      std::max(outputIndex[dim] + outputSize[dim], outputIndex[dim] + outputSize[dim] +
               m_Offset[dim]) - inputIndex[dim];
    }

  // Build the input requested region
  InputRegionType inputRequestedRegion;
  inputRequestedRegion.SetIndex(inputIndex);
  inputRequestedRegion.SetSize(inputSize);

  // Apply the radius
  inputRequestedRegion.PadByRadius(m_Radius);

  // Try to apply the requested region to the input image
  if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
    {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    }
  else
    {
    // Build an exception
    itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
	
  MaskImagePointer mask = const_cast<MaskImageType *>(this->GetMaskImage());
  if (mask)
    {
    mask->SetRequestedRegion( inputRequestedRegion );
    }
}

template <class TInputImage, class TOutputImage, class TMaskImage>
void
ScalarImageToTexturesMaskFilter<TInputImage, TOutputImage, TMaskImage>
::ThreadedGenerateData(const OutputRegionType& outputRegionForThread, int threadId)
{
  typename TMaskImage::ConstPointer mask = this->GetMaskImage();
  itk::ImageRegionConstIteratorWithIndex<MaskImageType> maskIt(mask, outputRegionForThread);

  // Retrieve the input and output pointers
  InputImagePointerType  inputPtr             =      const_cast<InputImageType *>(this->GetInput());
  OutputImagePointerType energyPtr            =      this->GetEnergyOutput();
  OutputImagePointerType entropyPtr           =      this->GetEntropyOutput();
  OutputImagePointerType correlationPtr       =      this->GetCorrelationOutput();
  OutputImagePointerType invDiffMomentPtr     =      this->GetInverseDifferenceMomentOutput();
  OutputImagePointerType inertiaPtr           =      this->GetInertiaOutput();
  OutputImagePointerType clusterShadePtr      =      this->GetClusterShadeOutput();
  OutputImagePointerType clusterProminencePtr =      this->GetClusterProminenceOutput();
  OutputImagePointerType haralickCorPtr       =      this->GetHaralickCorrelationOutput();

  // Fill default output buffers
  for (int i=0; i<8; i++)
  {
	OutputImagePointerType outPtr = this->GetOutput(i);
	outPtr->FillBuffer(0);
  }

  // Build output iterators
  itk::ImageRegionIteratorWithIndex<OutputImageType> energyIt(energyPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          entropyIt(entropyPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          correlationIt(correlationPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          invDiffMomentIt(invDiffMomentPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          inertiaIt(inertiaPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          clusterShadeIt(clusterShadePtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          clusterProminenceIt(clusterProminencePtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          haralickCorIt(haralickCorPtr, outputRegionForThread);

  // Go to begin
  maskIt.GoToBegin();
  energyIt.GoToBegin();
  entropyIt.GoToBegin();
  correlationIt.GoToBegin();
  invDiffMomentIt.GoToBegin();
  inertiaIt.GoToBegin();
  clusterShadeIt.GoToBegin();
  clusterProminenceIt.GoToBegin();
  haralickCorIt.GoToBegin();

  // Build the co-occurence matrix generator
  CoocurrenceMatrixGeneratorPointerType coOccurenceMatrixGenerator = CoocurrenceMatrixGeneratorType::New();
  coOccurenceMatrixGenerator->SetInput(inputPtr);
  coOccurenceMatrixGenerator->SetOffset(m_Offset);
  coOccurenceMatrixGenerator->SetNumberOfBinsPerAxis(m_NumberOfBinsPerAxis);
  coOccurenceMatrixGenerator->SetPixelValueMinMax(m_InputImageMinimum, m_InputImageMaximum);

  // Build the texture calculator
  TextureCoefficientsCalculatorPointerType texturesCalculator = TextureCoefficientsCalculatorType::New();

  // Set-up progress reporting
  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Iterate on outputs to compute textures
  while (!energyIt.IsAtEnd()
         && !entropyIt.IsAtEnd()
         && !correlationIt.IsAtEnd()
         && !invDiffMomentIt.IsAtEnd()
         && !inertiaIt.IsAtEnd()
         && !clusterShadeIt.IsAtEnd()
         && !clusterProminenceIt.IsAtEnd()
         && !haralickCorIt.IsAtEnd())
    {
		if (maskIt.Get() != 0)
		{
			// Compute the region on which co-occurence will be estimated
			typename InputRegionType::IndexType inputIndex = energyIt.GetIndex() - m_Radius;
			typename InputRegionType::SizeType  inputSize;

			// First, apply offset
			for (unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim)
			  {
			  inputSize[dim] = 2 * m_Radius[dim] + 1;
			  }

			// Build the input  region
			InputRegionType inputRegion;
			inputRegion.SetIndex(inputIndex);
			inputRegion.SetSize(inputSize);

			// Compute the co-occurence matrix
			coOccurenceMatrixGenerator->SetRegion(inputRegion);
			coOccurenceMatrixGenerator->Update();

			// Compute textures indices
			texturesCalculator->SetInput(coOccurenceMatrixGenerator->GetOutput());
			texturesCalculator->Update();

			// Fill outputs
			energyIt.Set(texturesCalculator->GetEnergy());
			entropyIt.Set(texturesCalculator->GetEntropy());
			correlationIt.Set(texturesCalculator->GetCorrelation());
			invDiffMomentIt.Set(texturesCalculator->GetInverseDifferenceMoment());
			inertiaIt.Set(texturesCalculator->GetInertia());
			clusterShadeIt.Set(texturesCalculator->GetClusterShade());
			clusterProminenceIt.Set(texturesCalculator->GetClusterProminence());
			haralickCorIt.Set(texturesCalculator->GetHaralickCorrelation());
		}

    // Update progress
    progress.CompletedPixel();

    // Increment iterators
	++maskIt;
    ++energyIt;
    ++entropyIt;
    ++correlationIt;
    ++invDiffMomentIt;
    ++inertiaIt;
    ++clusterShadeIt;
    ++clusterProminenceIt;
    ++haralickCorIt;
    }

}

} // End namespace otb

#endif
