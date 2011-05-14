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
#ifndef __otbScalarImageToTexturesFilter_txx
#define __otbScalarImageToTexturesFilter_txx

#include "otbScalarImageToTexturesFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace otb
{
template <class TInputImage, class TOutputImage, class TMaskImage>
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
::ScalarImageToTexturesFilter() : m_Radius(),
  m_Offset(),
  m_NumberOfBinsPerAxis(8),
  m_InputImageMinimum(0),
  m_InputImageMaximum(256),
  m_Feature("Energy")
{
  // There are 8 outputs corresponding to the 8 textures indices
  this->SetNumberOfOutputs(1);

  this->SetNthOutput(0, OutputImageType::New());

  // Create the 8 outputs
  //this->SetNthOutput(0, OutputImageType::New());
 /* this->SetNthOutput(1, OutputImageType::New());
  this->SetNthOutput(2, OutputImageType::New());
  this->SetNthOutput(3, OutputImageType::New());
  this->SetNthOutput(4, OutputImageType::New());
  this->SetNthOutput(5, OutputImageType::New());
  this->SetNthOutput(6, OutputImageType::New());
  this->SetNthOutput(7, OutputImageType::New());*/
}

template <class TInputImage, class TOutputImage, class TMaskImage>
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
::~ScalarImageToTexturesFilter()
{}

//template <class TInputImage, class TOutputImage, class TMaskImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
//::GetEnergyOutput()
//{
//  if (this->GetNumberOfOutputs() < 1)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(0));
//}

//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetEntropyOutput()
//{
//  if (this->GetNumberOfOutputs() < 2)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(1));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetCorrelationOutput()
//{
//  if (this->GetNumberOfOutputs() < 3)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(2));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetInverseDifferenceMomentOutput()
//{
//  if (this->GetNumberOfOutputs() < 4)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(3));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetInertiaOutput()
//{
//  if (this->GetNumberOfOutputs() < 5)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(4));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetClusterShadeOutput()
//{
//  if (this->GetNumberOfOutputs() < 6)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(5));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetClusterProminenceOutput()
//{
//  if (this->GetNumberOfOutputs() < 7)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(6));
//}
//
//template <class TInputImage, class TOutputImage>
//typename ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::OutputImageType *
//ScalarImageToTexturesFilter<TInputImage, TOutputImage>
//::GetHaralickCorrelationOutput()
//{
//  if (this->GetNumberOfOutputs() < 8)
//    {
//    return 0;
//    }
//  return static_cast<OutputImageType *>(this->GetOutput(7));
//}

// Validate feature input
template <class TInputImage, class TOutputImage, class TMaskImage>
void
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
::SetFeature(std::string f)
{
	if (f == "Energy" || f == "Entropy" || f == "Correlation" || f == "InverseDifferenceMoment" || f == "Intertia" || f == "ClusterShade" || f == "ClusterProminence" || f == "HaralickCorrelation")
	{
		m_Feature = f;
	} else {
		std::cout << "Feature not recognized. Default feature (Energy) set.";
		m_Feature = "Energy";
	}
}

template <class TInputImage, class TOutputImage, class TMaskImage>
void
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
::SetFeature(unsigned int i)
{
	if (i >=0 && i<8)
	{
		if		( i==0)		{ m_Feature = "Energy"; }
		else if ( i==1 )	{ m_Feature = "Correlation"; }
		else if ( i==2 )	{ m_Feature = "InverseDifferenceMoment"; }
		else if ( i==3 )	{ m_Feature = "Intertia"; }
		else if ( i==4 )	{ m_Feature = "ClusterShade"; }
		else if ( i==5 )	{ m_Feature = "Correlation"; }
		else if ( i==6 )	{ m_Feature = "ClusterProminence"; }
		else if ( i==7 )	{ m_Feature = "HaralickCorrelation"; }
	} else {
		std::cout << "Feature not recognized. Default feature (Energy) set.";
		m_Feature = "Energy";
	}
}

template <class TInputImage, class TOutputImage, class TMaskImage>
void
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
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
ScalarImageToTexturesFilter<TInputImage, TOutputImage, TMaskImage>
::ThreadedGenerateData(const OutputRegionType& outputRegionForThread, int threadId)
{
	typename TMaskImage::ConstPointer mask = this->GetMaskImage();
	itk::ImageRegionConstIteratorWithIndex<MaskImageType> maskIt(mask, outputRegionForThread);

  // Retrieve the input and output pointers
  InputImagePointerType  inputPtr             =      const_cast<InputImageType *>(this->GetInput());
  OutputImagePointerType outputPtr			  =		 this->GetOutput();
  
  /*OutputImagePointerType energyPtr            =      this->GetEnergyOutput();
  OutputImagePointerType entropyPtr           =      this->GetEntropyOutput();
  OutputImagePointerType correlationPtr       =      this->GetCorrelationOutput();
  OutputImagePointerType invDiffMomentPtr     =      this->GetInverseDifferenceMomentOutput();
  OutputImagePointerType inertiaPtr           =      this->GetInertiaOutput();
  OutputImagePointerType clusterShadePtr      =      this->GetClusterShadeOutput();
  OutputImagePointerType clusterProminencePtr =      this->GetClusterProminenceOutput();
  OutputImagePointerType haralickCorPtr       =      this->GetHaralickCorrelationOutput();*/

  // Build output iterators
  itk::ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread);
  /*itk::ImageRegionIteratorWithIndex<OutputImageType> energyIt(energyPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          entropyIt(entropyPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          correlationIt(correlationPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          invDiffMomentIt(invDiffMomentPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          inertiaIt(inertiaPtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          clusterShadeIt(clusterShadePtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          clusterProminenceIt(clusterProminencePtr, outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>          haralickCorIt(haralickCorPtr, outputRegionForThread);*/

  // Go to begin
  maskIt.GoToBegin();
  outputIt.GoToBegin();
  /*entropyIt.GoToBegin();
  correlationIt.GoToBegin();
  invDiffMomentIt.GoToBegin();
  inertiaIt.GoToBegin();
  clusterShadeIt.GoToBegin();
  clusterProminenceIt.GoToBegin();
  haralickCorIt.GoToBegin();*/

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
  while (!maskIt.IsAtEnd()
	     && !outputIt.IsAtEnd() 
		/*&& !energyIt.IsAtEnd()
         && !entropyIt.IsAtEnd()
         && !correlationIt.IsAtEnd()
         && !invDiffMomentIt.IsAtEnd()
         && !inertiaIt.IsAtEnd()
         && !clusterShadeIt.IsAtEnd()
         && !clusterProminenceIt.IsAtEnd()
         && !haralickCorIt.IsAtEnd()*/)
    {
		if (maskIt.Get() != 0)
		{
			// Compute the region on which co-occurence will be estimated
			typename InputRegionType::IndexType inputIndex = outputIt.GetIndex() - m_Radius;
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
			if      (m_Feature == "Energy")						{ outputIt.Set( texturesCalculator->GetEnergy() ); }
			else if (m_Feature == "Entropy")					{ outputIt.Set( texturesCalculator->GetEntropy() ); }
			else if (m_Feature == "Correlation")				{ outputIt.Set( texturesCalculator->GetCorrelation() ); }
			else if (m_Feature == "InverseDifferenceMoment")	{ outputIt.Set( texturesCalculator->GetInverseDifferenceMoment() ); }
			else if (m_Feature == "Intertia")					{ outputIt.Set( texturesCalculator->GetInertia() ); }
			else if (m_Feature == "ClusterShade")				{ outputIt.Set( texturesCalculator->GetClusterShade() ); }
			else if (m_Feature == "ClusterProminence")			{ outputIt.Set( texturesCalculator->GetClusterProminence() ); }
			else if (m_Feature == "HaralickCorrelation")		{ outputIt.Set( texturesCalculator->GetHaralickCorrelation() ); }

			/*energyIt.Set(texturesCalculator->GetEnergy());
		    entropyIt.Set(texturesCalculator->GetEntropy());
			correlationIt.Set(texturesCalculator->GetCorrelation());
			invDiffMomentIt.Set(texturesCalculator->GetInverseDifferenceMoment());
			inertiaIt.Set(texturesCalculator->GetInertia());
			clusterShadeIt.Set(texturesCalculator->GetClusterShade());
			clusterProminenceIt.Set(texturesCalculator->GetClusterProminence());
			haralickCorIt.Set(texturesCalculator->GetHaralickCorrelation());*/
		}

		// Update progress
		progress.CompletedPixel();

		// Increment iterators
		++maskIt;
		++outputIt;
		/*++energyIt;
	    ++entropyIt;
		++correlationIt;
		++invDiffMomentIt;
		++inertiaIt;
		++clusterShadeIt;
		++clusterProminenceIt;
		++haralickCorIt;*/
    }

}

} // End namespace otb

#endif
