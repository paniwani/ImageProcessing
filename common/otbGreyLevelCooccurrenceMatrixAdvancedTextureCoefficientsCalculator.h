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
#ifndef __otbGreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator_h
#define __otbGreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator_h

#include "itkHistogram.h"
#include "itkMacro.h"

namespace otb {

/** \class GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator
 *  \brief This class computes texture feature coefficients from a grey level
 * co-occurrence matrix.
 *
 * This class computes features that summarize image texture, given a grey level
 * co-occurrence matrix (generated by a ScalarImageToGreyLevelCooccurrenceMatrixGenerator
 * or related class).
 *
 * The features calculated are as follows (where \f$ g(i, j) \f$ is the element in
 * cell i, j of a a normalized GLCM):
 *
 * "Mean" \f$ = \sum_{i,j}i g(i,j) \f$
 *
 * "Sum of squares: Variance" \f$ = f_4 = \sum_{i,j}(i -mu)^2 g(i,j) \f$
 *
 * "Sum average" \f$ = f_6 = -\sum_{i}i g_{x+y}(i)
 *
 * "Sum Variance" \f$ = f_7 = \sum_{i}(i - f_8)^2 g_{x+y}(i) \f$
 *
 * "Sum Entropy" \f$= f_8 = -\sum_{i}g_{x+y}(i) log (g_{x+y}(i)) \f$
 *
 * "Difference variance" \f$ = f_10 = variance of g_{x-y}(i)
 *
 * "Difference entropy" \f$ = f_11 = -\sum_{i}g_{x-y}(i) log (g_{x-y}(i)) \f$
 *
 * "Information Measures of Correlation IC1" \f$ = f_12 = \frac{f_9 - HXY1}{H} \f$
 *
 * "Information Measures of Correlation IC2" \f$ = f_13 = \sqrt{1 - \exp{-2}|HXY2 - f_9|} \f$
 *
 * Above, \f$ \mu =  \f$ (weighted pixel average) \f$ = \sum_{i,j}i \cdot g(i, j) =
 * \sum_{i,j}j \cdot g(i, j) \f$ (due to matrix summetry), and
 *
 * \f$ \g_{x+y}(k) =  \sum_{i}\sum_{j}g(i)\f$ where \f$ i+j=k \f$ and \f$ k = 2,3,..,2N_[g}  \f$ and
 *
 * \f$ \g_{x-y}(k) =  \sum_{i}\sum_{j}g(i)\f$ where \f$ i-j=k \f$ and \f$ k = 0,1,..,N_[g}-1  \f$
 *
 * NOTA BENE: The input histogram will be forcably normalized!
 * This algorithm takes three passes through the input
 * histogram if the histogram was already normalized, and four if not.
 *
 * Web references:
 *
 * http://www.cssip.uq.edu.au/meastex/www/algs/algs/algs.html
 * http://www.ucalgary.ca/~mhallbey/texture/texture_tutorial.html
 *
 * Print references:
 *
 * Haralick, R.M., K. Shanmugam and I. Dinstein. 1973.  Textural Features for
 * Image Classification. IEEE Transactions on Systems, Man and Cybernetics.
 * SMC-3(6):610-620.
 *
 * Haralick, R.M. 1979. Statistical and Structural Approaches to Texture.
 * Proceedings of the IEEE, 67:786-804.
 *
 * R.W. Conners and C.A. Harlow. A Theoretical Comaprison of Texture Algorithms.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence,  2:204-222, 1980.
 *
 * R.W. Conners, M.M. Trivedi, and C.A. Harlow. Segmentation of a High-Resolution
 * Urban Scene using Texture  Operators. Computer Vision, Graphics and Image
 * Processing, 25:273-310,  1984.
 *
 * \sa ScalarImageToGreyLevelCooccurrenceMatrixGenerator
 * \sa GreyLevelCooccurrenceMatrixTextureCoefficientsCalculator
 * \sa ScalarImageTextureCalculator
 *
 */

template <typename THistogram>
class ITK_EXPORT GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator :
  public itk::Object
{
public:
  /** Standard typedefs */
  typedef GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator Self;
  typedef itk::Object                                                      Superclass;
  typedef itk::SmartPointer<Self>                                          Pointer;
  typedef itk::SmartPointer<const Self>                                    ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator, itk::Object);

  /** standard New() method support */
  itkNewMacro(Self);

  typedef THistogram                                    HistogramType;
  typedef typename HistogramType::Pointer               HistogramPointer;
  typedef typename HistogramType::ConstPointer          HistogramConstPointer;
  typedef typename HistogramType::MeasurementType       MeasurementType;
  typedef typename HistogramType::RelativeFrequencyType RelativeFrequencyType;
  typedef typename HistogramType::MeasurementVectorType MeasurementVectorType;
  typedef typename HistogramType::IndexType             IndexType;
  typedef typename HistogramType::AbsoluteFrequencyType         FrequencyType; //FIXME several type in the new stat framework

  /** Connects the GLCM histogram over which the features are going to be computed */
  itkSetObjectMacro(Histogram, HistogramType);
  itkGetObjectMacro(Histogram, HistogramType);

  itkGetMacro(Mean, double);
  itkGetMacro(Variance, double);
  itkGetMacro(SumAverage, double);
  itkGetMacro(SumVariance, double);
  itkGetMacro(SumEntropy, double);
  itkGetMacro(DifferenceEntropy, double);
  itkGetMacro(DifferenceVariance, double);
  itkGetMacro(IC1, double);
  itkGetMacro(IC2, double);

  /** Triggers the Computation of the histogram */
  void Compute(void);

protected:
  GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator() {};
  virtual ~GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  GreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

  HistogramPointer m_Histogram;

  void ComputeMean();
  double m_Mean, m_Variance, m_SumAverage, m_SumVariance, m_SumEntropy, m_DifferenceEntropy,
         m_DifferenceVariance, m_IC1, m_IC2;

  double ComputePS(long unsigned int k);
  double ComputePD(long unsigned int k);
};
} // end of namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbGreyLevelCooccurrenceMatrixAdvancedTextureCoefficientsCalculator.txx"
#endif

#endif
