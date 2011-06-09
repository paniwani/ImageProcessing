/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator_h
#define __itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator_h

#include "itkHistogram.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkSimpleDataObjectDecorator.h"

namespace itk {
namespace Statistics {

/** \class GreyLevelRunLengthMatrixTextureCoefficientsCalculator
*  \brief This class computes texture feature coefficients from a grey level
* run-length matrix.
*/

template< class THistogram >
class ITK_EXPORT GreyLevelRunLengthMatrixTextureCoefficientsCalculator : public ProcessObject
{
public:
  /** Standard typedefs */
  typedef GreyLevelRunLengthMatrixTextureCoefficientsCalculator     Self;
  typedef ProcessObject                          Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( GreyLevelRunLengthMatrixTextureCoefficientsCalculator, ProcessObject );

  /** standard New() method support */
  itkNewMacro( Self );

  typedef THistogram                                      HistogramType;
  typedef typename HistogramType::Pointer                 HistogramPointer;
  typedef typename HistogramType::ConstPointer            HistogramConstPointer;
  typedef typename HistogramType::MeasurementType         MeasurementType;
  typedef typename HistogramType::MeasurementVectorType   MeasurementVectorType;
  typedef typename HistogramType::IndexType               IndexType;
  typedef typename HistogramType::
    TotalAbsoluteFrequencyType                            FrequencyType;

  /** Triggers the Computation of the histogram */
  void Compute( void );

  /** Method to Set/Get the input Histogram */
  void SetInput ( const HistogramType * histogram );
  const HistogramType * GetInput() const;

  /** Smart Pointer type to a DataObject. */
  typedef DataObject::Pointer                   DataObjectPointer;

  /** Type of DataObjects used for scalar outputs */
  typedef SimpleDataObjectDecorator<MeasurementType>     MeasurementObjectType;

  /** Methods to return the short run emphasis. */
  MeasurementType GetShortRunEmphasis() const;
  const MeasurementObjectType* GetShortRunEmphasisOutput() const;

  /** Methods to return the long run emphasis. */
  MeasurementType GetLongRunEmphasis() const;
  const MeasurementObjectType* GetLongRunEmphasisOutput() const;

  /** Methods to return the grey level nonuniformity. */
  MeasurementType GetGreyLevelNonuniformity() const;
  const MeasurementObjectType* GetGreyLevelNonuniformityOutput() const;

  /** Methods to return the run length nonuniformity. */
  MeasurementType GetRunLengthNonuniformity() const;
  const MeasurementObjectType* GetRunLengthNonuniformityOutput() const;

  /** Methods to return the low grey level run emphasis. */
  MeasurementType GetLowGreyLevelRunEmphasis() const;
  const MeasurementObjectType* GetLowGreyLevelRunEmphasisOutput() const;

  /** Methods to return the high grey level run emphasis. */
  MeasurementType GetHighGreyLevelRunEmphasis() const;
  const MeasurementObjectType* GetHighGreyLevelRunEmphasisOutput() const;

  /** Methods to return the short run low grey level run emphasis. */
  MeasurementType GetShortRunLowGreyLevelEmphasis() const;
  const MeasurementObjectType* GetShortRunLowGreyLevelEmphasisOutput() const;

  /** Methods to return the short run high grey level run emphasis. */
  MeasurementType GetShortRunHighGreyLevelEmphasis() const;
  const MeasurementObjectType* GetShortRunHighGreyLevelEmphasisOutput() const;

  /** Methods to return the long run low grey level run emphasis. */
  MeasurementType GetLongRunLowGreyLevelEmphasis() const;
  const MeasurementObjectType* GetLongRunLowGreyLevelEmphasisOutput() const;

  /** Methods to return the long run high grey level run emphasis. */
  MeasurementType GetLongRunHighGreyLevelEmphasis() const;
  const MeasurementObjectType* GetLongRunHighGreyLevelEmphasisOutput() const;

  itkGetMacro( TotalNumberOfRuns, unsigned long );

  typedef enum
    {
    ShortRunEmphasis,
    LongRunEmphasis,
    GreyLevelNonuniformity,
    RunLengthNonuniformity,
    LowGreyLevelRunEmphasis,
    HighGreyLevelRunEmphasis,
    ShortRunLowGreyLevelEmphasis,
    ShortRunHighGreyLevelEmphasis,
    LongRunLowGreyLevelEmphasis,
    LongRunHighGreyLevelEmphasis
    }  RunLengthFeatureName;

  /** convenience method to access the run length values */
  MeasurementType GetFeature( RunLengthFeatureName name );

protected:
  GreyLevelRunLengthMatrixTextureCoefficientsCalculator();
  ~GreyLevelRunLengthMatrixTextureCoefficientsCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Make a DataObject to be used for output output. */
  virtual DataObjectPointer MakeOutput( unsigned int );


  void GenerateData();


private:
  GreyLevelRunLengthMatrixTextureCoefficientsCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned long                           m_TotalNumberOfRuns;

};

} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator.txx"
#endif

#endif
