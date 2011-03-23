/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHessianSmoothed3DToVesselnessMeasureImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 20:59:44 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHessianSmoothed3DToVesselnessMeasureImageFilter_txx
#define __itkHessianSmoothed3DToVesselnessMeasureImageFilter_txx

#include "itkHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vnl/vnl_math.h"

#include "itkCastImageFilter.h"

#define EPSILON  1e-03

namespace itk
{

/**
 * Constructor
 */
template < typename TPixel >
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::HessianSmoothed3DToVesselnessMeasureImageFilter()
{
  m_Alpha = 0.5;
  m_Beta  = 0.3;
  m_Gamma = 0.3;
  m_Eta   = 0.2;

  m_C = 10e-6;

  m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  m_SymmetricEigenValueFilter->SetDimension( ImageDimension );
  m_SymmetricEigenValueFilter->OrderEigenValuesBy( EigenAnalysisFilterType::FunctorType::OrderByValue );
  
  // By default, scale the vesselness measure by the largest eigenvalue
  m_ScaleVesselnessMeasure  = false;
  
}


template < typename TPixel >
void 
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::GenerateData()
{
  itkDebugMacro(
      << "HessianSmoothed3DToVesselnessMeasureImageFilter generating data ");

  m_SymmetricEigenValueFilter->SetInput( this->GetInput() );
  
  typename OutputImageType::Pointer output = this->GetOutput();

  m_Sigma = this->GetSigma();

  typedef typename EigenAnalysisFilterType::OutputImageType
                                            EigenValueImageType;

  m_SymmetricEigenValueFilter->Update();
  
  const typename EigenValueImageType::ConstPointer eigenImage = 
                    m_SymmetricEigenValueFilter->GetOutput();
  
  // walk the region of eigen values and get the vesselness measure
  EigenValueArrayType eigenValue;
  ImageRegionConstIterator<EigenValueImageType> it;
  it = ImageRegionConstIterator<EigenValueImageType>(
      eigenImage, eigenImage->GetRequestedRegion());
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());

  // Create input iterator
  ImageRegionConstIterator<InputImageType> input_it;
  input_it = ImageRegionConstIterator<InputImageType>(
      this->GetInput(), this->GetInput()->GetRequestedRegion());

  oit.GoToBegin();
  it.GoToBegin();
  input_it.GoToBegin();

  typedef Image<float, 3> FloatImageType;

  std::vector< FloatImageType::Pointer > LambdaImageVector(3);
  std::vector< FloatImageType::Pointer > FImageVector(6); // A, B, C1, C2, Cup, Rut

  std::vector< ImageRegionIterator<FloatImageType> > LambdaIterVector(3);
  std::vector< ImageRegionIterator<FloatImageType> > FIterVector(6);

  for (int i=0; i<3; i++)
  {
	  LambdaImageVector[i] = FloatImageType::New();
	  LambdaImageVector[i]->SetRegions(output->GetRequestedRegion());
	  LambdaImageVector[i]->Allocate();

	  LambdaIterVector[i] = ImageRegionIterator<FloatImageType>( LambdaImageVector[i], output->GetRequestedRegion() );
	  LambdaIterVector[i].GoToBegin();
  }

  for (int i=0; i<6; i++)
  {
	  FImageVector[i] = FloatImageType::New();
	  FImageVector[i]->SetRegions(output->GetRequestedRegion());
	  FImageVector[i]->Allocate();

	  FIterVector[i] = ImageRegionIterator<FloatImageType>( FImageVector[i], output->GetRequestedRegion() );
	  FIterVector[i].GoToBegin();
  }

  // write eigen values to text file
  //std::ofstream file;
  //file.open("ev.txt");

  while (!it.IsAtEnd())
    {
    // Get the eigen value
    eigenValue = it.Get();


    // Find the smallest eigenvalue
    double smallest = vnl_math_abs( eigenValue[0] );
    double Lambda1 = eigenValue[0];
 
    for ( unsigned int i=1; i <=2; i++ )
      {
      if ( vnl_math_abs( eigenValue[i] ) < smallest )
        {
        Lambda1 = eigenValue[i];
        smallest = vnl_math_abs( eigenValue[i] );
        }
      }

    // Find the largest eigenvalue
    double largest = vnl_math_abs( eigenValue[0] );
    double Lambda3 = eigenValue[0];
 
    for ( unsigned int i=1; i <=2; i++ )
      {
      if (  vnl_math_abs( eigenValue[i] ) > largest  )
        {
        Lambda3 = eigenValue[i];
        largest = vnl_math_abs( eigenValue[i] );
        }
      }


    //  find Lambda2 so that |Lambda1| < |Lambda2| < |Lambda3|
    double Lambda2 = eigenValue[0];

    for ( unsigned int i=0; i <=2; i++ )
      {
      if ( eigenValue[i] != Lambda1 && eigenValue[i] != Lambda3 )
        {
        Lambda2 = eigenValue[i];
        break;
        }
      }

	LambdaIterVector[0].Set(Lambda1);
	LambdaIterVector[1].Set(Lambda2);
	LambdaIterVector[2].Set(Lambda3);

	//file << Lambda1 << " " << Lambda2 << " " << Lambda3 << "\n";

    //if ( Lambda3 <= 100.0 || 
     //    vnl_math_abs( Lambda3 ) < EPSILON )


	if ( ( vnl_math_abs( Lambda3 ) >= 100 && input_it.Get() >= 200 ) ||
		 ( vnl_math_abs( Lambda3 ) <= 100 && input_it.Get() >= -250 && input_it.Get() <= 200 )
		) 
      {
      
   
      double Lambda1Abs = vnl_math_abs( Lambda1 );
      double Lambda2Abs = vnl_math_abs( Lambda2 );
      double Lambda3Abs = vnl_math_abs( Lambda3 );

	  double Ra = Lambda1Abs / vcl_sqrt( Lambda2Abs * Lambda3Abs );

	  double Fa = vcl_exp(  - vnl_math_sqr( Ra ) / (2 * vnl_math_sqr( m_Alpha ) ) );
		
	  double Rb = Lambda2Abs / Lambda3Abs;

	  double Fb = vcl_exp( - vnl_math_sqr( Rb - m_Gamma ) / (2 * vnl_math_sqr( m_Beta ) ) );

	  double Fc1 = (1.0 - vcl_exp( - vnl_math_sqr( Lambda1Abs / Lambda2Abs ) / (2 * vnl_math_sqr( m_Eta ) ) ) );

	  double Fc2 = (1.0 - vcl_exp( - vnl_math_sqr( Lambda2Abs / Lambda3Abs ) / (2 * vnl_math_sqr( m_Eta ) ) ) );

	   if(  m_ScaleVesselnessMeasure ) 
       {
			Fa *= Lambda3Abs;
			Fb *= Lambda3Abs;
			Fc1 *= Lambda3Abs;
			Fc2 *= Lambda3Abs;
       }

	  double Frut = Fa * Fb;

	  double Fcup =  Fc1 * Fc2;

	  double H = vnl_math_max( Frut, Fcup );

	  // A, B, C1, C2, Cup, Rut
	  FIterVector[0].Set(Fa);
	  FIterVector[1].Set(Fb);
	  FIterVector[2].Set(Fc1);
	  FIterVector[3].Set(Fc2);
	  FIterVector[4].Set(Fcup);
	  FIterVector[5].Set(Frut);

      if(  m_ScaleVesselnessMeasure ) 
        {
        oit.Set( static_cast< OutputPixelType >(
                                     Lambda3Abs*H ) );
        }
      else
        {
        oit.Set( static_cast< OutputPixelType >( H ) );
        }
      } else {
		oit.Set( NumericTraits< OutputPixelType >::Zero );
	  }

    ++it;
    ++oit;
	++input_it;

	for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
	for (int i=0; i<6; i++) { ++FIterVector[i]; }

    }

    typedef itk::ImageFileWriter< FloatImageType >  WriterType;

    for (int i=0; i<3; i++)
	{
		WriterType::Pointer writer = WriterType::New();
		WriterType::SetGlobalWarningDisplay(false);
	   
		std::stringstream ss;
		ss << "Lambda_sigma_" << m_Sigma << "_" << i << ".nii";
		writer->SetFileName(ss.str().c_str());
		writer->SetInput( LambdaImageVector[i] );
		std::cout<<"Writing: "<<ss.str()<<std::endl;
		writer->Update();
		writer.~SmartPointer();
	}

	std::vector<std::string> fs(6);
	fs[0]="a";
	fs[1]="b";
	fs[2]="c1";
	fs[3]="c2";
	fs[4]="cup";
	fs[5]="rut";

	for (int i=0; i<6; i++)
	{
		WriterType::Pointer writer = WriterType::New();
		WriterType::SetGlobalWarningDisplay(false);
	   
		std::stringstream ss;
		ss << "F_" << fs[i] << "_" << m_Sigma << ".nii";
		writer->SetFileName(ss.str().c_str());
		writer->SetInput( FImageVector[i] );
		std::cout<<"Writing: "<<ss.str()<<std::endl;
		writer->Update();
		writer.~SmartPointer();
	}

	typedef itk::CastImageFilter< OutputImageType, FloatImageType >  CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(output);
	castFilter->Update();

	WriterType::Pointer writer = WriterType::New();
	
	std::stringstream ss;
	ss << "H_sigma_" << m_Sigma << ".nii";
	writer->SetFileName(ss.str().c_str());
	writer->SetInput( castFilter->GetOutput() );
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();

 //file.close();
    
}

template < typename TPixel >
void
HessianSmoothed3DToVesselnessMeasureImageFilter< TPixel >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta:  " << m_Beta  << std::endl;
  os << indent << "Gamma: " << m_Gamma << std::endl;
  os << indent << "Eta: "   << m_Eta   << std::endl;

  os << indent << "C: " << m_C << std::endl;
}


} // end namespace itk
  
#endif
