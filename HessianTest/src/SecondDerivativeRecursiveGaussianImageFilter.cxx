/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SecondDerivativeRecursiveGaussianImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2005-08-27 01:46:02 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

/*-----------------------------------------------
STRATEGY

1. Threshold at 200 HU 
2. Compute second derivatives using Gaussian convolution
	- Note: Set sigma to 1 voxel
3. Construct Hessian matrix
4. Solve for eigenvalues and eigenvectors
5. Use eigenvalues in F_cup and F_rut enhancement functions
6. Plot

-----------------------------------------------*/

#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkImage.h"
#include <string>
#include <itkImageRegionIteratorWithIndex.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix.h>
#include <stdlib.h>
#include <vector>

double functionA(double lambda[], double alpha);
double functionB(double lambda[], double beta, double gamma);
double functionC(double ev1, double ev2, double eta);
double f_rut(double lambda[], double alpha, double beta, double gamma);
double f_cup(double lambda[], double eta);
//int compare (const void * a, const void * b);
//int compare (const double a, const double b);
bool compare (double a, double b);

int main(int argc, char * argv [] )
{
	/*
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "SecondDerivativeRecursiveGaussianImageFilter input colon_mask colon_mask_shaded " << std::endl;
		system("PAUSE");
		return EXIT_FAILURE;
	}
	*/

	typedef float            PixelType;
	typedef float            OutputPixelType;

	const unsigned int  Dimension = 3;

	typedef itk::Image< PixelType,       Dimension >  ImageType;
	typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;
	typedef itk::ImageRegionIteratorWithIndex<ImageType>	IteratorType;

	typedef itk::ImageFileReader< ImageType       >   ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >   WriterType;

	typedef itk::ImageDuplicator< OutputImageType >   DuplicatorType;

	typedef itk::RecursiveGaussianImageFilter< 
								  ImageType, 
								  ImageType >  FilterType;

	ReaderType::Pointer  reader  = ReaderType::New();
	WriterType::Pointer  writer  = WriterType::New();

	DuplicatorType::Pointer duplicator  = DuplicatorType::New();

	//Read input file
	std::cout << "Reading input file." << std::endl;
	reader->SetFileName( argv[1] );

	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem reading the input file" << std::endl;
		std::cerr << excp << std::endl;
		system("PAUSE");
		return EXIT_FAILURE;
	}

	ImageType::Pointer input = reader->GetOutput();

	/*
	//Read colon mask
	std::cout << "Reading colon mask..." << std::endl;
	reader->SetFileName( argv[2] );

	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem reading the colon mask" << std::endl;
		std::cerr << excp << std::endl;
		system("PAUSE");
		return EXIT_FAILURE;
	}

	ImageType::Pointer colon_mask = reader->GetOutput();
	
	//Read colon mask shaded
	std::cout << "Reading colon mask shaded..." << std::endl;
	reader->SetFileName( argv[3] );

	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem reading the colon mask shaded" << std::endl;
		std::cerr << excp << std::endl;
		system("PAUSE");
		return EXIT_FAILURE;
	}

	ImageType::Pointer colon_mask_shaded = reader->GetOutput();
	*/

	//Setup second derivative computation

	std::cout << "Computing second derivatives." << std::endl;

	FilterType::Pointer ga = FilterType::New();
	FilterType::Pointer gb = FilterType::New();
	FilterType::Pointer gc = FilterType::New();

	ga->SetDirection( 0 );
	gb->SetDirection( 1 );
	gc->SetDirection( 2 );

	ga->SetZeroOrder();
	gb->SetZeroOrder();
	gc->SetSecondOrder();

	//Set all sigma values
	const int numOfSigma = 1;
	float sigma_ar[numOfSigma][3] = {  {.8203, .8203, 1} };

	//Set params
	const int numOfAlpha = 1;
	double alpha[numOfAlpha] = {.5};
	
	double beta = 0.3;
	double gamma = 0.3;
	double eta = 0.2;

	/*float sigma_ar[13][3] = {	{.8203, .8203, 2.5 },	//anisotropic gaussian
								{.0104, .0104, 1},
								{1.6406, 1.6406, 5},
								{.0208, .0208, 2},
								{1.6406, 1.6406, 2.5},
								{.0208, .0208, 1},
	{5,5,5},	//isotropic
	{1,1,1},
	{2,2,2},
	{2.5,2.5,2.5},
	{1.6406,1.6406,1.6406},
	{.8203,.8203,.8203},
	{.0104,.0104,.0104}}; */

	//float sigma_ar[] = {.0104, .0208, 1, 2, .8203, 1.6406, 2.5, 5};	//isotropic gaussian
	

	std::cout << "Allocating images." << std::endl;

	ImageType::Pointer F_rut_image = ImageType::New();
	F_rut_image->SetRegions(input->GetRequestedRegion());
	F_rut_image->SetSpacing(input->GetSpacing());
	F_rut_image->Allocate();

	ImageType::Pointer F_a_image = ImageType::New();
	F_a_image->SetRegions(input->GetRequestedRegion());
	F_a_image->SetSpacing(input->GetSpacing());
	F_a_image->Allocate();

	ImageType::Pointer F_b_image = ImageType::New();
	F_b_image->SetRegions(input->GetRequestedRegion());
	F_b_image->SetSpacing(input->GetSpacing());
	F_b_image->Allocate();


	/*
	ImageType::Pointer F_cup_image = ImageType::New();
	F_cup_image->SetRegions(input->GetRequestedRegion());
	F_cup_image->SetSpacing(input->GetSpacing());
	F_cup_image->Allocate();
	*/

	IteratorType input_iter(input, input->GetRequestedRegion() );
	//IteratorType colon_mask_iter(colon_mask, colon_mask->GetRequestedRegion() );
	//IteratorType colon_mask_shaded_iter(colon_mask_shaded, colon_mask_shaded->GetRequestedRegion() );
	IteratorType F_rut_iter(F_rut_image, F_rut_image->GetRequestedRegion() );
	//IteratorType F_cup_iter(F_cup_image, F_cup_image->GetRequestedRegion() );
	IteratorType F_a_iter(F_a_image, F_a_image->GetRequestedRegion() );
	IteratorType F_b_iter(F_b_image, F_b_image->GetRequestedRegion() );
	
	for (int n=0; n < numOfAlpha ; n++) {
	for (int k=0; k < numOfSigma ; k++) {
		ga->SetSigma( sigma_ar[k][0] );
		gb->SetSigma( sigma_ar[k][1] );
		gc->SetSigma( sigma_ar[k][2] );
		
		ga->SetInput( input );
		gb->SetInput( ga->GetOutput() );
		gc->SetInput( gb->GetOutput() );

		duplicator->SetInputImage( gc->GetOutput() );

		gc->Update(); 
		duplicator->Update();

		ImageType::Pointer Izz = duplicator->GetOutput();

		/*
		writer->SetInput( Izz );
		outputFileName = outputPrefix + "-Izz.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/

		gc->SetDirection( 1 );  // gc now works along Y
		gb->SetDirection( 2 );  // gb now works along Z

		gc->Update();
		duplicator->Update();

		ImageType::Pointer Iyy = duplicator->GetOutput();

		/*
		writer->SetInput( Iyy );
		outputFileName = outputPrefix + "-Iyy.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/

		gc->SetDirection( 0 );  // gc now works along X
		ga->SetDirection( 1 );  // ga now works along Y

		gc->Update();
		duplicator->Update();

		ImageType::Pointer Ixx = duplicator->GetOutput();

		/*
		writer->SetInput( Ixx );
		outputFileName = outputPrefix + "-Ixx.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/


		ga->SetDirection( 0 );
		gb->SetDirection( 1 );
		gc->SetDirection( 2 );

		ga->SetZeroOrder();
		gb->SetFirstOrder();
		gc->SetFirstOrder();

		gc->Update();
		duplicator->Update();

		ImageType::Pointer Iyz = duplicator->GetOutput();

		/*
		writer->SetInput( Iyz );
		outputFileName = outputPrefix + "-Iyz.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/

		ga->SetDirection( 1 );
		gb->SetDirection( 0 );
		gc->SetDirection( 2 );

		ga->SetZeroOrder();
		gb->SetFirstOrder();
		gc->SetFirstOrder();

		gc->Update();
		duplicator->Update();

		ImageType::Pointer Ixz = duplicator->GetOutput();

		/*
		writer->SetInput( Ixz );
		outputFileName = outputPrefix + "-Ixz.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/

		ga->SetDirection( 2 );
		gb->SetDirection( 0 );
		gc->SetDirection( 1 );

		ga->SetZeroOrder();
		gb->SetFirstOrder();
		gc->SetFirstOrder();

		gc->Update();
		duplicator->Update();

		ImageType::Pointer Ixy = duplicator->GetOutput();

		/*
		writer->SetInput( Ixy );
		outputFileName = outputPrefix + "-Ixy.hdr";
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();
		*/

		//Iterate through images
		IteratorType Ixx_iter(Ixx, Ixx->GetRequestedRegion() );
		IteratorType Iyy_iter(Iyy, Iyy->GetRequestedRegion() );
		IteratorType Izz_iter(Izz, Izz->GetRequestedRegion() );
		IteratorType Ixy_iter(Ixy, Ixy->GetRequestedRegion() );
		IteratorType Ixz_iter(Ixz, Ixz->GetRequestedRegion() );
		IteratorType Iyz_iter(Iyz, Iyz->GetRequestedRegion() );
		
		for(input_iter.GoToBegin(), Ixx_iter.GoToBegin(), Iyy_iter.GoToBegin(), Izz_iter.GoToBegin(), Ixy_iter.GoToBegin(), Ixz_iter.GoToBegin(), Iyz_iter.GoToBegin(), F_rut_iter.GoToBegin(), F_a_iter.GoToBegin(), F_b_iter.GoToBegin();//, F_cup_iter.GoToBegin();//, colon_mask_iter.GoToBegin();//, colon_mask_shaded_iter.GoToBegin(); 
			!input_iter.IsAtEnd() && !Ixx_iter.IsAtEnd() && !Iyy_iter.IsAtEnd() && !Izz_iter.IsAtEnd() && !Ixy_iter.IsAtEnd() && !Ixz_iter.IsAtEnd() && !Iyz_iter.IsAtEnd() && !F_rut_iter.IsAtEnd() && !F_a_iter.IsAtEnd() && !F_b_iter.IsAtEnd();//!F_cup_iter.IsAtEnd();// && !colon_mask_iter.IsAtEnd();// && !colon_mask_shaded_iter.IsAtEnd();
			++input_iter, ++Ixx_iter, ++Iyy_iter, ++Izz_iter, ++Ixy_iter, ++Ixz_iter, ++Iyz_iter, ++F_rut_iter, ++F_a_iter, ++F_b_iter) {// ++F_cup_iter){//, ++colon_mask_iter) {// ++colon_mask_shaded_iter) {
	
			//if (input_iter.Get() < 200 && input_iter.Get() > -500) {	//threshold
			//if (colon_mask_shaded_iter.Get() == 128) {	//look only in shaded region
			//if (colon_mask_iter.Get() == 40) {	//look only inside the colon

			if (input_iter.Get() == 40) { //input is the colon mask, look only inside
			
				//Construct Hessian matrix using symmetry
				vnl_matrix_fixed<float,3,3> H;
				H[0][0] = Ixx_iter.Get();	H[0][1] = Ixy_iter.Get();	H[0][2] = Ixz_iter.Get();
				H[1][0] = Ixy_iter.Get();	H[1][1] = Iyy_iter.Get();	H[1][2] = Iyz_iter.Get();
				H[2][0] = Ixz_iter.Get();	H[2][1] = Iyz_iter.Get();	H[2][2] = Izz_iter.Get();

				//Solve eigensystem
				vnl_symmetric_eigensystem<float> es(H);
				double lambda[3] = {es.get_eigenvalue(2), es.get_eigenvalue(1), es.get_eigenvalue(0)};
				
				//Sort eigenvalues in ascending absolute value
				std::vector<double> lambda_vector (lambda, lambda+3);
				std::vector<double>::iterator it;
				sort(lambda_vector.begin(), lambda_vector.end(), compare);
				
				int i=0;
				for (it=lambda_vector.begin(); it!=lambda_vector.end(); ++it) {
						lambda[i++] = *it;
				}

				if (lambda[2] > 0) {

					//std::cout << "lambda: " << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << std::endl;
					
					//Solve functions
					double fa = functionA(lambda, alpha[n]);
					//double fb = functionB(lambda, beta, gamma);

					F_a_iter.Set(fa);
					//F_b_iter.Set(fb);
					//F_rut_iter.Set(fa*fb);

					//F_rut_iter.Set(f_rut(lambda, alpha, beta, gamma));
					//F_cup_iter.Set(f_cup(lambda, eta));
				} else {
					F_a_iter.Set(0);
				}
				
			} else {
				//F_rut_iter.Set(0);
				//F_cup_iter.Set(0);
				F_a_iter.Set(0);
				//F_b_iter.Set(0);
			}
		}

		//Write images
		
		
		std::cout << "Writing images." << std::endl;
		
		/*
		
		writer->SetInput( F_rut_image );
		std::stringstream ss;
		ss << "F_rut" << "_sigmaX" << sigma_ar[k][0] << "_sigmaY" << sigma_ar[k][1] << "_sigmaZ" << sigma_ar[k][2] << "_alpha " << alpha[k]".hdr";
		writer->SetFileName(ss.str().c_str());
		writer->Update();

		*/

		writer->SetInput( F_a_image );
		std::stringstream ss;
		ss << "F_a" << "_sigmaX" << sigma_ar[k][0] << "_sigmaY" << sigma_ar[k][1] << "_sigmaZ" << sigma_ar[k][2] << "_alpha " << alpha[n] << ".hdr";
		writer->SetFileName(ss.str().c_str());//
		writer->Update();


		
		/*
		std::cout << "Writing F_cup" << std::endl;

		ss.str("");
		ss << "F_cup" << "_sigmaX" << sigma_ar[k][0] << "_sigmaY" << sigma_ar[k][1] << "_sigmaZ" << sigma_ar[k][2] << ".hdr";

		writer->SetInput( F_cup_image );
		writer->SetFileName(ss.str().c_str());
		writer->Update();
		*/

	}
	}

	system("PAUSE");
	return EXIT_SUCCESS;
}

double functionA(double lambda[], double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs(lambda[1]*lambda[2]));
	//std::cout << "Ra: " << Ra << std::endl;
	return exp(-vnl_math_sqr(Ra)/(2*vnl_math_sqr(alpha)));
}

double functionB(double lambda[], double beta, double gamma) {
	double Rb;
	Rb = abs(lambda[1]/lambda[2]);
	//std::cout << "Rb: " << Rb << std::endl;
	return exp(-vnl_math_sqr(Rb - gamma)/(2*vnl_math_sqr(beta)));
}

double functionC(double ev1, double ev2, double eta) {
	double Rc;
	Rc = abs(ev1)/abs(ev2);
	return 1.0 - exp(-vnl_math_sqr(Rc)/(2*vnl_math_sqr(eta)));
}

double f_rut(double lambda[], double alpha, double beta, double gamma) {
	double fa = functionA(lambda, alpha);
	double fb = functionB(lambda, beta, gamma);
	//std::cout << "fa: " << fa << std::endl;
	//std::cout << "fb: " << fb << std::endl;
	
	double frut = fa*fb;
	//std::cout << "frut " << frut << std::endl;

	return frut;
	//return functionA(lambda, alpha)*functionB(lambda, beta, gamma);
}

double f_cup(double lambda[], double eta) {
	return functionC(lambda[0], lambda[1], eta)*functionC(lambda[1], lambda[2], eta);
}

//int compare (const void * a, const void * b) {
//	return ( abs(*(double*)a) - abs(*(double*)b) );
//}

bool compare (double a, double b) {
	return ( abs(a) < abs(b) );
}