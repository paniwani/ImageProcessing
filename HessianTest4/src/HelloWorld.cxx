#include <stdlib.h>
#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include <vector>

typedef								  float				O_PixelType;				//output type
typedef itk::Image<  O_PixelType,  3 >  O_ImageType;

double functionA(std::vector<double>& lambda, double alpha);
double functionB(std::vector<double>& lambda, double beta, double gamma);
double functionC(double ev1, double ev2, double eta);
double f_rut(std::vector<double>& lambda, double alpha, double beta, double gamma);
double f_cup(std::vector<double>& lambda, double eta);
bool compare (double a, double b);
void WriteITK(O_ImageType::Pointer image, const char * name);

int main( int argc, char *argv[] )
{
	// INPUTS
	char *  in_file_name      =					argv[1];
	char *  in_colon_mask	  =					argv[2];

	// SET PARAMETERS
	/*
	double sigma[3] = {1,2,3};
	double alpha[1] = {.5};
	double beta[3] = {.1,.3,.5};
	double gamma = 0.3;
	double eta[1] = {0.2};
	*/

	/*
	double sigma[3] = {1,2,3};//0.8203,1.6406};
	double alpha[4] = {.1,.25,.5,1};
	double beta[3] = {.1,.3,.5};
	double gamma = 0.3;
	double eta[3] = {0.1,0.2,0.5};
	*/

	double sigma[2] = {2.0,3.0};//0.8203,1.6406};
	double alpha[1] = {0.3};
	double beta[1] = {0.3};
	double gamma = {0.3};
	double eta[1] = {0.2} ;

	// Flag to use input or mask in hessian
	bool useMask = false;

	std::stringstream outputPrefix;
	outputPrefix << "output_" << "realInput_allColon_" << in_file_name;

	// ITK Pixel Types
	typedef                                 unsigned char             F_PixelType;
	typedef itk::SymmetricSecondRankTensor< float,       3 >  T_PixelType;
	typedef itk::Vector<                    float,       3 >  E_PixelType;
	typedef itk::Vector<                    E_PixelType, 3 > EV_PixelType;
	
	typedef								    float				M_PixelType;				//mask type


	// ITK Image Types
	typedef itk::Image<  F_PixelType,  3 >  F_ImageType;
	typedef itk::Image<  T_PixelType,  3 >  T_ImageType;
	typedef itk::Image<  M_PixelType,  3 >  M_ImageType;

	// Read Input
	typedef itk::ImageFileReader< F_ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( in_file_name );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	F_ImageType::Pointer input = reader->GetOutput();

	
	// Read colon mask
	typedef itk::ImageFileReader< M_ImageType > M_ReaderType;
	M_ReaderType::Pointer M_reader = M_ReaderType::New();
	M_reader->SetFileName( in_colon_mask );
	try
	{
		M_reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	M_ImageType::Pointer colon_mask = M_reader->GetOutput();

	for (int count_sigma = 0; count_sigma < sizeof(sigma)/sizeof(double); count_sigma++) {
	for (int count_alpha = 0; count_alpha < sizeof(alpha)/sizeof(double); count_alpha++) {
	for (int count_beta = 0; count_beta < sizeof(beta)/sizeof(double); count_beta++) {
	for (int count_eta = 0; count_eta < sizeof(eta)/sizeof(double); count_eta++) {

	// ITK Hessian
	T_ImageType::Pointer hess_image;
	if (useMask) {
		typedef itk::HessianRecursiveGaussianImageFilter < M_ImageType, T_ImageType > HessianFilterTypeM;
		HessianFilterTypeM::Pointer hessian = HessianFilterTypeM::New();
		hessian->SetSigma( sigma[count_sigma] );
		hessian->SetInput(colon_mask);
		try {
			hessian->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cerr << excp << std::endl;
		}  
		hess_image = hessian->GetOutput();
	} else {
		typedef itk::HessianRecursiveGaussianImageFilter < F_ImageType, T_ImageType > HessianFilterTypeF;
		HessianFilterTypeF::Pointer hessian = HessianFilterTypeF::New();
		hessian->SetSigma( sigma[count_sigma] );
		hessian->SetInput(input);
		try {
			hessian->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cerr << excp << std::endl;
		}  
		hess_image = hessian->GetOutput();
	}

	// Get eigen values/vectors
	typedef itk::SymmetricEigenAnalysis
	< T_PixelType, E_PixelType, EV_PixelType > EigAnalysisType;
	EigAnalysisType eig;
	eig.SetDimension( 3 );
	eig.SetOrderEigenMagnitudes( true );
	eig.SetOrderEigenValues( true );

	typedef itk::ImageRegionIteratorWithIndex< T_ImageType > T_IteratorType;
	T_IteratorType hess_it( hess_image, hess_image->GetLargestPossibleRegion() );

	// Input iterator
	typedef itk::ImageRegionIteratorWithIndex< F_ImageType > F_IteratorType;
	F_IteratorType input_it(input, input->GetLargestPossibleRegion() );

	
	// Colon mask iterators
	typedef itk::ImageRegionIteratorWithIndex< M_ImageType > M_IteratorType;
	M_IteratorType colon_mask_it(colon_mask, colon_mask->GetLargestPossibleRegion() );
	//M_IteratorType colon_mask_dilated_it(colon_mask_dilated, colon_mask_dilated->GetLargestPossibleRegion() );

	// Eigenvalue images
	O_ImageType::Pointer lambda1 = O_ImageType::New();
	lambda1->SetRegions(input->GetLargestPossibleRegion() );
	lambda1->Allocate();
	lambda1->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer lambda2 = O_ImageType::New();
	lambda2->SetRegions(input->GetLargestPossibleRegion() );
	lambda2->Allocate();
	lambda2->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer lambda3 = O_ImageType::New();
	lambda3->SetRegions(input->GetLargestPossibleRegion() );
	lambda3->Allocate();
	lambda3->SetSpacing(input->GetSpacing());
	
	O_ImageType::Pointer Fa_image = O_ImageType::New();
	Fa_image->SetRegions(input->GetLargestPossibleRegion());
	Fa_image->Allocate();
	Fa_image->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer Fb_image = O_ImageType::New();
	Fb_image->SetRegions(input->GetLargestPossibleRegion());
	Fb_image->Allocate();
	Fb_image->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer Fc1_image = O_ImageType::New();
	Fc1_image->SetRegions(input->GetLargestPossibleRegion());
	Fc1_image->Allocate();
	Fc1_image->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer Fc2_image = O_ImageType::New();
	Fc2_image->SetRegions(input->GetLargestPossibleRegion());
	Fc2_image->Allocate();
	Fc2_image->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer Frut_image = O_ImageType::New();
	Frut_image->SetRegions(input->GetLargestPossibleRegion());
	Frut_image->Allocate();
	Frut_image->SetSpacing(input->GetSpacing());

	O_ImageType::Pointer Fcup_image = O_ImageType::New();
	Fcup_image->SetRegions(input->GetLargestPossibleRegion());
	Fcup_image->Allocate();
	Fcup_image->SetSpacing(input->GetSpacing());


	// Define iterators
	typedef itk::ImageRegionIteratorWithIndex< O_ImageType > O_IteratorType;
	O_IteratorType lambda1_it( lambda1, lambda1->GetLargestPossibleRegion() );
	O_IteratorType lambda2_it( lambda2, lambda2->GetLargestPossibleRegion() );
	O_IteratorType lambda3_it( lambda3, lambda3->GetLargestPossibleRegion() );
	O_IteratorType Fa_it( Fa_image, Fa_image->GetLargestPossibleRegion() );
	O_IteratorType Fb_it( Fb_image, Fb_image->GetLargestPossibleRegion() );
	O_IteratorType Fc1_it( Fc1_image, Fc1_image->GetLargestPossibleRegion() );
	O_IteratorType Fc2_it( Fc2_image, Fc2_image->GetLargestPossibleRegion() );
	O_IteratorType Frut_it( Frut_image, Frut_image->GetLargestPossibleRegion() );
	O_IteratorType Fcup_it( Fcup_image, Fcup_image->GetLargestPossibleRegion() );

	int count = 0;
  
  for ( hess_it.GoToBegin(), lambda1_it.GoToBegin(), lambda2_it.GoToBegin(), lambda3_it.GoToBegin(), input_it.GoToBegin(), colon_mask_it.GoToBegin(), Fa_it.GoToBegin(), Fb_it.GoToBegin(), Fc1_it.GoToBegin(), Fc2_it.GoToBegin(), Frut_it.GoToBegin(), Fcup_it.GoToBegin(); 
	  !hess_it.IsAtEnd() && !lambda1_it.IsAtEnd() && !lambda2_it.IsAtEnd() && !lambda3_it.IsAtEnd() && !input_it.IsAtEnd()  && !colon_mask_it.IsAtEnd() && !Fa_it.IsAtEnd() && !Fb_it.IsAtEnd() && !Fc1_it.IsAtEnd() && !Fc2_it.IsAtEnd() && !Frut_it.IsAtEnd() && !Fcup_it.IsAtEnd();
	  ++hess_it, ++lambda1_it, ++lambda2_it, ++lambda3_it, ++input_it, ++colon_mask_it, ++Fa_it, ++Fb_it, ++Fc1_it, ++Fc2_it, ++Frut_it, ++Fcup_it)
    {
		if (useMask && colon_mask_it.Get() < 0) {
			lambda1_it.Set(0);
			lambda2_it.Set(0);
			lambda3_it.Set(0);
			Fa_it.Set(0);
			Fb_it.Set(0);
			Fc1_it.Set(0);
			Fc2_it.Set(0);
			Frut_it.Set(0);
			Fcup_it.Set(0);
		} else {
			T_ImageType::IndexType index;
			T_ImageType::PointType point;
			T_PixelType            hess_matrix;
			E_PixelType            gradient_values;
			E_PixelType            eigen_values;
			EV_PixelType           eigen_vectors;

			index = hess_it.GetIndex();
			hess_image->TransformIndexToPhysicalPoint( index, point );
			hess_matrix = hess_it.Get();
			eig.ComputeEigenValuesAndVectors( hess_matrix, eigen_values, eigen_vectors );
			
			//Sort eigenvalues in ascending absolute value using a vector
			std::vector<double> ev_vector(eigen_values.Begin(), eigen_values.End());
			sort(ev_vector.begin(), ev_vector.end(), compare);

			/*
			std::cout 
				<< "Point: (" << point[0] << "," << point[1] << "," << point[2] << ") "
				<< ev_vector[0] << " " <<  ev_vector[1] << " " << ev_vector[2] << std::endl;
			 */

			lambda1_it.Set(ev_vector[0]);
			lambda2_it.Set(ev_vector[1]);
			lambda3_it.Set(ev_vector[2]);
			
			/*// Test eigenvalues
			if ((abs(ev_vector[0]) < 0.5) && (ev_vector[2]>0)  && (abs(ev_vector[2]) > abs(3*ev_vector[1]))) {//( abs(ev_vector[2]) > abs(5*ev_vector[1] && (abs(ev_vector[2] < abs(10*ev_vector[1]))))) {
				Fa_it.Set(1);
			} else {
				Fa_it.Set(0);
			}
			*/

			double fa = functionA(ev_vector,alpha[count_alpha]);
			double fb = functionB(ev_vector,beta[count_beta],gamma);
			double fc1 = functionC(ev_vector[0], ev_vector[1], eta[count_eta]);
			double fc2 = functionC(ev_vector[1], ev_vector[2], eta[count_eta]);

			Fa_it.Set(fa);
			Fb_it.Set(fb);
			Fc1_it.Set(fc1);
			Fc2_it.Set(fc2);
			Frut_it.Set(fa*fb);
			Fcup_it.Set(fc1*fc2);
		}
    }

  std::cout << "Writing files." << std::endl;
  std::stringstream ss;
 
  ss << "lambda1_sigma" << sigma[count_sigma] << ".hdr";
  WriteITK(lambda1, ss.str().c_str());

  ss.str("");
  ss << "lambda2_sigma" << sigma[count_sigma] << ".hdr";
  WriteITK(lambda2, ss.str().c_str());

  ss.str("");
  ss << "lambda3_sigma" << sigma[count_sigma] << ".hdr";
  WriteITK(lambda3, ss.str().c_str());
 
  ss.str("");
  ss << "Fa_sigma" << sigma[count_sigma] << "_alpha_" << alpha[count_alpha] << ".hdr";
  WriteITK(Fa_image, ss.str().c_str());

  ss.str("");
  ss << "Fb_sigma" << sigma[count_sigma] << "_beta_" << beta[count_beta] << "_gamma_" << gamma << ".hdr";
  WriteITK(Fb_image, ss.str().c_str());

  WriteITK(Fc1_image, "Fc1.hdr");
  WriteITK(Fc2_image, "Fc2.hdr");
  WriteITK(Fcup_image, "Fcup.hdr");
  WriteITK(Frut_image, "Frut.hdr"); 

  // Move images to output folder
  std::stringstream folder;
  folder << outputPrefix.str() << "_sigma_" << sigma[count_sigma] << "_alpha_" << alpha[count_alpha] << "_beta_" << beta[count_beta] << "_gamma_" << gamma << "_eta_" << eta[count_eta];

  ss.str("");
  ss << "mkdir " << folder.str();
  system(ss.str().c_str());
  
  ss.str("");
  ss << "move lambda*.hdr " << folder.str();
  system(ss.str().c_str());

  ss.str("");
  ss << "move F*.hdr " << folder.str();
  std::cout << ss.str() << std::endl;
  system(ss.str().c_str());

  ss.str("");
  ss << "move lambda*.img " << folder.str();
  system(ss.str().c_str());

  ss.str("");
  ss << "move F*.img " << folder.str();
  system(ss.str().c_str());

  }
  }
  }
  }
  
  system("PAUSE");
  return EXIT_SUCCESS;
}

bool compare (double a, double b) {
	return ( abs(a) < abs(b) );
}

double functionA(std::vector<double>& lambda, double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs((lambda[1]*lambda[2])));
	//Ra = abs(lambda[0])*sqrt(abs(lambda[1]))/sqrt(abs(lambda[2]));	//new Ra to weigh lambda 3 heavily

	//std::cout << "Ra: " << Ra << std::endl;
	return exp(-vnl_math_sqr(Ra)/(2*vnl_math_sqr(alpha)));
}

double functionB(std::vector<double>& lambda, double beta, double gamma) {
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

double f_rut(std::vector<double>& lambda, double alpha, double beta, double gamma) {
	//double fa = functionA(lambda, alpha);
	//double fb = functionB(lambda, beta, gamma);
	//std::cout << "fa: " << fa << std::endl;
	//std::cout << "fb: " << fb << std::endl;
	
	//double frut = fa*fb;
	//std::cout << "frut " << frut << std::endl;

	//return frut;
	return functionA(lambda, alpha)*functionB(lambda, beta, gamma);
}

double f_cup(std::vector<double>& lambda, double eta) {
	return functionC(lambda[0], lambda[1], eta)*functionC(lambda[1], lambda[2], eta);
}

void WriteITK(O_ImageType::Pointer image, const char * name) {
	typedef itk::ImageFileWriter< O_ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(name);
	writer->SetInput( image );

	try {
		writer->Update();
	} catch( itk::ExceptionObject & excp ) {
		std::cerr << excp << std::endl;
	}  
}