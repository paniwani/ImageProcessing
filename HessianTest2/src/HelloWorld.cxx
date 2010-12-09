#include "itkImage.h"
#include <iostream>
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSymmetricSecondRankTensor.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkFixedArray.h>
#include <stdlib.h>
#include <vector>

double functionA(std::vector<double>& lambda, double alpha);
double functionB(std::vector<double>& lambda, double beta, double gamma);
double functionC(double ev1, double ev2, double eta);
double f_rut(std::vector<double>& lambda, double alpha, double beta, double gamma);
double f_cup(std::vector<double>& lambda, double eta);
bool compare (double a, double b);

int main( int argc, char *argv[] )
{
	if ( argc < 2 ) {
		std::cerr << "Usage: " << std::endl;
		std::cerr << "HessianTest2 input outputImage" << std::endl;
		system("PAUSE");
		return EXIT_FAILURE;
	}

	typedef float														PixelType;
	typedef float														OutputPixelType;	
	const unsigned int													Dimension = 3;
    typedef itk::Image< PixelType, Dimension >							ImageType;
	typedef itk::Image< OutputPixelType, Dimension >					OutputImageType;
	typedef itk::HessianRecursiveGaussianImageFilter< ImageType >		HessianFilterType;
	typedef itk::SymmetricSecondRankTensor< double >					TensorType;
	typedef itk::Image<TensorType, Dimension>							TensorImageType;
	typedef itk::ImageFileReader< ImageType       >						ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >						WriterType;
	typedef itk::ImageRegionIteratorWithIndex<ImageType>				IteratorType;
	typedef itk::ImageRegionIteratorWithIndex<TensorImageType>			IteratorTensorType;
	typedef itk::FixedArray<double>										EigenValueArrayType;
	
	ReaderType::Pointer  reader  = ReaderType::New();
	WriterType::Pointer  writer  = WriterType::New();

	//Set params
	double alpha = 0.1;
	double beta = 0.3;
	double gamma = 0.3;
	double eta = 0.2;
	
	//Read input image
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

	//double sigma_ar[] = {.0104, .0208, 1, 2, .8203, 1.6406, 2.5, 5};
	const int numOfSigma = 1;
	double sigma_ar[numOfSigma] = {.8203};
	HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
	hessianFilter->SetInput( reader->GetOutput() );

	//Allocate cup and rut images
	std::cout << "Allocating F_rut and F_cup" << std::endl;

	ImageType::Pointer F_rut_image = ImageType::New();
	F_rut_image->SetRegions(input->GetRequestedRegion());
	F_rut_image->Allocate();
	F_rut_image->SetSpacing(input->GetSpacing());

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
	F_cup_image->Allocate();
	F_cup_image->SetSpacing(input->GetSpacing());
	*/

	//Iterate through input image and fcup/frut images
	IteratorType input_iter(input, input->GetRequestedRegion() );
	//IteratorType colon_mask_iter(colon_mask, colon_mask->GetRequestedRegion() );
	//IteratorType colon_mask_shaded_iter(colon_mask_shaded, colon_mask_shaded->GetRequestedRegion() );
	IteratorType F_rut_iter(F_rut_image, F_rut_image->GetRequestedRegion() );
	//IteratorType F_cup_iter(F_cup_image, F_cup_image->GetRequestedRegion() );
	IteratorType F_a_iter(F_a_image, F_a_image->GetRequestedRegion() );
	IteratorType F_b_iter(F_b_image, F_b_image->GetRequestedRegion() );

	TensorImageType::Pointer tensorImage;
	TensorType tensor;
	EigenValueArrayType ev;
	
	for (int i=0; i < numOfSigma; i++) {
		//Apply Hessian Filter
		hessianFilter->SetSigma(.08);//sigma_ar[i]);
		hessianFilter->Update();

		//Get tensor image
		tensorImage = hessianFilter->GetOutput();
		IteratorTensorType tensor_iter( tensorImage, input->GetRequestedRegion() );

		std::cout << "Computing eigenvalues" << std::endl;

		for (	input_iter.GoToBegin(), tensor_iter.GoToBegin(), F_rut_iter.GoToBegin(), F_a_iter.GoToBegin(), F_b_iter.GoToBegin();
				!input_iter.IsAtEnd() && !tensor_iter.IsAtEnd() && !F_rut_iter.IsAtEnd() && !F_a_iter.IsAtEnd() && !F_b_iter.IsAtEnd();
				++input_iter, ++tensor_iter, ++F_rut_iter, ++F_a_iter, ++F_b_iter)
		{
			//if (input_iter.Get() < 200 && input_iter.Get() > -500) {	//threshold
				
			if (input_iter.Get() == 40) { //input is colon mask, look only inside
				//Solve eigensystem
				tensor = tensor_iter.Get();
				tensor.ComputeEigenValues(ev);
				
				//Sort eigenvalues in ascending absolute value using a vector
				std::vector<double> ev_vector(ev.Begin(), ev.End());
				sort(ev_vector.begin(), ev_vector.end(), compare);

				//std::cout << "lambda: " << ev_vector[0] << ", " << ev_vector[1] << ", " << ev_vector[2] << std::endl;

				if (ev_vector[2] > 0) {
					//Solve functions
					//F_rut_iter.Set(f_rut(ev_vector, alpha, beta, gamma));
					//F_cup_iter.Set(f_cup(ev_vector, eta));
					F_a_iter.Set(functionA(ev_vector, alpha));
				} else {
					F_a_iter.Set(0);
				}
			
			} else {
				//F_rut_iter.Set(0);
				//F_cup_iter.Set(0);
				F_a_iter.Set(0);
			}		
		}

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
		ss << "F_a" << "_sigma" << sigma_ar[i] << "_alpha " << alpha << ".hdr";
		writer->SetFileName(ss.str().c_str());
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

	system("PAUSE");
    return EXIT_SUCCESS;
}

bool compare (double a, double b) {
	return ( abs(a) < abs(b) );
}

double functionA(std::vector<double>& lambda, double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs(lambda[1]*lambda[2]));
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
	double fa = functionA(lambda, alpha);
	double fb = functionB(lambda, beta, gamma);
	//std::cout << "fa: " << fa << std::endl;
	//std::cout << "fb: " << fb << std::endl;
	
	double frut = fa*fb;
	//std::cout << "frut " << frut << std::endl;

	return frut;
	//return functionA(lambda, alpha)*functionB(lambda, beta, gamma);
}

double f_cup(std::vector<double>& lambda, double eta) {
	return functionC(lambda[0], lambda[1], eta)*functionC(lambda[1], lambda[2], eta);
}