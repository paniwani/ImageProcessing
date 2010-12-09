/*
Neil Panjwani
11/15/2010
Computes Hessian Response across sigma scales
*/
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
#include "itkOrientImageFilter.h"
#include "NIH_ITK_Utility.h"
#include <CIS_Image_Processing_Algo_3D.h>
#include <CIS_Array_Image3D.h>
#include <string>

typedef itk::Image< float, 3 > ImageType;

double fnA(std::vector<double>& lambda, double alpha);
double fnB(std::vector<double>& lambda, double beta, double gamma);
double fnC(double ev1, double ev2, double eta);
double fnRut(std::vector<double>& lambda, double alpha, double beta, double gamma);
double fnCup(std::vector<double>& lambda, double eta);
bool compare (double a, double b);
void WriteITK(ImageType::Pointer image, const char * name);
void WriteCIS(ImageType::Pointer image, std::string str);

int main(int argc, char * argv[])
{
	/****************INPUTS**********************
	argv[1]: File name of input image
	argv[2]: File name of output hessian response
	********************************************/
	
	char * in_file_name = argv[1];

	// ITK Pixel Types
	typedef itk::SymmetricSecondRankTensor< float, 3 > TensorPixelType;
	typedef itk::Vector< float, 3> VectorPixelType;
	typedef itk::Vector< VectorPixelType, 3> EigenMatrixPixelType;


	// ITK Image Types
	typedef itk::Image< TensorPixelType,  3 > TensorImageType;

	// Read Input
	typedef itk::ImageFileReader< ImageType > ReaderType;
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

	ImageType::Pointer input = reader->GetOutput();

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	/*
	// Select tagged regions only
	
	IteratorType input_it( input, input->GetLargestPossibleRegion() );

	for ( input_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it)
	{
		if (input_it.Get() < 200) { input_it.Set(0); }
	}

	WriteITK(input,"tagged.hdr");
	*/

	ImageType::SpacingType spacing = input->GetSpacing();

	// Define parameters
	//double sigma[2]		= {0.8, 1.6};
	double sigma;
	double alpha		= .5;
	double beta			= .3;
	double gamma		= .3;
	double eta			= .2;

	int count = 1;
	// Loop through parameters
	//for (int count_sigma = 0; count_sigma < sizeof(sigma)/sizeof(double); count_sigma++)	{
	
	for (double sigmaFactor = .5; sigmaFactor <= 8 ; sigmaFactor+=.5) {
		std::cout << "Loop " << count << " of 16" << std::endl; 

		// Set sigma
		sigma = sigmaFactor*spacing[0];

		// Setup images
		ImageType::Pointer H = ImageType::New();
		H->SetRegions(input->GetLargestPossibleRegion() );
		H->Allocate();
		H->SetSpacing(input->GetSpacing());

		IteratorType H_it( H, H->GetLargestPossibleRegion() );
		
		// Compute Hessian
		std::cout << "Computing Hessian" << std::endl;
		typedef itk::HessianRecursiveGaussianImageFilter < ImageType, TensorImageType > HessianFilterType;
		HessianFilterType::Pointer hessian = HessianFilterType::New();
		hessian->SetSigma( sigma );
		hessian->SetInput(input);
		try 
		{
			hessian->Update();
		} catch( itk::ExceptionObject & excp ) 
		{
			std::cerr << excp << std::endl;
		}  
		TensorImageType::Pointer hess_image = hessian->GetOutput();
		
		// Get eigenvalues
		std::cout << "Computing Eigenvalues" << std::endl;
		typedef itk::SymmetricEigenAnalysis< TensorPixelType, VectorPixelType, EigenMatrixPixelType > EigAnalysisType;
		EigAnalysisType eig;
		eig.SetDimension( 3 );
		eig.SetOrderEigenMagnitudes( true );
		eig.SetOrderEigenValues( true );

		// Iterators
		typedef itk::ImageRegionIteratorWithIndex< TensorImageType > TensorIteratorType;
		TensorIteratorType hess_it( hess_image, hess_image->GetLargestPossibleRegion() );
		
		std::cout << "Computing Response" << std::endl;

		for (	hess_it.GoToBegin(), H_it.GoToBegin();
				!hess_it.IsAtEnd() && !H_it.IsAtEnd();
				++hess_it, ++H_it )
		{
			// Box region
			ImageType::IndexType index = H_it.GetIndex();
			if (index[0] > 65 && index[0] < 430 && index[1] < 512-160 && index[1] > 512-380)
			{
				//std::cout << index[0] << " " << index[1] << " " << index[2] << std::endl;

				// Compute eigenvalues
				VectorPixelType eigen_values;
				eig.ComputeEigenValues( hess_it.Get(), eigen_values );

				// Sort eigenvalues
				std::vector<double> lambda(eigen_values.Begin(), eigen_values.End());
				sort(lambda.begin(), lambda.end(), compare);

				// Compute hessian response
				double frut = fnRut(lambda, alpha, beta, gamma);
				double fcup = fnCup(lambda, eta);
		
				// Set max value across scales
				double max = frut>fcup ? frut : fcup;
				//if (max > H_it.Get()) {	H_it.Set(max); }
				H_it.Set(max);
			} else {
				H_it.Set(0);
			}
		}
		/*
		// Orient all input images into RAS orientation
		itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
		orienter->UseImageDirectionOn();
		orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR);
		orienter->SetInput(H);
		try {
			orienter->Update();
		} catch( itk::ExceptionObject & excp ) 
		{
			std::cerr << excp << std::endl;
		}*/

		// Write H response
		std::stringstream ss;
		ss << "H_sigmaFactor_" << sigmaFactor << ".hdr";
		WriteCIS(H, ss.str());
		//WriteITK(H, ss.str().c_str());

		std::cout << "End Loop " << count++ << " of 16" << std::endl << std::endl;
	}

	system("pause");
	return 0;
}

bool compare (double a, double b) {
	return ( abs(a) < abs(b) );
}

double fnA(std::vector<double>& lambda, double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs((lambda[1]*lambda[2])));
	return exp(-vnl_math_sqr(Ra)/(2*vnl_math_sqr(alpha)));
}

double fnB(std::vector<double>& lambda, double beta, double gamma) {
	double Rb;
	Rb = abs(lambda[1])/abs(lambda[2]);
	return exp(-vnl_math_sqr(Rb - gamma)/(2*vnl_math_sqr(beta)));
}

double fnC(double ev1, double ev2, double eta) {
	double Rc;
	Rc = abs(ev1)/abs(ev2);
	return 1.0 - exp(-vnl_math_sqr(Rc)/(2*vnl_math_sqr(eta)));
}

double fnRut(std::vector<double>& lambda, double alpha, double beta, double gamma) {
	if (lambda[2] > 0)
	{
		return fnA(lambda, alpha)*fnB(lambda, beta, gamma);
	} else {
		return 0;
	}
}

double fnCup(std::vector<double>& lambda, double eta) {
	if (lambda[2] > 0)
	{
		return fnC(lambda[0], lambda[1], eta)*fnC(lambda[1], lambda[2], eta);
	} else {
		return 0;
	}
}

void WriteITK(ImageType::Pointer image, const char * name) {
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(name);
	writer->SetInput(image);

	try  
	{
		writer->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	} 
}

void WriteCIS(ImageType::Pointer image, std::string str) {
	//Use CIS to write image
	CIS_Array_Image3D_short *output_CIS;
	output_CIS = new CIS_Array_Image3D_short();

	//Convert ITK image back into CIS image
	std::cout << "Converting ITK to CIS..." << std::endl;
	Copy_ITKImage_to_CISImage(image, output_CIS);

	char * cstr = new char [str.size()+1];
	strcpy(cstr, str.c_str());

	Write_Analyze_File(cstr, *output_CIS);
	std::cout << "Done." << std::endl;
}