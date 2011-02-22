#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkImageRegionIterator.h";

int main( int argc, char * argv[] )
{

	typedef    float    InputPixelType;
	typedef    float    OutputPixelType;

	typedef itk::Image< InputPixelType,  2 >   InputImageType;
	typedef itk::Image< OutputPixelType, 2 >   OutputImageType;

	typedef itk::ImageFileReader< InputImageType >  ReaderType;

	typedef itk::ImageRegionIterator< OutputImageType > IteratorType;

	typedef itk::GradientAnisotropicDiffusionImageFilter<
			   InputImageType, OutputImageType >  FilterType;
	FilterType::Pointer filter = FilterType::New();

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( "C:/GitProjects/GradientAnisotropicDiffusionImageFilter/build64-3.16.0/mr10_316_13s.i0355_1_input_sample2.hdr" );

	filter->SetInput( reader->GetOutput() );

	const unsigned int numberOfIterations = 5;

	const double       timeStep = 0.125;

	filter->SetNumberOfIterations( numberOfIterations );
	filter->SetTimeStep( timeStep );



	const double       conductance[7] = {0,.5,1,2,5,10,100};

	for (int th = 200; th <= 800; th+=100) {
	for (int i = 0; i < 7; i++) {

		filter->SetConductanceParameter( conductance[i] );

		filter->Update();

		IteratorType it( filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );

		for ( it.GoToBegin(); !it.IsAtEnd() ; ++it ) 
		{
			if ( it.Get() > th )	{ it.Set( 1 ); }
			else					{ it.Set( 0 ); }
		}

		//typedef unsigned char WritePixelType;
		//typedef itk::Image< WritePixelType, 2 > WriteImageType;
		typedef itk::Image< OutputPixelType, 2 > WriteImageType;
		
		/*typedef itk::RescaleIntensityImageFilter< 
				   OutputImageType, WriteImageType > RescaleFilterType;

		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
		rescaler->SetOutputMinimum(   0 );
		rescaler->SetOutputMaximum( 255 );*/

		typedef itk::ImageFileWriter< WriteImageType >  WriterType;

		WriterType::Pointer writer = WriterType::New();
		
		std::stringstream ss;
		ss << "sample2_threshold_" << th << "_iter_" << numberOfIterations << "_time_0.125_conductance_" << conductance[i] << ".hdr";

		writer->SetFileName( ss.str() );

		//rescaler->SetInput( filter->GetOutput() );
		writer->SetInput( filter->GetOutput() );
		writer->GlobalWarningDisplayOff();
		writer->Update();

	}
	}


	system("pause");
	return EXIT_SUCCESS;
}

