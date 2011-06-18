// a test routine for regional extrema using flooding
#include "itkRegionalMaximaImageFilter.h"
#include "itkHConvexImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

int main(int, char * argv[])
{
  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );

  typedef itk::RegionalMaximaImageFilter< IType, IType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetFullyConnected( atoi(argv[1]) );
  filter->SetFlatIsMaxima( atoi(argv[2]) );
  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();


  // produce the same output with other filters
  typedef itk::HConvexImageFilter< IType, IType > ConvexType;
  ConvexType::Pointer convex = ConvexType::New();
  convex->SetInput( reader->GetOutput() );
  convex->SetFullyConnected( atoi(argv[1]) );
  convex->SetHeight( 1 );

  // convex gives maxima with value=1 and others with value=0
  // rescale the image so we have maxima=255 other=0
  typedef itk::RescaleIntensityImageFilter< IType, IType > RescaleType;
  RescaleType::Pointer rescale = RescaleType::New();
  rescale->SetInput( convex->GetOutput() );
  rescale->SetOutputMaximum( 255 );
  rescale->SetOutputMinimum( 0 );

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( rescale->GetOutput() );
  writer2->SetFileName( argv[5] );
  writer2->Update();

  return 0;
}



// const unsigned int Dimension = 2;
// #include <itkImage.h> 						
// #include <iostream> 		
// #include <utils2.h>
// #include <itkValuedRegionalMinimaImageFilter.h>
// #include <itkValuedRegionalMaximaImageFilter.h>
 												
// int main(int argc, char * argv[])				
// { 		
	// ByteImageType::Pointer input = ReadITK <ByteImageType> ("cthead1.png");

	// typedef itk::ValuedRegionalMinimaImageFilter<ByteImageType,ByteImageType> MinimaFilterType;
	// MinimaFilterType::Pointer minimaFilter = MinimaFilterType::New();
	// minimaFilter->SetInput(input);
	// minimaFilter->Update();
	// WriteITK <ByteImageType> (minimaFilter->GetOutput(),"minima.nii");

	// typedef itk::ValuedRegionalMaximaImageFilter<ByteImageType,ByteImageType> MaximaFilterType;
	// MaximaFilterType::Pointer maximaFilter = MaximaFilterType::New();
	// maximaFilter->SetInput(input);
	// maximaFilter->Update();
	// WriteITK <ByteImageType> (maximaFilter->GetOutput(),"maxima.nii");

	//system("pause");
	// return 0; 									
// } 												
