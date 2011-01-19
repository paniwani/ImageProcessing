#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <iostream>
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkNeighborhoodIterator.h>
#include <time.h>
#include <itkFlatStructuringElement.h>

typedef short													PixelType;
typedef int														ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	
typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;

typedef itk::NeighborhoodIterator< ScalarImageType3D > NeighborhoodIteratorType;

typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;

typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
typedef itk::BinaryMedianImageFilter< ScalarImageType3D, ScalarImageType3D> BinaryMedianFilterType;


template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );
int writeCount = 1;
 
int main() 
{ 
	// Read input
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/2_airThresholdSmooth.hdr" );
	reader->Update();

	ScalarImageType3D::Pointer input = reader->GetOutput();

	// Set the clock
	clock_t start = clock();
	
	/*ImageType3D::RegionType region = input->GetLargestPossibleRegion();
	ImageType3D::IndexType endIndex = region.GetIndex();
	ImageType3D::IndexType startIndex = region.GetIndex();	
	endIndex[0]+=(region.GetSize()[0]-1);
	endIndex[1]+=(region.GetSize()[1]-1);
	endIndex[2]+=(region.GetSize()[2]-1);

	ScalarImageType3D::Pointer output = ScalarImageType3D::New();
	output->SetRegions( region );
	output->Allocate();

	ScalarIteratorType inputIt( input, region );
	ScalarIteratorType outputIt( output, region );

	for(inputIt.GoToBegin(), outputIt.GoToBegin();
		!inputIt.IsAtEnd() && !outputIt.IsAtEnd();
		++inputIt, ++outputIt) {

		outputIt.Set( inputIt.Get() );

	}

	for(inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt){
		ScalarImageType3D::IndexType index = inputIt.GetIndex();
		
		if ( inputIt.Get() == 1) {
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							//for (int k=-1;k<=1;k++) {
							//	if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
							ScalarImageType3D::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]};
									
									if (input->GetPixel( neighbor_index ) == 0) {
										output->SetPixel( neighbor_index, 1 );
									}
				
								//	}
								//}
							
						}
					}
				}
			}
		}
	}

	*/

	// Create binary ball structuring element
	//typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
	typedef itk::FlatStructuringElement< 3 > FlatStructuringElementType;
	FlatStructuringElementType flatse;
    flatse.SetRadius( 1 );
    //flatse.CreateStructuringElement();

	// Dilate
	typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, FlatStructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( input );
	dilateFilter->SetKernel( flatse );
	dilateFilter->SetDilateValue( 1 );

	dilateFilter->Update();

	WriteITK <ScalarImageType3D> ( dilateFilter->GetOutput() , "segmentedColon_4px.hdr");



	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);

	WriteITK <ScalarImageType3D> ( dilateFilter->GetOutput(), "output.hdr");
	

	system("pause"); 
	return 0; 
} 


template< class T >
void WriteITK (typename T::Pointer image , std::string ss )
{
	typedef itk::ImageFileWriter< T >	WriterType;
	WriterType::Pointer writer = WriterType::New();

	std::stringstream ss2;
	ss2 << writeCount++ << "_" << ss;
	writer->SetFileName(ss2.str().c_str());

	writer->SetInput(image);
	writer->GlobalWarningDisplayOff();
	writer->ReleaseDataFlagOn();
	std::cerr<<"Writing: "<< ss <<std::endl;
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return;
	}
}
