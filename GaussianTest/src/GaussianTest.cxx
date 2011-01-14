#include "itkImage.h" 
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include "itkConnectedThresholdImageFilter.h"

template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );
int writeCount = 1;

typedef short													PixelType;
typedef int														ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	

typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;
 
int main() 
{ 
	// Read airmask
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/4_airThresholdWithoutLungs.hdr" );
	reader->Update();
	ScalarImageType3D::Pointer airMask = reader->GetOutput();
	ScalarIteratorType airMaskIt( airMask, airMask->GetLargestPossibleRegion() );

	WriteITK <ScalarImageType3D> (airMask, "airMask.hdr");

	// Read input CT
	typedef itk::ImageFileReader< ImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();

	reader2->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/1_input.hdr" );
	reader2->Update();
	ImageType3D::Pointer input = reader2->GetOutput();

	// Set regions
	ImageType3D::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType3D::IndexType endIndex = fullRegion.GetIndex();
	ImageType3D::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);


	/*typedef itk::DiscreteGaussianImageFilter< ScalarImageType3D, ScalarImageType3D > DiscreteGaussianFilterType;
	DiscreteGaussianFilterType::Pointer discreteGaussianFilter = DiscreteGaussianFilterType::New();
	discreteGaussianFilter->SetInput( input );
	discreteGaussianFilter->SetVariance( 1 );
	//discreteGaussianFilter->SetMaximumKernelWidth();
	discreteGaussianFilter->Update();

	WriteITK <ScalarImageType3D> (discreteGaussianFilter->GetOutput(), "gaussian.hdr");*/

	// Create binary ball structuring element
	typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 4 );
    structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( airMask );
	dilateFilter->SetKernel( structuringElement );
	dilateFilter->SetDilateValue( 1 );

	try
	{
		dilateFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer dilateAir = dilateFilter->GetOutput();
	ScalarIteratorType dilateAirIt( dilateAir, dilateAir->GetLargestPossibleRegion() );

	WriteITK <ScalarImageType3D> ( dilateAir , "airMaskDilated.hdr");

	// Apply region growing filter to detect tagged regions >= 200
	typedef itk::ConnectedThresholdImageFilter< ImageType3D, ScalarImageType3D > ConnectedThresholdFilterType;
	ConnectedThresholdFilterType::Pointer connectedThresholdFilter = ConnectedThresholdFilterType::New();
	connectedThresholdFilter->SetLower( 200 );
	//connectedThresholdFilter->SetUpper( max );
	connectedThresholdFilter->SetReplaceValue( 1 );

	// Set edges of air as seeds
	for (	airMaskIt.GoToBegin();
			!airMaskIt.IsAtEnd();
			++airMaskIt) 
	{	
		if (airMaskIt.Get() == 1) {	// at air

			ImageType3D::IndexType index = airMaskIt.GetIndex();
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							for (int k=-1;k<=1;k++) {
								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
									ImageType3D::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};

									if ( airMask->GetPixel( neighborIndex ) == 0 )	// if neighbor is non-air
									{
										airMask->SetPixel( index, 2); // mark as edge
										connectedThresholdFilter->AddSeed( index );
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// Set dilated air region as seed in connected threshold
	for (	airMaskIt.GoToBegin(), dilateAirIt.GoToBegin();
			!airMaskIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
			++airMaskIt, ++dilateAirIt		) 
	{	
		if ( airMaskIt.Get() == 0 && dilateAirIt.Get() == 1 ) 
		{
			ImageType3D::IndexType index = airMaskIt.GetIndex();
			airMask->SetPixel( index , 2);
			connectedThresholdFilter->AddSeed( index );
		}
	}

	// Write updated air mask to show region growing seeds
	WriteITK <ScalarImageType3D> ( airMask, "airMaskWithSeeds.hdr");

	connectedThresholdFilter->SetInput( input );

	try
	{
		connectedThresholdFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer tagged = connectedThresholdFilter->GetOutput();
	//ScalarIteratorType taggedIt( tagged, tagged->GetLargestPossibleRegion() );

	WriteITK <ScalarImageType3D> ( tagged, "taggedRegionGrowing.hdr");

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
