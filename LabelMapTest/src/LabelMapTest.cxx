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
#include "itkOrientImageFilter.h"

template< class T >
typename T::Pointer AllocateNewImage( typename T::RegionType region, typename T::SpacingType spacing );

typedef short													PixelType;
typedef int														ScalarPixelType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	

typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;

template< class T >
void WriteITK ( typename T::Pointer image , std::string ss );
int writeCount = 1;
 
int main() 
{ 
	// Read input
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( "C:/ImageData/mr10_092_13p.i0344_airThreshold_relabel.hdr" );

	try
	{
		reader->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer input = reader->GetOutput();
	ScalarImageType3D::RegionType region = input->GetLargestPossibleRegion();
	ScalarImageType3D::SizeType size = region.GetSize();
	ScalarImageType3D::SpacingType spacing = input->GetSpacing();
	ScalarIteratorType inputIt( input, input->GetLargestPossibleRegion() );

	//WriteITK <ScalarImageType3D> (input, "input.hdr");

	/*itk::OrientImageFilter< ScalarImageType3D, ScalarImageType3D>::Pointer orienter = itk::OrientImageFilter< ScalarImageType3D, ScalarImageType3D>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS);
    orienter->SetInput(input);
    orienter->Update();
	input = orienter->GetOutput();
	
	WriteITK <ScalarImageType3D> (input, "inputOrient.hdr");*/

	// Convert label image to label map
	typedef unsigned long LabelType;	
	typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
	typedef itk::LabelMap< LabelObjectType > LabelMapType;

	typedef itk::LabelImageToShapeLabelMapFilter< ScalarImageType3D, LabelMapType > 	ConverterType;
	ConverterType::Pointer converter = ConverterType::New();
	converter->SetInput( input );
	//converter->SetComputePerimeter( true );
	converter->Update();

	// Setup text file
	std::ofstream myfile;
	myfile.open("shapeParams.txt");

	// Get attributes of label map
	LabelMapType::Pointer labelMap = converter->GetOutput();
	myfile << "#\tSize\tRegionElongation\tSizeRegionRatio\tBinaryPrincipalMoments\tBinaryPrincipalAxes\tBinaryElongation\tEquivalentRadius\tEquivalentPerimeter\tEquivalentEllipsoidSize\tBinaryFlatness\n";

	float lungs[4] = {180, 150, 340, 170};
	float tol = .05;

	for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
	{
		const LabelObjectType * lo = labelMap->GetLabelObject( label );
		
		// Centroid is returned [X, Z, Y]
		myfile << label << "\t" << lo->GetPhysicalSize() << "\t";
		//std::cout << label << "\t" << lo->Centroid() << "\t";
		myfile << lo->GetRegionElongation() << "\t";
		myfile << lo->GetSizeRegionRatio() << "\t";
		//myfile << lo->GetBinaryPrincipalMoments() << "\t";
		//myfile << lo->GetBinaryPrincipalAxes() << "\t";
		myfile << lo->GetBinaryElongation() << "\t";
		myfile << lo->GetEquivalentRadius() << "\t";
		myfile << lo->GetEquivalentPerimeter() << "\t";
		//myfile << lo->GetEquivalentEllipsoidSize() << "\t";
		myfile << lo->GetBinaryFlatness() << "\n";



		//std::cout << labelObject->GetCentroid()[0] / spacing[0] << ", ";
		//std::cout << size[1] - (labelObject->GetCentroid()[2] / spacing[1]) << "]\t" << std::endl;
		//<< labelObject->GetCentroid() << "\t" << labelObject->GetRoundness() << std::endl;

		/*
		float centroidX = labelObject->GetCentroid()[0] / spacing[0];
		float centroidY = size[1] - (lo->GetCentroid()[2] / spacing[1]);

		for (int i=0; i < 2; i++) {
			if ( (centroidX > (1-tol)*lungs[2*i]) && (centroidX < (1+tol)*lungs[2*i]) && (centroidY > (1-tol)*lungs[2*i+1]) && (centroidY < (1+tol)*lungs[2*i+1]) ) {
				std::cout << label << "\t" << lo->GetPhysicalSize() << "\t[" << centroidX << ", " << centroidY << "]" << std::endl;
			}
		}
		*/
	}


	//// convert label map back into image
	//typedef itk::LabelMapToBinaryImageFilter <LabelMapType, ScalarImageType3D> LabelMapConvertType;
	//LabelMapConvertType::Pointer mapconvert = LabelMapConvertType::New();
	//mapconvert->SetInput( labelMap );
	//mapconvert->Update();
	//WriteITK <ScalarImageType3D> (mapconvert->GetOutput(), "mapconvert.hdr");
	
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

// Allocate new image with proper region, spacing, and orientation
template< class T >
typename T::Pointer AllocateNewImage(typename T::RegionType region, typename T::SpacingType spacing ) {
	typename T::Pointer image = typename T::New();
	image->SetRegions( region );
	image->SetSpacing( spacing );
	image->Allocate();

	itk::OrientImageFilter< T, T>::Pointer orienter = itk::OrientImageFilter< T, T>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS);
    orienter->SetInput(image);
    orienter->Update();
	image = orienter->GetOutput();
	return image;
}
