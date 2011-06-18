const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkAdaptiveOtsuThresholdImageFilterModified.h>
#include <itkImageRegionMultidimensionalSplitter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>

#include "itkVector.h"
#include "itkPointSet.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRandomNonRepeatingConstIteratorWithIndex.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
 												
int main(int argc, char * argv[])				
{ 					
	/*unsigned int splineOrder = 3;
	unsigned int numberOfLevels = 3;
	unsigned int numberOfControlPoints = 20;*/

	///////////////////////////////////////////////
	// Load image
	///////////////////////////////////////////////
	ImageType::Pointer input = ReadITK <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm/mr10_092_13p_i0140.dcm");
	ImageType::RegionType region = input->GetLargestPossibleRegion();
		
	///////////////////////////////////////////////
	// Segment colon
	///////////////////////////////////////////////
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ImageType> (input,colonRegion);
	//Mask <ImageType,ByteImageType> (input,colon,-1025);

	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	/////////////////////////////////////////////////
	//// Setup otsu calculator and compute global threshold
	/////////////////////////////////////////////////
	//typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	//OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	//otsuCalculator->SetImage( input );
	//otsuCalculator->SetMinMax(true);
	//otsuCalculator->SetHistogramMin(-300);
	//otsuCalculator->SetHistogramMax(1500);
	//otsuCalculator->SetNumberOfHistogramBins(255);
	//otsuCalculator->Compute();
	//PixelType otsuGlobal = otsuCalculator->GetThreshold();

	//std::cout << "Otsu global: " << otsuGlobal << std::endl;

	/////////////////////////////////////////////////
	//// Divide image into sub regions and generate pointset
	/////////////////////////////////////////////////
	//typedef itk::ImageRegionMultidimensionalSplitter<Dimension> SplitterType;
	//SplitterType::Pointer splitter = SplitterType::New();
	//
	//unsigned int requestedNumOfSplits = 5000;
	//unsigned int numOfSplits = splitter->GetNumberOfSplits(region,requestedNumOfSplits);	
	//
	//typedef ImageType::PointType      InputPointType;
	//typedef InputPointType::CoordRepType   InputCoordType;
	//
	//typedef itk::Vector< InputCoordType, 1 >         VectorType;
	//typedef itk::Image< VectorType, Dimension > VectorImageType;
	//typedef VectorImageType::PixelType VectorPixelType;

	//typedef itk::PointSet< VectorPixelType, Dimension >      PointSetType;
	//typedef PointSetType::Pointer                   PointSetPointer;
	//typedef PointSetType::PointType                 PointSetPointType;
	//typedef PointSetType::PointsContainerPointer    PointsContainerPointer;
	//typedef PointSetType::PointDataContainerPointer PointDataContainerPointer;
	//typedef PointSetType::PointDataContainer        PointDataContainer;

	//

	//PointSetType::Pointer pointSet = PointSetType::New();
	//PointsContainerPointer pointsContainer = pointSet->GetPoints();
	//pointsContainer->Reserve( numOfSplits );

	//PointDataContainerPointer pointDataContainer = PointDataContainer::New();
	//pointDataContainer->Reserve( numOfSplits );
	//pointSet->SetPointData( pointDataContainer );

	//// For each subregion, compute local threshold and add it to pointset
	//for (int i=0; i<numOfSplits; i++)
	//{
	//	// Compute center pixel
	//	ImageType::RegionType localRegion = splitter->GetSplit(i,numOfSplits,region);
	//	
	//	ImageType::IndexType localIndex = localRegion.GetIndex();
	//	ImageType::SizeType localSize = localRegion.GetSize();
	//	
	//	ImageType::IndexType centerIndex;

	//	for (int j=0; j<Dimension; j++)
	//	{
	//		centerIndex[j] = localIndex[j] + (long) (localSize[j]/2);
	//	}

	//	// Add index to point set
	//	PointSetPointType point;

	//	input->TransformIndexToPhysicalPoint(centerIndex,point);

	//	pointsContainer->SetElement(i,point);

	//	// Compute local otsu
	//	otsuCalculator->SetRegion(localRegion);
	//	otsuCalculator->Compute();
	//	PixelType otsuLocal = otsuCalculator->GetThreshold();

	//	// Validate threshold bounds
	//	if (otsuLocal < 150 || otsuLocal > 650)
	//	{
	//		otsuLocal = otsuGlobal;
	//	}

	//	// Add threshold data to pointset
	//	VectorPixelType V;
	//	V[0] = static_cast<InputCoordType>( otsuLocal );
	//	pointDataContainer->SetElement( i, V );
	//}

	//// bspline
	//typedef itk::BSplineScatteredDataPointSetToImageFilter< PointSetType,VectorImageType > SDAFilterType;
	//typedef SDAFilterType::Pointer SDAFilterPointer;

	//typedef itk::Image< InputCoordType, Dimension > CoordImageType;
	//typedef CoordImageType::Pointer        CoordImagePointer;
	//typedef itk::VectorIndexSelectionCastImageFilter< VectorImageType,ByteImageType > IndexFilterType;
	//typedef IndexFilterType::Pointer IndexFilterPointer;

	//SDAFilterType::ArrayType ncps;
	//ncps.Fill( numberOfControlPoints );

	//SDAFilterPointer filter = SDAFilterType::New();
	//filter->SetSplineOrder( splineOrder );
	//filter->SetNumberOfControlPoints( ncps );
	//filter->SetNumberOfLevels( numberOfLevels );

	//// Define the parametric domain.
	//filter->SetOrigin( input->GetOrigin() ); 
	//filter->SetSpacing( input->GetSpacing() );
	//filter->SetSize( region.GetSize() );
	//filter->SetInput( pointSet );
	//filter->Update();

	//IndexFilterPointer componentExtractor = IndexFilterType::New();
	//componentExtractor->SetInput( filter->GetOutput() );
	//componentExtractor->SetIndex( 0 );
	//componentExtractor->Update();
	//ByteImageType::Pointer thresholdImage = componentExtractor->GetOutput();

	//WriteITK <ByteImageType> (thresholdImage,"thresholdImage.nii");


	/////////////////////////////////////////////////
	//// Get local threshold
	/////////////////////////////////////////////////

	typedef itk::AdaptiveOtsuThresholdImageFilterModified<ImageType,ByteImageType> AdaptiveThresholdFilterType;
	AdaptiveThresholdFilterType::Pointer thresholder = AdaptiveThresholdFilterType::New();
	thresholder->SetInput(input);
	thresholder->SetInsideValue(255);
	
	/*
	thresholder->SetOutsideValue(0);
	thresholder->SetNumberOfHistogramBins(255);
	thresholder->SetNumberOfControlPoints(30);
	thresholder->SetNumberOfLevels(4);
	thresholder->SetNumberOfSamples(5000);
	*/

	thresholder->Update();

	WriteITK <ByteImageType> (thresholder->GetOutput(),"threshold.nii");

	// View the threshold image
	WriteITK <ByteImageType> (thresholder->GetThresholdImage(),"thresholdImage.nii");



	//// View the pointset
	//typedef itk::PointSetToImageFilter<AdaptiveThresholdFilterType::PointSetType, AdaptiveThresholdFilterType::VectorImageType>  ConverterType;
	//ConverterType::Pointer converter = ConverterType::New();
	//converter->SetInput( thresholder->GetPointSet() );
	//
	//typedef itk::VectorIndexSelectionCastImageFilter< AdaptiveThresholdFilterType::VectorImageType,ImageType > IndexFilterType;
	//IndexFilterType::Pointer indexFilter = IndexFilterType::New();
	//indexFilter->SetInput( converter->GetOutput() );
	//indexFilter->SetIndex(0);
	//
	//indexFilter->Update();
	//WriteITK <ImageType> (indexFilter->GetOutput(),"pointset.nii");

	
	return 0;
}