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
	// Start clock
	//clock_t start = clock();

	// Read input
	typedef itk::ImageFileReader< ImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName( "C:/ImageData/mr10_092_13p.i0344.hdr" );

	try
	{
		reader->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ImageType3D::Pointer input = reader->GetOutput();
	ImageType3D::RegionType region = input->GetLargestPossibleRegion();

	ImageType3D::IndexType endIndex = region.GetIndex();
	ImageType3D::IndexType startIndex = region.GetIndex();	
	endIndex[0]+=(region.GetSize()[0]-1);
	endIndex[1]+=(region.GetSize()[1]-1);
	endIndex[2]+=(region.GetSize()[2]-1);

	ScalarImageType3D::SizeType size = region.GetSize();
	ImageType3D::SpacingType spacing = input->GetSpacing();
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );

	WriteITK <ImageType3D> (input, "input.hdr");


	// Test images
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();
	//reader2->SetFileName( "trial1data/4_airThresholdWithoutLungs.hdr" );
	reader2->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/8_stoolThresholdWithoutBoneSmoothed.hdr" );
	reader2->Update();
	ScalarImageType3D::Pointer stoolThresholdSmooth = reader2->GetOutput();

	ReaderType2::Pointer reader3 = ReaderType2::New();
	reader3->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/4_airThresholdWithoutLungs.hdr" );
	reader3->Update();
	ScalarImageType3D::Pointer airThreshold = reader3->GetOutput();

	// Fast Dilate Test

	/*NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	NeighborhoodIteratorType nit( radius, test, test->GetLargestPossibleRegion() );

	clock_t start = clock();

	ScalarImageType3D::Pointer dilateTest = ScalarImageType3D::New();
	dilateTest->SetRegions( region );
	dilateTest->SetSpacing( spacing );
	dilateTest->Allocate();
	ScalarIteratorType dilateTestIt( dilateTest, region );
	
	for (nit.GoToBegin(), dilateTestIt.GoToBegin();
		!nit.IsAtEnd() && !dilateTestIt.IsAtEnd();
		++nit, ++dilateTestIt ) {

		dilateTestIt.Set( nit.GetCenterPixel() );

		if ( nit.GetCenterPixel() == 1 ) { // white
			for (int i = 0; i < nit.Size(); i++) {
				if (nit.GetPixel(i) == 0) {	// black neighbor
					ScalarImageType3D::IndexType neighbor = nit.GetIndex(i);

					if ( neighbor[0]>=startIndex[0] && neighbor[0]<=endIndex[0] &&
						 neighbor[1]>=startIndex[1] && neighbor[1]<=endIndex[1] &&
						 neighbor[2]>=startIndex[2] && neighbor[2]<=endIndex[2] ) {	// within bounds
					
						//std::cout << nit.GetIndex() << std::endl;
						//std::cout << nit.GetIndex(i) << std::endl;
						dilateTest->SetPixel( neighbor , 1 );	// expand white pixels
						break;
					}
				}
			}
		}
	}*/

	/*
	
	for (testIt.GoToBegin(), dilateTestIt.GoToBegin();
		!testIt.IsAtEnd() && !dilateTestIt.IsAtEnd();
		++testIt, ++dilateTestIt ) {
			dilateTestIt.Set( testIt.Get() );

			ScalarImageType3D::IndexType index = testIt.GetIndex();
	
			if ( testIt.Get() == 1 ) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										
										ScalarImageType3D::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]+k};
										if ( test->GetPixel(neighbor_index) == 0) {
											dilateTest->SetPixel( neighbor_index, 1 );
										}
									}
								}
							}
						}
					}
				}
			}
	}

	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	WriteITK <ScalarImageType3D> (dilateTest, "test4.hdr");
	*/

	/*ReaderType2::Pointer reader3 = ReaderType2::New();
	reader3->SetFileName( "trial1data/7_stoolThresholdWithoutBone.hdr" );
	reader3->Update();
	ScalarImageType3D::Pointer stoolThreshold = reader3->GetOutput();
	stoolThreshold->SetSpacing( spacing );*/
	
	//// Threshold to detect air
	//ScalarImageType3D::Pointer airThreshold = ScalarImageType3D::New();
	//airThreshold->SetRegions( region );
	//airThreshold->SetSpacing( spacing );
	//airThreshold->Allocate();
	ScalarIteratorType airThresholdIt( airThreshold, region );

	//for (	inputIt.GoToBegin(), airThresholdIt.GoToBegin();
	//		!inputIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
	//		++inputIt, ++airThresholdIt		) 
	//{
	//	if ( inputIt.Get() <= -600 )	{ airThresholdIt.Set( 1 ); }
	//	else							{ airThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> (airThreshold, "airThreshold.hdr");

	//// Run connected component filter to remove lungs from air threshold
	//ConnectedComponentFilterType::Pointer connectedComponentAirFilter = ConnectedComponentFilterType::New();
	//connectedComponentAirFilter->SetInput( airThreshold );

	//try
	//{
	//	connectedComponentAirFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	////WriteITK <ScalarImageType3D> ( connectedComponentAirFilter->GetOutput(), "connectedComponentAir.hdr");

	//// Relabel components
	//RelabelFilterType::Pointer relabelAirFilter = RelabelFilterType::New();
	//relabelAirFilter->SetInput( connectedComponentAirFilter->GetOutput() );
	//relabelAirFilter->SetMinimumObjectSize( 5000 );

	//try
	//{
	//	relabelAirFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ScalarImageType3D::Pointer relabelAir = relabelAirFilter->GetOutput();

	//std::cout << "Number of objects: " << relabelAirFilter->GetNumberOfObjects() << std::endl;

	//WriteITK <ScalarImageType3D> ( relabelAir, "relabelAirThreshold.hdr");
	//
	//// Convert label image to label map
	//typedef unsigned long LabelType;	
	//typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
	//typedef itk::LabelMap< LabelObjectType > LabelMapType;

	//typedef itk::LabelImageToShapeLabelMapFilter< ScalarImageType3D, LabelMapType > 	ConverterType;
	//ConverterType::Pointer converter = ConverterType::New();
	//converter->SetInput( relabelAir );
	//converter->Update();
	//LabelMapType::Pointer labelMap = converter->GetOutput();

	//std::cout << "#\tPhysicalSize\tCentroid\tModifiedCentroid" << std::endl;

	//// Find labels for lungs and store
	//std::vector< unsigned int > lungLabel;

	//for( unsigned int label=2; label<=labelMap->GetNumberOfLabelObjects(); label++ )	//ignore labels 0 and 1 (body and external air)
	//{
	//	const LabelObjectType * lo = labelMap->GetLabelObject( label );
	//	
	//	float centroidX = lo->GetCentroid()[0] / spacing[0];
	//	float centroidY = lo->GetCentroid()[1] / spacing[1];
	//	float ps = lo->GetPhysicalSize();
	//	std::cout << label << "\t" << lo->GetPhysicalSize() << "\t" << lo->GetCentroid() << "\t[" << centroidX << ", " << centroidY << "]\t";

	//	if ( centroidY > size[1]/2 && ps > 100000) {
	//		std::cout << "FOUND LUNG";
	//		lungLabel.push_back( label );
	//	}

	//	std::cout << std::endl;
	//}

	//// Remove external air and lungs from air threshold
	//ScalarIteratorType relabelAirIt( relabelAir, region );

	//for (	relabelAirIt.GoToBegin(), airThresholdIt.GoToBegin();
	//		!relabelAirIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
	//		++relabelAirIt, ++airThresholdIt		) 
	//{
	//	if ( relabelAirIt.Get() == 1 )		{ airThresholdIt.Set( 0 ); }

	//	for (int i = 0; i < lungLabel.size() ; i++ ) {
	//		if ( relabelAirIt.Get() == lungLabel[i] ) {
	//			airThresholdIt.Set( 0 );
	//		}
	//	}
	//}

	//WriteITK <ScalarImageType3D> (airThreshold, "airThresholdWithoutLungs.hdr");


	//// Threshold to detect stool
	//ScalarImageType3D::Pointer stoolThreshold = ScalarImageType3D::New();
	//stoolThreshold->SetRegions( region );
	//stoolThreshold->Allocate();
	//ScalarIteratorType stoolThresholdIt( stoolThreshold, region );

	//for (	inputIt.GoToBegin(), stoolThresholdIt.GoToBegin();
	//		!inputIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
	//		++inputIt, ++stoolThresholdIt		) 
	//{
	//	if ( inputIt.Get() >= 200 )		{ stoolThresholdIt.Set( 1 ); }
	//	else							{ stoolThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> (stoolThreshold, "stoolThreshold.hdr");

	//// Run connected component filter to remove bone from stool threshold
	//ConnectedComponentFilterType::Pointer connectedComponentStoolFilter = ConnectedComponentFilterType::New();
	//connectedComponentStoolFilter->SetInput( stoolThreshold );

	//try
	//{
	//	connectedComponentStoolFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	////WriteITK <ScalarImageType3D> ( connectedComponentStoolFilter->GetOutput(), "connectedComponentStool.hdr");

	//// Relabel components
	//RelabelFilterType::Pointer relabelStoolFilter = RelabelFilterType::New();
	//relabelStoolFilter->SetInput( connectedComponentStoolFilter->GetOutput() );
	//relabelStoolFilter->SetMinimumObjectSize( 300000 );

	//try
	//{
	//	relabelStoolFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ScalarImageType3D::Pointer relabelStool = relabelStoolFilter->GetOutput();

	//std::cout << "Number of objects: " << relabelStoolFilter->GetNumberOfObjects() << std::endl;

	//WriteITK <ScalarImageType3D> ( relabelStool, "relabelStoolThreshold.hdr");

	//// Remove bone from stool threshold
	//ScalarIteratorType relabelStoolIt( relabelStool, region );
	//for (	relabelStoolIt.GoToBegin(), stoolThresholdIt.GoToBegin();
	//		!relabelStoolIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
	//		++relabelStoolIt, ++stoolThresholdIt		) 
	//{
	//	if ( relabelStoolIt.Get() == 1 )		{ stoolThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> ( stoolThreshold, "stoolThresholdWithoutBone.hdr");

	//// Remove noise from stool threshold image using median filter
	//BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	//
	//ScalarImageType3D::SizeType indexRadius;
	//indexRadius[0] = 1;
	//indexRadius[1] = 1;
	//indexRadius[2] = 1;
	//
	//medianFilter->SetRadius( indexRadius );
	//medianFilter->SetInput( stoolThreshold );
	//medianFilter->SetForegroundValue( 1 );
	//medianFilter->Update();

	//ScalarImageType3D::Pointer stoolThresholdSmooth = ScalarImageType3D::New();
	//stoolThresholdSmooth->SetRegions( region );
	//stoolThresholdSmooth->Allocate();
	//stoolThresholdSmooth = medianFilter->GetOutput();
	ScalarIteratorType stoolThresholdSmoothIt( stoolThresholdSmooth, region );

	//WriteITK <ScalarImageType3D> ( stoolThresholdSmooth, "stoolThresholdWithoutBoneSmoothed.hdr");

	/*// Add air stool boundaries to air image if stool is nearby
	ScalarImageType3D::Pointer airThresholdWithBoundaries = ScalarImageType3D::New();
	airThresholdWithBoundaries->SetRegions( region );
	airThresholdWithBoundaries->Allocate();
	ScalarIteratorType airThresholdWithBoundariesIt( airThresholdWithBoundaries, region );

	// Initialize air boundaries image to 0
	for (	airThresholdWithBoundariesIt.GoToBegin(); !airThresholdWithBoundariesIt.IsAtEnd();	++airThresholdWithBoundariesIt) 
	{
		airThresholdWithBoundariesIt.Set( 0 );
	}

	for (	stoolThresholdSmoothIt.GoToBegin(), airThresholdIt.GoToBegin();
			!stoolThresholdSmoothIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++stoolThresholdSmoothIt, ++airThresholdIt		) 
	{

			ScalarImageType3D::IndexType index = airThresholdIt.GetIndex();
	
			if ( airThresholdIt.Get() == 1 ) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										
										ScalarImageType3D::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]+k};
										
										if ( stoolThresholdSmooth->GetPixel(neighbor_index) == 1) {
											airThresholdWithBoundaries->SetPixel( neighbor_index, 1 );
										}
									}
								}
							}
						}
					}
				}
			}
	}

	WriteITK <ScalarImageType3D> ( airThresholdWithBoundaries, "airWithStoolBoundaries_1.hdr");
	*/

	// Median filter on air
	BinaryMedianFilterType::Pointer medianFilterAir = BinaryMedianFilterType::New();
	
	ScalarImageType3D::SizeType indexRadiusAir;
	indexRadiusAir[0] = 1;
	indexRadiusAir[1] = 1;
	indexRadiusAir[2] = 1;
	
	medianFilterAir->SetRadius( indexRadiusAir );
	medianFilterAir->SetInput( airThreshold );
	medianFilterAir->SetForegroundValue( 1 );
	medianFilterAir->Update();

	ScalarImageType3D::Pointer airThresholdSmooth = ScalarImageType3D::New();
	airThresholdSmooth->SetRegions( region );
	airThresholdSmooth->Allocate();
	airThresholdSmooth = medianFilterAir->GetOutput();
	ScalarIteratorType airThresholdSmoothIt( airThresholdSmooth, region );

	WriteITK <ScalarImageType3D> ( airThresholdSmooth, "airThresholdSmooth.hdr");

	// Dilate Air
	// Create binary ball structuring element
	StructuringElementType structuringElement;
    structuringElement.SetRadius( 1 );
    structuringElement.CreateStructuringElement();
	
	BinaryDilateFilterType::Pointer dilateAirFilter = BinaryDilateFilterType::New();
	dilateAirFilter->SetInput( airThresholdSmooth );
	dilateAirFilter->SetKernel( structuringElement );
	dilateAirFilter->SetDilateValue( 1 );
	dilateAirFilter->Update();

	ScalarImageType3D::Pointer dilateAir = ScalarImageType3D::New();
	dilateAir->SetRegions( region );
	dilateAir->Allocate();
	dilateAir = dilateAirFilter->GetOutput();

	ScalarIteratorType dilateAirIt( dilateAir, region );

	WriteITK <ScalarImageType3D> ( dilateAir, "dilateAir_1px.hdr");

	// Add stool regions to air image
	for (	stoolThresholdSmoothIt.GoToBegin(), dilateAirIt.GoToBegin();
			!stoolThresholdSmoothIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
			++stoolThresholdSmoothIt, ++dilateAirIt) 
	{
		if ( stoolThresholdSmoothIt.Get() == 1) { 
			dilateAirIt.Set( 1 ); 
		}
	}

	WriteITK <ScalarImageType3D> ( dilateAir, "taggedAll.hdr");

	/*// Add stool regions and air boundaries to air threshold image
	for (	stoolThresholdSmoothIt.GoToBegin(), airThresholdIt.GoToBegin(), airThresholdWithBoundariesIt.GoToBegin();
			!stoolThresholdSmoothIt.IsAtEnd() && !airThresholdIt.IsAtEnd() &&  !airThresholdWithBoundariesIt.IsAtEnd();
			++stoolThresholdSmoothIt, ++airThresholdIt, ++airThresholdWithBoundariesIt	) 
	{
		if ( stoolThresholdSmoothIt.Get() == 1 || airThresholdWithBoundariesIt.Get() == 1) { 
			airThresholdIt.Set( 1 ); 
		}
	}

	WriteITK <ScalarImageType3D> ( airThreshold, "taggedAllWithBoundaries_post.hdr");
	*/
	
	/*
	// Dilate air and stool using 3x3x3 region
	clock_t start = clock();

	ScalarImageType3D::Pointer dilateAirStool = ScalarImageType3D::New();
	dilateAirStool->SetRegions( region );
	dilateAirStool->SetSpacing( spacing );
	dilateAirStool->Allocate();
	ScalarIteratorType dilateAirStoolIt( dilateAirStool, region );

	// Copy image
	for (dilateAirStoolIt.GoToBegin(), airThresholdIt.GoToBegin();
		!dilateAirStoolIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
		++dilateAirStoolIt, ++airThresholdIt ) {
		
			dilateAirStoolIt.Set( airThresholdIt.Get() );
	}

	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	start = clock();

	// Dilate by 1 pixel in each dimension
	for (dilateAirStoolIt.GoToBegin(); !dilateAirStoolIt.IsAtEnd(); ++dilateAirStoolIt ) {

			ScalarImageType3D::IndexType index = dilateAirStoolIt.GetIndex();
	
			if ( dilateAirStoolIt.Get() == 1 ) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										
										ScalarImageType3D::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]+k};
										
										if ( dilateAirStool->GetPixel(neighbor_index) == 0) {
											dilateAirStool->SetPixel( neighbor_index, 1 );
										}
									}
								}
							}
						}
					}
				}
			}
	}

	WriteITK <ScalarImageType3D> ( dilateAirStool , "airThresholdDilated_1_new.hdr");
	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	*/

	// Run connected component to remove erroneus regions
	ConnectedComponentFilterType::Pointer connectedComponentAirStoolFilter = ConnectedComponentFilterType::New();
	connectedComponentAirStoolFilter->SetInput( dilateAir );

	try
	{
		connectedComponentAirStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	//WriteITK <ScalarImageType3D> ( connectedComponentAirStoolFilter->GetOutput(), "connectedComponentAirStool.hdr");

	// Relabel components
	RelabelFilterType::Pointer relabelAirStoolFilter = RelabelFilterType::New();
	relabelAirStoolFilter->SetInput( connectedComponentAirStoolFilter->GetOutput() );
	relabelAirStoolFilter->SetMinimumObjectSize( 5000 );

	try
	{
		relabelAirStoolFilter->Update();
	} catch ( itk::ExceptionObject &excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}

	ScalarImageType3D::Pointer relabelAirStool = relabelAirStoolFilter->GetOutput();

	std::cout << "Number of objects: " << relabelAirStoolFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelAirStool, "relabelAirStool.hdr");

	// Remove erroneous regions
	ScalarIteratorType relabelAirStoolIt( relabelAirStool, region );

	for (	relabelAirStoolIt.GoToBegin(), dilateAirIt.GoToBegin();
			!relabelAirStoolIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
			++relabelAirStoolIt, ++dilateAirIt		) 
	{
		if ( relabelAirStoolIt.Get() != 1 ) { dilateAirIt.Set( 0 ); }
	}

	WriteITK <ScalarImageType3D> ( dilateAir , "colonClean.hdr");

	clock_t start = clock();

	// Create binary ball structuring element
	StructuringElementType structuringElement2;
    structuringElement2.SetRadius( 4 );
    structuringElement2.CreateStructuringElement();

	// Dilate
	typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
	BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	dilateFilter->SetInput( dilateAir );
	dilateFilter->SetKernel( structuringElement2 );
	dilateFilter->SetDilateValue( 1 );
	dilateFilter->Update();

	WriteITK <ScalarImageType3D> ( dilateFilter->GetOutput() , "segmentedColon_4px.hdr");

	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);

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
