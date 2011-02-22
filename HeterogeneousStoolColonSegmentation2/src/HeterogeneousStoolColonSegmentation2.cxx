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
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
<<<<<<< HEAD
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkTextOutput.h"
 
typedef int														PixelType;
=======
 
typedef short													PixelType;
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1
typedef int														ScalarPixelType;
typedef itk::GDCMImageIO										ImageIOType;

typedef itk::Image< PixelType, 2 >								ImageType2D; 
typedef itk::Image< PixelType, 3 >								ImageType3D;
<<<<<<< HEAD
typedef itk::Image< float, 3>									FloatImageType3D;
=======
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

typedef itk::Image< ScalarPixelType, 2 >						ScalarImageType2D;
typedef itk::Image< ScalarPixelType, 3 >						ScalarImageType3D;	

typedef itk::ImageRegionIteratorWithIndex< ImageType3D >		IteratorType;
<<<<<<< HEAD
typedef itk::ImageRegionIteratorWithIndex< FloatImageType3D >	FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;
typedef itk::ImageRegionIterator< ScalarImageType3D >			ScalarIteratorTypeWithoutIndex;
=======
typedef itk::ImageRegionIteratorWithIndex< ScalarImageType3D >	ScalarIteratorType;
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

typedef itk::NeighborhoodIterator< ScalarImageType3D > NeighborhoodIteratorType;

typedef itk::ConnectedComponentImageFilter< ScalarImageType3D, ScalarImageType3D > ConnectedComponentFilterType;
typedef itk::RelabelComponentImageFilter< ScalarImageType3D, ScalarImageType3D > RelabelFilterType;

typedef itk::BinaryBallStructuringElement< ScalarImageType3D::PixelType, 3> StructuringElementType;
typedef itk::BinaryDilateImageFilter< ScalarImageType3D, ScalarImageType3D, StructuringElementType > BinaryDilateFilterType;
typedef itk::BinaryMedianImageFilter< ScalarImageType3D, ScalarImageType3D> BinaryMedianFilterType;

typedef itk::ImageSeriesReader< ImageType3D >        ImageSeriesReaderType;
typedef itk::RegularExpressionSeriesFileNames RegexFileNamesType;

typedef unsigned long LabelType;	
typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
typedef itk::LabelMap< LabelObjectType > LabelMapType;

typedef itk::LabelImageToShapeLabelMapFilter< ScalarImageType3D, LabelMapType > 	LabelImageToShapeLabelMapFilterType;
<<<<<<< HEAD
typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType3D, FloatImageType3D >  DiffusionFilterType;

struct point{
	int intensity;
	int size;
	int order;
};

typedef struct point ptype;
=======

>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

template< class T >
void WriteITK ( typename T::Pointer image , std::string prefix, std::string ss );
int writeCount = 1;
<<<<<<< HEAD
ImageType3D::Pointer ReadDicom( std::string path );
std::vector<std::string> explode( const std::string &delimiter, const std::string &str);
template <class ObjectType>
void printRefCount(ObjectType a, char* message);
void Relabel( ScalarImageType3D::Pointer cc, int minSize );
bool compare(ptype a, ptype b);
bool compareintensity(ptype a, ptype b);


int main( int argc, char* argv[] )
{
	itk::OutputWindow::SetInstance(itk::TextOutput::New());
	// Verify the number of parameters in the command line
	/*
=======

ImageType3D::Pointer ReadDicom( std::string path );

std::vector<std::string> explode( const std::string &delimiter, const std::string &str);

int main( int argc, char* argv[] )
{
	// Verify the number of parameters in the command line
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1
	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory\n";
		system("pause");
		return EXIT_FAILURE;
<<<<<<< HEAD
	}*/

	std::cout << "Starting the clock." << std::endl;
	clock_t start = clock();

	////std::string dataset = argv[1];

	//std::string dataset = "C:/GitProjects/HeterogeneousStoolColonSegmentation2/data/mr10-uncleansed-training/mr10_316_13s.i0355/dcm";

	//std::vector<std::string> datasetArr = explode( "/", dataset );
	//std::string ds = datasetArr[ datasetArr.size() - 2 ];
	//std::cout << "Dataset: " << ds << std::endl;

	//// Read and write dicom input
	//ImageType3D::Pointer input = ReadDicom( dataset );

	//// Get region info
	//ImageType3D::RegionType region = input->GetLargestPossibleRegion();
	//ImageType3D::IndexType endIndex = region.GetIndex();
	//ImageType3D::IndexType startIndex = region.GetIndex();	
	//endIndex[0]+=(region.GetSize()[0]-1);
	//endIndex[1]+=(region.GetSize()[1]-1);
	//endIndex[2]+=(region.GetSize()[2]-1);

	//ScalarImageType3D::SizeType size = region.GetSize();
	//ImageType3D::SpacingType spacing = input->GetSpacing();
	

	

	// Test images
	
	typedef itk::ImageFileReader< ImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();
	reader2->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/data/mr10_316_13s.i0355_1_input_sample.hdr" );
	//reader2->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/mr10_316_13s.i0355_5_stoolThreshold.hdr" );
	reader2->Update();
	ScalarImageType3D::Pointer input = reader2->GetOutput();
	ScalarImageType3D::RegionType region = input->GetLargestPossibleRegion();
	std::string ds = "";
	
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );

	WriteITK <ImageType3D> (input, ds, "input.hdr");
	

	//// Threshold to detect air
	//std::cout << "Thresholding air" << std::endl;
	//ScalarImageType3D::Pointer airThreshold = ScalarImageType3D::New();
	//airThreshold->SetRegions( region );
	//airThreshold->SetSpacing( spacing );
	//airThreshold->Allocate();
	//ScalarIteratorType airThresholdIt( airThreshold, region );

	//for (	inputIt.GoToBegin(), airThresholdIt.GoToBegin();
	//		!inputIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
	//		++inputIt, ++airThresholdIt		) 
	//{
	//	if ( inputIt.Get() <= -600 )	{ airThresholdIt.Set( 1 ); }
	//	else							{ airThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> (airThreshold, ds, "airThreshold.hdr");

	//// Run connected component on air to detect lungs
	//std::cout << "Detecting lungs" << std::endl;
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	//ccFilter->SetInput( airThreshold );
	//ccFilter->Update();

	//ScalarImageType3D::Pointer relabelAir = ccFilter->GetOutput();
	//
	////WriteITK <ScalarImageType3D> ( relabelAir, ds, "connectedComponentAir.hdr");

	////// Relabel components
	////RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	////relabelFilter->SetInput( ccFilter->GetOutput() );

	////// Set the minimum size to detect lungs
	////relabelFilter->SetMinimumObjectSize( 20000 );
	//std::cout << "Relabeling air" << std::endl;
	//Relabel( relabelAir, 20000);
	//WriteITK <ScalarImageType3D> ( relabelAir, ds, "relabelAirThreshold.hdr");
	//
	//// Convert label image to label map
	//LabelImageToShapeLabelMapFilterType::Pointer convertLabelMapFilter = LabelImageToShapeLabelMapFilterType::New();
	//convertLabelMapFilter->SetInput( relabelAir );
	//convertLabelMapFilter->Update();
	//LabelMapType::Pointer labelMapAir = convertLabelMapFilter->GetOutput();

	//std::cout << "Segmenting lung components" << std::endl;

	//// Output lung data to text file
	//std::ofstream myfile;
	//std::string textFileName = ds + "_data.txt";
	//myfile.open( textFileName.c_str() );
	//myfile << "Dataset: " << ds << "\n";

	//// Cutoff percentages to detect centroid location of lungs
	//float cutoffY = 0.6;
	//float cutoffZ = 0.15;

	//myfile << "Size in pixel units: [" << size[0] << ", " << size[1] << ", " << size[2] << "]\n";
	//myfile << "Size in physical units: [" << size[0]*spacing[0] << ", " << size[1]*spacing[1] << ", " << size[2]*spacing[2] << "]\n";
	//myfile << "Y-axis lung cutoff in pixel units: " << cutoffY*size[1] << "\n";
	//myfile << "Z-axis lung cutoff in pixel units: " << cutoffZ*size[2] << "\n\n\n";

	//myfile << "#\tPhysicalSize\tCentroid\tModifiedCentroid\n";

	//// Find labels for lungs and store
	//std::vector< unsigned int > lungLabel;

	//for( unsigned int label=1; label<=labelMapAir->GetNumberOfLabelObjects(); label++ )	//ignore labels 0 and 1 (body and external air)
	//{
	//	const LabelObjectType * lo = labelMapAir->GetLabelObject( label );

	//	// Use pixel units, view in MIPAV
	//	float centroidX = lo->GetCentroid()[0] / spacing[0];
	//	float centroidY = size[1] - (lo->GetCentroid()[1] / spacing[1]); // flipped in MIPAV
	//	float centroidZ = lo->GetCentroid()[2] / spacing[2];

	//	float ps = lo->GetPhysicalSize();
	//	myfile << label << "\t" << lo->GetPhysicalSize() << "\t" << lo->GetCentroid() << "\t[" << centroidX << ", " << centroidY << ", " << centroidZ << "]\t";

	//	// Set physical size threshold that is proportional to the z centroid
	//	float propPS = 68000*100*centroidZ/size[2] - 186000;

	//	if ( centroidY < cutoffY*size[1] && centroidZ < cutoffZ*size[2] && ( lo->GetPhysicalSize() > propPS ) ) {
	//		myfile << "FOUND LUNG";
	//		lungLabel.push_back( label );
	//	}

	//	myfile << "\n";
	//}

	//myfile.close();

	//std::cout << "Removing lungs" << std::endl;
=======
	}

	// Start clock
	clock_t start = clock();

	std::string dataset = argv[1];

	//std::string dataset = "C:/GitProjects/ReadWriteDicom/data/mr10-uncleansed/mr10_316_13s.i0355/dcm";
	//std::string dataset = "C:/GitProjects/ReadWriteDicom/data/mr10-uncleansed/mr10_342_13p.i0396/dcm";
	//std::string dataset = "C:/GitProjects/ReadWriteDicom/data/mr10-uncleansed/mr10_358_13s.i0340/dcm";
	//std::string dataset = "C:/GitProjects/ReadWriteDicom/data/mr10-uncleansed/mr10_351_13p.i0448/dcm";
	//std::string dataset = "C:/GitProjects/ReadWriteDicom/data/mr10-uncleansed/mr10_379_08s.i0589/dcm";
	

	std::vector<std::string> datasetArr = explode( "/", dataset );
	std::string ds = datasetArr[ datasetArr.size() - 2 ];

	std::cout << "Dataset: " << ds << std::endl;

	// Read and write dicom input
	ImageType3D::Pointer input = ReadDicom( dataset );

	// Read input
	//typedef itk::ImageFileReader< ImageType3D > ReaderType;
	//ReaderType::Pointer reader = ReaderType::New();

	////reader->SetFileName( "C:/ImageData/mr10_092_13p.i0344.hdr" );
	//reader->SetFileName( "C:/ImageData/mr10_047_08s_i0526_1-150.hdr" );
	//reader->Update();
	//ImageType3D::Pointer input = reader->GetOutput();

	ImageType3D::RegionType region = input->GetLargestPossibleRegion();

	ImageType3D::IndexType endIndex = region.GetIndex();
	ImageType3D::IndexType startIndex = region.GetIndex();	
	endIndex[0]+=(region.GetSize()[0]-1);
	endIndex[1]+=(region.GetSize()[1]-1);
	endIndex[2]+=(region.GetSize()[2]-1);

	ScalarImageType3D::SizeType size = region.GetSize();
	ImageType3D::SpacingType spacing = input->GetSpacing();
	IteratorType inputIt( input, input->GetLargestPossibleRegion() );

	WriteITK <ImageType3D> (input, ds, "input.hdr");

	/*// Test images
	typedef itk::ImageFileReader< ScalarImageType3D > ReaderType2;
	ReaderType2::Pointer reader2 = ReaderType2::New();
	//reader2->SetFileName( "trial1data/4_airThresholdWithoutLungs.hdr" );
	reader2->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/3_relabelAirThreshold.hdr" );
	reader2->Update();
	ScalarImageType3D::Pointer relabelAir = reader2->GetOutput();*/
	

	/*
	ReaderType2::Pointer reader3 = ReaderType2::New();
	reader3->SetFileName( "C:/GitProjects/HeterogeneousStoolColonSegmentation2/build64-3.16.0/4_airThresholdWithoutLungs.hdr" );
	reader3->Update();
	ScalarImageType3D::Pointer airThreshold = reader3->GetOutput();

	ReaderType2::Pointer reader3 = ReaderType2::New();
	reader3->SetFileName( "trial1data/7_stoolThresholdWithoutBone.hdr" );
	reader3->Update();
	ScalarImageType3D::Pointer stoolThreshold = reader3->GetOutput();
	stoolThreshold->SetSpacing( spacing );
	*/
	
	// Threshold to detect non-air and isolate colon body
	ScalarImageType3D::Pointer airThreshold = ScalarImageType3D::New();
	airThreshold->SetRegions( region );
	airThreshold->SetSpacing( spacing );
	airThreshold->Allocate();
	ScalarIteratorType airThresholdIt( airThreshold, region );

	for (	inputIt.GoToBegin(), airThresholdIt.GoToBegin();
			!inputIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++inputIt, ++airThresholdIt		) 
	{
		if ( inputIt.Get() <= -600 )	{ airThresholdIt.Set( 0 ); }
		else							{ airThresholdIt.Set( 1 ); }
	}

	WriteITK <ScalarImageType3D> (airThreshold, ds, "airThreshold.hdr");

	// Run connected component filter to detect colon body
	ConnectedComponentFilterType::Pointer connectedComponentNonAirFilter = ConnectedComponentFilterType::New();
	connectedComponentNonAirFilter->SetInput( airThreshold );
	connectedComponentNonAirFilter->Update();
	//WriteITK <ScalarImageType3D> ( connectedComponentNonAirFilter->GetOutput(), ds, "connectedComponentNonAir.hdr");

	// Relabel components
	RelabelFilterType::Pointer relabelNonAirFilter = RelabelFilterType::New();
	relabelNonAirFilter->SetInput( connectedComponentNonAirFilter->GetOutput() );

	// Set the minimum size as 20% of all pixels to detect colon body
	relabelNonAirFilter->SetMinimumObjectSize( 0.2*size[0]*size[1]*size[2] );
	relabelNonAirFilter->Update();

	ScalarImageType3D::Pointer relabelNonAir = relabelNonAirFilter->GetOutput();
	ScalarIteratorType relabelNonAirIt( relabelNonAir, region );
	std::cout << "Number of nonAir objects: " << relabelNonAirFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelNonAir, ds, "relabelNonAirThreshold.hdr");
	
	// Convert label image to label map
	LabelImageToShapeLabelMapFilterType::Pointer converter = LabelImageToShapeLabelMapFilterType::New();
	converter->SetInput( relabelNonAir );
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();

	// Get region of body, label 1
	const LabelObjectType * lo = labelMap->GetLabelObject( 1 );
	ScalarImageType3D::RegionType bodyRegion = lo->GetRegion();

	// Reset air threshold to segment only body
	for (	relabelNonAirIt.GoToBegin(), airThresholdIt.GoToBegin();
			!relabelNonAirIt.IsAtEnd() && !airThresholdIt.IsAtEnd();
			++relabelNonAirIt, ++airThresholdIt		) 
	{
		if ( airThresholdIt.Get() == 0 ) {
			airThresholdIt.Set( 1 );
		} else {
			airThresholdIt.Set( 0 );
		}

		if ( relabelNonAirIt.Get() != 1 ) { // convert outer non-air to air
			airThresholdIt.Set( 1 );
		}
	}

	WriteITK <ScalarImageType3D> ( airThreshold, ds, "airThresholdBody.hdr");

	// Run connected component on air to detect lungs
	ConnectedComponentFilterType::Pointer connectedComponentAirFilter = ConnectedComponentFilterType::New();
	connectedComponentAirFilter->SetInput( airThreshold );
	connectedComponentAirFilter->Update();
	//WriteITK <ScalarImageType3D> ( connectedComponentAirFilter->GetOutput(), ds, "connectedComponentAir.hdr");

	// Relabel components
	RelabelFilterType::Pointer relabelAirFilter = RelabelFilterType::New();
	relabelAirFilter->SetInput( connectedComponentAirFilter->GetOutput() );

	// Set the minimum size to detect lungs
	relabelAirFilter->SetMinimumObjectSize( 5000 );
	relabelAirFilter->Update();

	ScalarImageType3D::Pointer relabelAir = relabelAirFilter->GetOutput();
	std::cout << "Number of air objects: " << relabelAirFilter->GetNumberOfObjects() << std::endl;

	WriteITK <ScalarImageType3D> ( relabelAir, ds, "relabelAirThreshold.hdr");
	
	// Convert label image to label map
	LabelImageToShapeLabelMapFilterType::Pointer converterAir = LabelImageToShapeLabelMapFilterType::New();
	converterAir->SetInput( relabelAir );
	converterAir->Update();
	LabelMapType::Pointer labelMapAir = converterAir->GetOutput();

	// Output lung data to text file
	std::ofstream myfile;
	std::string textFileName = ds + "_data.txt";
	myfile.open( textFileName.c_str() );
	myfile << "Dataset: " << ds << "\n";

	// Cutoff percentages to detect centroid location of lungs
	float cutoffY = 0.6;
	float cutoffZ = 0.15;

	myfile << "Size in pixel units: [" << size[0] << ", " << size[1] << ", " << size[2] << "]\n";
	myfile << "Size in physical units: [" << size[0]*spacing[0] << ", " << size[1]*spacing[1] << ", " << size[2]*spacing[2] << "]\n";
	myfile << "Y-axis lung cutoff in pixel units: " << cutoffY*size[1] << "\n";
	myfile << "Z-axis lung cutoff in pixel units: " << cutoffZ*size[2] << "\n\n\n";

	myfile << "#\tPhysicalSize\tCentroid\tModifiedCentroid\tInsideBody\n";

	// Find labels for lungs and store
	std::vector< unsigned int > lungLabel;

	for( unsigned int label=1; label<=labelMapAir->GetNumberOfLabelObjects(); label++ )	//ignore labels 0 and 1 (body and external air)
	{
		const LabelObjectType * lo = labelMapAir->GetLabelObject( label );

		// Use pixel units, view in MIPAV
		float centroidX = lo->GetCentroid()[0] / spacing[0];
		float centroidY = size[1] - (lo->GetCentroid()[1] / spacing[1]); // flipped in MIPAV
		float centroidZ = lo->GetCentroid()[2] / spacing[2];

		float ps = lo->GetPhysicalSize();
		myfile << label << "\t" << lo->GetPhysicalSize() << "\t" << lo->GetCentroid() << "\t[" << centroidX << ", " << centroidY << ", " << centroidZ << "]\t" << bodyRegion.IsInside( lo->GetRegion() ) << "\t";

		// Set physical size threshold that is proportional to the z centroid
		float propPS = 68000*100*centroidZ/size[2] - 186000;

		if ( centroidY < cutoffY*size[1] && centroidZ < cutoffZ*size[2] && bodyRegion.IsInside( lo->GetRegion() ) && ( lo->GetPhysicalSize() > propPS ) ) {
			myfile << "FOUND LUNG";
			lungLabel.push_back( label );
		}

		myfile << "\n";
	}

	myfile.close();
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

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

	//WriteITK <ScalarImageType3D> (airThreshold, ds, "airThresholdWithoutLungs.hdr");

<<<<<<< HEAD
	//std::cout << "Thresholding stool" << std::endl;

	// input -> diffusion -> threshold -> connected component -> relabel
	float conductance[1] = {0.5};
	int numIter = 5;
	std::stringstream ss;

	for (int th=420; th < 500; th += 20) {
		ScalarImageType3D::Pointer stoolThreshold = ScalarImageType3D::New();
		stoolThreshold->SetRegions( region );
		stoolThreshold->Allocate();
		ScalarIteratorType stoolThresholdIt( stoolThreshold, region );
		
		for (	inputIt.GoToBegin(), stoolThresholdIt.GoToBegin();
					!inputIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
					++inputIt, ++stoolThresholdIt		) 
		{
			if ( inputIt.Get() >= th )		{ stoolThresholdIt.Set( 1 ); }
			else							{ stoolThresholdIt.Set( 0 ); }
		}

		ss.str("");
		ss << "th_" << th << "_inputThreshold.hdr";
		WriteITK <ScalarImageType3D> (stoolThreshold, ds, ss.str());

		//for (int i=0; i<1; i++) {
			/*
			// Perform anisotropic diffusion
			DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
			diffusionFilter->SetInput( input );
			diffusionFilter->SetNumberOfIterations( numIter );
			diffusionFilter->SetTimeStep( 0.25 );
			diffusionFilter->SetConductanceParameter( conductance[i] );
			diffusionFilter->Update();
			FloatImageType3D::Pointer diffuse = diffusionFilter->GetOutput();
			FloatIteratorType diffuseIt ( diffuse, region );
			ss.str("");
			ss << "th_" << th << "_conductance_" << conductance[i] << "_iter_" << numIter << "_diffusion.hdr";
			WriteITK < FloatImageType3D > (diffuse, ds, ss.str() );
			*/

		/*
			// Threshold to detect stool
			for (	inputIt.GoToBegin(), stoolThresholdIt.GoToBegin();
					!inputIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
					++inputIt, ++stoolThresholdIt		) 
			{
				if ( inputIt.Get() >= th )	{ stoolThresholdIt.Set( 1 ); }
				else							{ stoolThresholdIt.Set( 0 ); }
			}

			ss.str("");
			ss << "th_" << th << "_conductance_" << conductance[i] << "_iter_" << numIter <<  "_inputThreshold.hdr";
			WriteITK <ScalarImageType3D> (stoolThreshold, ds, ss.str());
			*/

			//Reuse pipeline filters
			//ccFilter.~SmartPointer();

			// Run connected component filter to remove bone from stool threshold
			std::cout << "Detecting bone" << std::endl;
			ccFilter = ConnectedComponentFilterType::New();
			ccFilter->SetInput( stoolThreshold );
			ccFilter->Update();
			ScalarImageType3D::Pointer relabelStool = ccFilter->GetOutput();

			//WriteITK <ScalarImageType3D> ( relabelStool, ds, "connectedComponentStool.hdr");

			// Relabel components
			/*RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
			relabelFilter->SetInput( relabelStool );
			relabelFilter->SetMinimumObjectSize( 300000 );
			relabelFilter->Update();*/
			Relabel( relabelStool, 0 );

			// Isolate top 10 components
			ScalarIteratorType rsIt( relabelStool, relabelStool->GetLargestPossibleRegion() );
			for ( rsIt.GoToBegin(); !rsIt.IsAtEnd(); ++rsIt )
			{
				if ( rsIt.Get() > 10 ) { rsIt.Set( 0 ); }
			}
			
			ss.str("");
			ss << "th_" << th << "_relabel.hdr";
			WriteITK <ScalarImageType3D> ( relabelStool, ds, ss.str());

			//// Remove bone from stool threshold
			//std::cout << "Removing bone" << std::endl;
			//ScalarIteratorType relabelStoolIt( relabelStool, region );
			//for (	relabelStoolIt.GoToBegin(), stoolThresholdIt.GoToBegin();
			//		!relabelStoolIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
			//		++relabelStoolIt, ++stoolThresholdIt		) 
			//{
			//	if ( relabelStoolIt.Get() == 1 )		{ stoolThresholdIt.Set( 0 ); }
			//}

			//ss.str("");
			//ss << "stoolThresholdWithoutBone_" << th << ".hdr";
			//WriteITK <ScalarImageType3D> ( stoolThreshold, ds, ss.str());

		//}
	}

	//// Reuse median filter
	//medianFilter.~SmartPointer();

	//// Median filter on air
	//medianFilter = BinaryMedianFilterType::New();
	//
	//medianFilter->SetRadius( indexRadius );
	//medianFilter->SetInput( airThreshold );
	//medianFilter->SetForegroundValue( 1 );
	//medianFilter->Update();

	//ScalarImageType3D::Pointer airThresholdSmooth = ScalarImageType3D::New();
	//airThresholdSmooth = medianFilter->GetOutput();
=======

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

	//WriteITK <ScalarImageType3D> (stoolThreshold, ds, "stoolThreshold.hdr");

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

	////WriteITK <ScalarImageType3D> ( connectedComponentStoolFilter->GetOutput(), ds, "connectedComponentStool.hdr");

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

	//WriteITK <ScalarImageType3D> ( relabelStool, ds, "relabelStoolThreshold.hdr");

	//// Remove bone from stool threshold
	//ScalarIteratorType relabelStoolIt( relabelStool, region );
	//for (	relabelStoolIt.GoToBegin(), stoolThresholdIt.GoToBegin();
	//		!relabelStoolIt.IsAtEnd() && !stoolThresholdIt.IsAtEnd();
	//		++relabelStoolIt, ++stoolThresholdIt		) 
	//{
	//	if ( relabelStoolIt.Get() == 1 )		{ stoolThresholdIt.Set( 0 ); }
	//}

	//WriteITK <ScalarImageType3D> ( stoolThreshold, ds, "stoolThresholdWithoutBone.hdr");

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
	//ScalarIteratorType stoolThresholdSmoothIt( stoolThresholdSmooth, region );

	//WriteITK <ScalarImageType3D> ( stoolThresholdSmooth, ds, "stoolThresholdWithoutBoneSmoothed.hdr");

	//// Median filter on air
	//BinaryMedianFilterType::Pointer medianFilterAir = BinaryMedianFilterType::New();
	//
	//ScalarImageType3D::SizeType indexRadiusAir;
	//indexRadiusAir[0] = 1;
	//indexRadiusAir[1] = 1;
	//indexRadiusAir[2] = 1;
	//
	//medianFilterAir->SetRadius( indexRadiusAir );
	//medianFilterAir->SetInput( airThreshold );
	//medianFilterAir->SetForegroundValue( 1 );
	//medianFilterAir->Update();

	//ScalarImageType3D::Pointer airThresholdSmooth = ScalarImageType3D::New();
	//airThresholdSmooth->SetRegions( region );
	//airThresholdSmooth->Allocate();
	//airThresholdSmooth = medianFilterAir->GetOutput();
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1
	//ScalarIteratorType airThresholdSmoothIt( airThresholdSmooth, region );

	//WriteITK <ScalarImageType3D> ( airThresholdSmooth, ds, "airThresholdSmooth.hdr");

	//// Dilate Air
	//// Create binary ball structuring element
	//StructuringElementType structuringElement;
 //   structuringElement.SetRadius( 1 );
 //   structuringElement.CreateStructuringElement();
	//
	//BinaryDilateFilterType::Pointer dilateAirFilter = BinaryDilateFilterType::New();
	//dilateAirFilter->SetInput( airThresholdSmooth );
	//dilateAirFilter->SetKernel( structuringElement );
	//dilateAirFilter->SetDilateValue( 1 );
	//dilateAirFilter->Update();

	//ScalarImageType3D::Pointer dilateAir = ScalarImageType3D::New();
	//dilateAir->SetRegions( region );
	//dilateAir->Allocate();
	//dilateAir = dilateAirFilter->GetOutput();

	//ScalarIteratorType dilateAirIt( dilateAir, region );

	//WriteITK <ScalarImageType3D> ( dilateAir, ds, "dilateAir_1px.hdr");

	//// Add stool regions to air image
	//for (	stoolThresholdSmoothIt.GoToBegin(), dilateAirIt.GoToBegin();
	//		!stoolThresholdSmoothIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
	//		++stoolThresholdSmoothIt, ++dilateAirIt) 
	//{
	//	if ( stoolThresholdSmoothIt.Get() == 1) { 
	//		dilateAirIt.Set( 1 ); 
	//	}
	//}

	//WriteITK <ScalarImageType3D> ( dilateAir, ds, "taggedAll.hdr");

	//// Run connected component to remove erroneus regions
<<<<<<< HEAD
	//ccFilter.~SmartPointer();
	//ccFilter = ConnectedComponentFilterType::New();
	//ccFilter->SetInput( dilateAir );
	//ccFilter->Update();
	//ScalarImageType3D::Pointer relabelAirStool = ccFilter->GetOutput();

	////WriteITK <ScalarImageType3D> ( ccFilter->GetOutput(), ds, "connectedComponentAirStool.hdr");

	//// Relabel components
	///*RelabelFilterType::Pointer relabelAirStoolFilter = RelabelFilterType::New();
	//relabelAirStoolFilter->SetInput( connectedComponentAirStoolFilter->GetOutput() );
	//relabelAirStoolFilter->SetMinimumObjectSize( 5000 );
	//relabelAirStoolFilter->Update();*/
	//Relabel( relabelAirStool, 5000 );
=======
	//ConnectedComponentFilterType::Pointer connectedComponentAirStoolFilter = ConnectedComponentFilterType::New();
	//connectedComponentAirStoolFilter->SetInput( dilateAir );

	//try
	//{
	//	connectedComponentAirStoolFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	////WriteITK <ScalarImageType3D> ( connectedComponentAirStoolFilter->GetOutput(), ds, "connectedComponentAirStool.hdr");

	//// Relabel components
	//RelabelFilterType::Pointer relabelAirStoolFilter = RelabelFilterType::New();
	//relabelAirStoolFilter->SetInput( connectedComponentAirStoolFilter->GetOutput() );
	//relabelAirStoolFilter->SetMinimumObjectSize( 5000 );

	//try
	//{
	//	relabelAirStoolFilter->Update();
	//} catch ( itk::ExceptionObject &excep )
	//{
	//	std::cerr << "Exception caught !" << std::endl;
	//	std::cerr << excep << std::endl;
	//}

	//ScalarImageType3D::Pointer relabelAirStool = relabelAirStoolFilter->GetOutput();

	//std::cout << "Number of objects: " << relabelAirStoolFilter->GetNumberOfObjects() << std::endl;
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

	//WriteITK <ScalarImageType3D> ( relabelAirStool, ds, "relabelAirStool.hdr");

	//// Remove erroneous regions
	//ScalarIteratorType relabelAirStoolIt( relabelAirStool, region );

	//for (	relabelAirStoolIt.GoToBegin(), dilateAirIt.GoToBegin();
	//		!relabelAirStoolIt.IsAtEnd() && !dilateAirIt.IsAtEnd();
	//		++relabelAirStoolIt, ++dilateAirIt		) 
	//{
	//	if ( relabelAirStoolIt.Get() != 1 ) { dilateAirIt.Set( 0 ); }
	//}

<<<<<<< HEAD
	////WriteITK <ScalarImageType3D> ( dilateAir , ds, "colonClean.hdr");
=======
	//WriteITK <ScalarImageType3D> ( dilateAir , ds, "colonClean.hdr");
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

	//// Create binary ball structuring element
	//StructuringElementType structuringElement2;
 //   structuringElement2.SetRadius( 10 );
 //   structuringElement2.CreateStructuringElement();

	//// Dilate
	//BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
	//dilateFilter->SetInput( dilateAir );
	//dilateFilter->SetKernel( structuringElement2 );
	//dilateFilter->SetDilateValue( 1 );
	//dilateFilter->Update();

	//WriteITK <ScalarImageType3D> ( dilateFilter->GetOutput() , ds, "segmentedColon_10px.hdr");

	printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
<<<<<<< HEAD

	system("pause");
=======
	//system("pause");
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1
	return 0;
}

template< class T >
void WriteITK (typename T::Pointer image , std::string prefix, std::string ss )
{
	typedef itk::ImageFileWriter< T >	WriterType;
	WriterType::Pointer writer = WriterType::New();

	std::stringstream ss2;
<<<<<<< HEAD
	//ss2 << prefix << "_" << writeCount++ << "_" << ss;
	ss2 << ss;
=======
	ss2 << prefix << "_" << writeCount++ << "_" << ss;
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

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

ImageType3D::Pointer ReadDicom( std::string path )
{	
	// Create reader
	ImageSeriesReaderType::Pointer reader = ImageSeriesReaderType::New();
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	RegexFileNamesType::Pointer fit = RegexFileNamesType::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
<<<<<<< HEAD
	//names.erase( names.begin()+130, names.end() );
=======
	//names.erase( names.begin()+150, names.end() );
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1

    reader->SetFileNames( names );
    reader->Update();

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType3D,ImageType3D>::Pointer orienter =     itk::OrientImageFilter<ImageType3D,ImageType3D>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI);
    orienter->SetInput(reader->GetOutput());
    orienter->Update();

    return orienter->GetOutput();

}

std::vector<std::string> explode( const std::string &delimiter, const std::string &str)
{
	std::vector<std::string> arr;

    int strleng = str.length();
    int delleng = delimiter.length();
    if (delleng==0)
        return arr;//no change

    int i=0; 
    int k=0;
    while( i<strleng )
    {
        int j=0;
        while (i+j<strleng && j<delleng && str[i+j]==delimiter[j])
            j++;
        if (j==delleng)//found delimiter
        {
            arr.push_back(  str.substr(k, i-k) );
            i+=delleng;
            k=i;
        }
        else
        {
            i++;
        }
    }
    arr.push_back(  str.substr(k, i-k) );
    return arr;
<<<<<<< HEAD
}

template <class ObjectType>
void printRefCount(ObjectType a, char* message)
{
	std::cout << message << std::endl;
	std::cout << "Reference Count is: " << a->GetReferenceCount() - 1 << std::endl;
}

void Relabel( ScalarImageType3D::Pointer cc, int minSize ) {
	ScalarIteratorTypeWithoutIndex it(cc, cc->GetLargestPossibleRegion());
	int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		count.push_back(a);
	}
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			++count[it.Get() - 1].size;
		}
	}
	
	sort(count.begin(), count.end(), compare);

	for(int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareintensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			if ( count[it.Get() - 1].size > minSize ) {
				it.Set(count[it.Get() - 1].order);
			} else {
				it.Set( 0 );
			}
		}
	}
}

bool compare(ptype a, ptype b)
{
	if(a.size > b.size)
	{
		return true;
	}
	else return false;
}

bool compareintensity(ptype a, ptype b)
{
	if(a.intensity < b.intensity)
	{
		return true;
	}
	else return false;
=======
>>>>>>> 7fd8874934cff62e22e9c5b170ab4f0e900d5ad1
}