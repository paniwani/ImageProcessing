/*
Optimize Voxel Edge Classification by fixing Air to Tissue transitions
which are shown as Tissue-Air to Stool-Tissue
*/

#include "itkImage.h"
#include <iostream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include "itkChamferDistanceTransformImageFilter.h"
#include <vector>
#include <fstream>
#include "pcre.h"
#include <itkGradientImageFilter.h>
#include <itkContinuousIndex.h>
#include <time.h>
#include "itkDiscreteGaussianImageFilter.h"

typedef itk::Image< float, 3 > ImageType;
typedef itk::BSplineInterpolateImageFunction<ImageType, float> InterpolationType;
typedef itk::CovariantVector< float, 3> CovariantVectorType;
typedef itk::Image< CovariantVectorType, 3 > CovariantImageType;
typedef itk::ContinuousIndex< float, 3 > ContinuousIndexType;

void FindVoxelsByGradient(ImageType::Pointer voxelEdge, ImageType::IndexType index, ImageType::IndexType startIndex, ImageType::IndexType endIndex, CovariantVectorType grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector);
void FindVoxelsByGradient2(ImageType::Pointer voxelEdge, ImageType::IndexType index, ImageType::IndexType startIndex, ImageType::IndexType endIndex, CovariantVectorType grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector);

void ReadITK(ImageType::Pointer &image, char * fileName);
void WriteITK(ImageType::Pointer image, std::string ss, int count);
double round(float d);
bool MatchVoxels(std::string inputString, char * regex);

const double PI = 3.141592;

int main(int argc, char * argv[])
{

	if( argc < 3 ) 
    { 
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "  inputImage voxelEdgeFile  " << std::endl;
		return EXIT_FAILURE;
    }

	//Start the clock
	time_t tic;
	tic = time(NULL);
	std::cout << "Clock started." << std::endl;

	// PARAMETERS
	int numOfVoxels = 10; // to move in direction of gradient

	char * inFileName = argv[1];
	char * veFileName = argv[2];

	// IO vars
	int count = 1;
	std::stringstream ss;
	std::string suffix = "";
	double dtheta = PI/3;

	// Pattern matching
	char *regex = "^2*[5-7]+1+$";

	// Setup text file
	std::ofstream myfile;
	ss.str("");
	ss << "output" << suffix << ".txt";
	myfile.open (ss.str().c_str());
   
	// Read input image
	ImageType::Pointer inputFull = ImageType::New();
	ReadITK(inputFull, inFileName);
	myfile << "Input: " << inFileName << "\n";

	// Read voxel image
	ImageType::Pointer voxelEdgeFull = ImageType::New();
	ReadITK(voxelEdgeFull, veFileName);

	myfile << "VoxelEdge Input: " << veFileName << "\n\n";
	myfile << "Number of voxels to move in direction of gradient: " << numOfVoxels << "\n";
	myfile << "Pattern: " << regex << "\n\n";

	// Get subregions
	ImageType::RegionType region = inputFull->GetLargestPossibleRegion();
	//ImageType::SizeType offset_size = {80, 40, 1};
	//ImageType::IndexType offset_index = {210,512-360,0};
	ImageType::SpacingType spacing = inputFull->GetSpacing();
	//ImageType::RegionType region;
	//region.SetSize(offset_size);
	//region.SetIndex(offset_index);

	ImageType::IndexType endIndex = region.GetIndex();
	ImageType::IndexType startIndex = region.GetIndex();	
	endIndex[0]+=(region.GetSize()[0]-1);
	endIndex[1]+=(region.GetSize()[1]-1);
	endIndex[2]+=(region.GetSize()[2]-1);

	inputFull->SetRequestedRegion( region );
	inputFull->Update();

	ImageType::Pointer input = ImageType::New();
	input->SetRegions( region );
	input->Allocate();
	input->SetSpacing( spacing );

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType inputFullIt( inputFull, region );
	IteratorType inputIt ( input, region );
	
	for(inputFullIt.GoToBegin() , inputIt.GoToBegin() ; !inputFullIt.IsAtEnd() && !inputIt.IsAtEnd(); ++inputFullIt, ++inputIt) {
		inputIt.Set(inputFullIt.Get());
	}


	ImageType::Pointer voxelEdge = ImageType::New();
	voxelEdge->SetRegions( region );
	voxelEdge->Allocate();
	voxelEdge->SetSpacing( spacing );

	ImageType::Pointer voxelEdgeUpdate = ImageType::New();
	voxelEdgeUpdate->SetRegions( region );
	voxelEdgeUpdate->Allocate();
	voxelEdgeUpdate->SetSpacing( spacing );

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType voxelEdgeFullIt( voxelEdgeFull, region );
	IteratorType voxelEdgeIt ( voxelEdge, region );
	IteratorType voxelEdgeUpdateIt ( voxelEdgeUpdate, region );

	
	for(voxelEdgeFullIt.GoToBegin() , voxelEdgeIt.GoToBegin(), voxelEdgeUpdateIt.GoToBegin();
		!voxelEdgeFullIt.IsAtEnd() && !voxelEdgeIt.IsAtEnd() && !voxelEdgeUpdateIt.IsAtEnd(); 
		++voxelEdgeFullIt, ++voxelEdgeIt, ++voxelEdgeUpdateIt) 
	{
		voxelEdgeIt.Set(voxelEdgeFullIt.Get());
		voxelEdgeUpdateIt.Set(voxelEdgeFullIt.Get());
	}

	// Write new input region
	ss.str("");
	ss << "input_subregion.hdr";
	WriteITK(input, ss.str(), count++);	

	// Write new voxel edge region
	ss.str("");
	ss << "voxelEdge_subregion.hdr";
	WriteITK(voxelEdge, ss.str(), count++);	

	// Calculate gradient
	typedef itk::GradientImageFilter< ImageType > GradientFilterType;
	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput(input);

	try  
	{
		gradientFilter->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	}

	CovariantImageType::Pointer gradient = gradientFilter->GetOutput();

	// Output gradient vector information to text file
	typedef itk::ImageRegionIteratorWithIndex< CovariantImageType > CovariantIteratorType;
	CovariantIteratorType gradientIt( gradient, region );

	// Make binary air mask
	std::cout << "Making air mask" << std::endl;
	
	ImageType::Pointer airMask = ImageType::New();
	airMask->SetRegions( region );
	airMask->Allocate();
	airMask->SetSpacing( spacing );
	
	IteratorType airMaskIt( airMask, region );

	for (	voxelEdgeIt.GoToBegin(), airMaskIt.GoToBegin();
			!voxelEdgeIt.IsAtEnd() && !airMaskIt.IsAtEnd();
			++voxelEdgeIt, ++airMaskIt)
	{
		if (voxelEdgeIt.Get() == 2) //air
		{
			airMaskIt.Set(1);
		} else {
			airMaskIt.Set(0);
		}
	}

	ss.str("");
	ss << "airMask.hdr";
	WriteITK(airMask, ss.str(), count++);

	// Find edges of air using binary mask
	for (	airMaskIt.GoToBegin(); !airMaskIt.IsAtEnd(); ++airMaskIt)
	{	
		if (airMaskIt.Get() == 1) {	// at air

			ImageType::IndexType index = airMaskIt.GetIndex();
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							for (int k=-1;k<=1;k++) {
								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
									ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};
									if ( airMask->GetPixel( neighborIndex ) == 0 )	// if neighbor is non-air
									{
										airMask->SetPixel( index, 2); // mark as edge
									}
								}
							}
						}
					}
				}
			}
		}
	}

	ss.str("");
	ss << "airMaskEdges.hdr";
	WriteITK( airMask, ss.str(), count++);


	/*
	// Find edges of air mask using sobel edge detection
	std::cout << "Finding edges of air mask" << std::endl;

	typedef itk::SobelEdgeDetectionImageFilter<ImageType, ImageType> SobelFilterType;
	SobelFilterType::Pointer sobelFilter = SobelFilterType::New();
	sobelFilter->SetInput( airMask );

	try  
	{
		sobelFilter->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	}

	ImageType::Pointer edge = sobelFilter->GetOutput();

	ss.str("");
	ss << "sobel_edge_air.hdr";
	WriteITK(edge, ss.str(), count++);

	// Find voxels of air mask touching edge
	IteratorType edgeIt( edge, region );
	for (	edgeIt.GoToBegin(), airMaskIt.GoToBegin();
			!edgeIt.IsAtEnd() && !airMaskIt.IsAtEnd();
			++edgeIt, ++airMaskIt)
	{
		ImageType::IndexType index = edgeIt.GetIndex();
		float edgeValue=edgeIt.Get();
		
		if (edgeValue>0) {
			
			for(int i=-1;i<=1;i++) {
				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-1;j<=1;j++) {
						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							for (int k=-1;k<=1;k++) {
								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
									ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};
									if ( airMask->GetPixel( neighborIndex ) == 1 )	// air
									{
										airMask->SetPixel( neighborIndex, 2); //edge
									}
								}
							}
						}
					}
				}
			}
		}
	}

	ss.str("");
	ss << "airMaskEdges.hdr";
	WriteITK( airMask, ss.str(), count++);

	*/

	/*

	// Compute chamfer distance to stool
	ImageType::Pointer chamferStool = ImageType::New();
	chamferStool->SetRegions(input->GetLargestPossibleRegion() );
	chamferStool->Allocate();
	chamferStool->SetSpacing(input->GetSpacing());
	
	IteratorType chamferStool_it( chamferStool, chamferStool->GetLargestPossibleRegion() );

	for (	voxelEdgeIt.GoToBegin(), chamferStool_it.GoToBegin();
			!voxelEdgeIt.IsAtEnd() && !chamferStool_it.IsAtEnd();
			++voxelEdgeIt, ++chamferStool_it)
	{
		if (voxelEdgeIt.Get() == 1)	//stool
		{
			chamferStool_it.Set(1);	
		} else {
			chamferStool_it.Set(0);
		}
	}

	ss.str("");
	ss << "stoolMask.hdr";
	WriteITK( chamferStool, ss.str(), count++);

	typedef itk::ChamferDistanceTransformImageFilter<ImageType, ImageType> ChamferDistanceFilterType;

	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamferStool);
	int chamfer_weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);//?

	try  
	{
		chamfer_filter->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	}
	
	chamferStool = chamfer_filter->GetOutput();
	//chamfer_filter.~SmartPointer();

	ss.str("");
	ss << "chamferStool_distance.hdr";
	WriteITK( chamferStool, ss.str(), count++);
	*/

	/*
	// Get cubic bspline image interpolation
    InterpolationType::Pointer inputInterpolator = InterpolationType::New();
    inputInterpolator->SetSplineOrder(3);
    inputInterpolator->SetInputImage(input);
	*/

	std::cout << "Number of voxels to move in direction of gradient: " << numOfVoxels << std::endl;

	std::cout << "Finding voxels by gradient" << std::endl;

	int countMatches = 0;
	int countEdges = 0;

	for (	airMaskIt.GoToBegin(), gradientIt.GoToBegin(), inputIt.GoToBegin(); 
			!airMaskIt.IsAtEnd() && !gradientIt.IsAtEnd() && !inputIt.IsAtEnd(); 
			++airMaskIt, ++gradientIt, ++inputIt)
	{
		
		if (airMaskIt.Get() == 2) { // at edge
			countEdges++;
			
			//double theta = -dtheta;

			CovariantImageType::IndexType idx = gradientIt.GetIndex();
			CovariantVectorType gradient = gradientIt.Get();
			gradient.Normalize();

			//for (int k=0; k < 3; k++)
			//{

				CovariantVectorType grad = gradient;

				// Shift grad vector
				//grad[0] = grad[0]*cos(theta) - grad[1]*sin(theta);
				//grad[1] = grad[0]*sin(theta) + grad[1]*cos(theta);

				// Output gradient from filter to text file
				myfile << "Index: (" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t" << inputIt.Get() << "\t";
				myfile << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\t";

				std::vector<ImageType::IndexType> indexVector; //store indices of voxels
				FindVoxelsByGradient(voxelEdge, idx, startIndex, endIndex, grad, numOfVoxels, indexVector);
				
				ss.str("");
				for (int i=0; i < indexVector.size() ; i++)
				{
					ss << voxelEdge->GetPixel( indexVector[i] );
				}

				std::string voxelString = ss.str();

				myfile << ss.str();

				ImageType::IndexType idx2 = idx;
				idx2[1] = 511-idx2[1];
				
				bool isMatch = MatchVoxels( voxelString, regex );

				if (isMatch)
				{
					myfile << "\tMATCH";
					countMatches++;

					// Convert TA/TS to SA transitions
					for (int i=0; i < voxelString.length(); i++)
					{
						if (voxelString.compare(i,1,"6") == 0 || voxelString.compare(i,1,"7") == 0)
						{
							voxelEdgeUpdate->SetPixel( indexVector[i], 5 ); // StoolAir
						}
					}
				}

				myfile << "\n";

				//theta += dtheta;
			//}
		}
	}
	
	myfile << "\n";
	myfile << "Number of air edge voxels patterns: " << countEdges << "\n";
	myfile << "Number of matched patterns: " << countMatches << " (" << 100*countMatches/countEdges << "%)\n";

	ss.str("");
	ss << "voxelEdgeUpdate" << suffix << ".hdr";
	WriteITK( voxelEdgeUpdate, ss.str(), count++);

	myfile.close();

	/*
	// Enforce that all voxels surrounded by a separate class should turn into that class
	for (	voxelEdgeUpdateIt.GoToBegin(); !voxelEdgeUpdateIt.IsAtEnd(); ++voxelEdgeUpdateIt)
	{	
		ImageType::IndexType index = voxelEdgeUpdateIt.GetIndex();
		
		float voxel = voxelEdgeUpdateIt.Get();
		float neighbor = voxel;
		int neighborCount = 1;
		bool sameNeighbors = true;

		for(int i=-1;i<=1;i++) {
			if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0] && sameNeighbors) {
				for (int j=-1;j<=1;j++) {
					if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]  && sameNeighbors) {
						for (int k=-1;k<=1;k++) {
							if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
								
								ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};
								
								// Check if all neighbors are the same
								if ( neighborCount == 1 ) {
									neighbor = voxelEdgeUpdate->GetPixel( neighborIndex ); // Set first neighbor
								} else {
									if ( voxelEdgeUpdate->GetPixel( neighborIndex ) != neighbor ) {
										sameNeighbors = false;
		 								break;
									}
								}
								
								neighborCount++;
							}
						}
					}
				}
			}
		}

		if (sameNeighbors && (voxel != neighbor)) {
			voxelEdgeUpdate->SetPixel(index, neighbor);
			std::cout << "Index: (" << index[0] << ", " << 511-index[1] << ", " << index[2] << ")\n";
		}
	}

	ss.str("");
	ss << "voxelEdgeUpdateNeighbors" << suffix << ".hdr";
	WriteITK( voxelEdgeUpdate, ss.str(), count++);
	*/

	//End the clock
	time_t toc;
	toc = time(NULL);
	int time, min, sec;
	time = difftime(toc,tic);
	time = time%3600;
	min = time/60;
	time = time%60;
	sec = time;
	std::cout << "Elapsed: " << min << " min " << sec << " sec" << std::endl;

	system("pause");
	return 0;
}

bool MatchVoxels(std::string inputString, char *regex) {
	const int ovecount = 6;
	pcre *re;
	const char *error;
	int erroffset;
	int ovector[ovecount];
	int rc;

	//char *regex = "^2*5*6+5*7+5*1*$";
	char *data;
	data = new char[inputString.length()+1];
	strcpy(data, inputString.c_str());

	re = pcre_compile(
	regex, /* the pattern */
	0, /* default options */
	&error, /* for error message */
	&erroffset, /* for error offset */
	NULL); /* use default character table */
	if (! re)
	{
		fprintf(stderr,"PCRE compilation failed at expression offset%d: %s\n", erroffset, error);
		return 1;
	}

	rc = pcre_exec(
	re, /* the compiled pattern */
	NULL, /* no extra data - we didn't study the pattern */
	data, /* the subject string */
	strlen(data), /* the length of the subject */
	0, /* start at offset 0 in the subject */
	0, /* default options */
	ovector, /* output vector for substring information */
	ovecount); /* number of elements in the output
	vector */

	if (rc < 0)
	{
		switch(rc)
		{
			case PCRE_ERROR_NOMATCH:
				//printf("No match found in text\n");
				break;
			default:
				//printf("Match error %d\n", rc);
				break;

			return false;
		}
	} else {
		return true;
	}

	/*
	// Print all matches
	std::cout << "Number of matches: " << rc << std::endl;
	for (int i=0; i < rc; i++)
	{
		printf("Match %d: %.*s\n", i, ovector[2*i + 1] - ovector[2*i], data + ovector[2*i]);
	}
	*/

}

void FindVoxelsByGradient(ImageType::Pointer voxelEdge, ImageType::IndexType index, ImageType::IndexType startIndex, ImageType::IndexType endIndex, CovariantVectorType grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector) {
	//Calculates indices of voxels in direction of gradient
	
	ContinuousIndexType offsetIndex;

	// Space increment
	double ds = 0.1;

	int i=0;
	ImageType::PixelType voxel = voxelEdge->GetPixel( index );

	while ( (i < numOfVoxels) && (voxel != 1) && (voxel != 3) ) // stop finding once you hit stool or tissue
	{
		// Reset new offset index
		offsetIndex[0] = index[0];
		offsetIndex[1] = index[1];
		offsetIndex[2] = index[2];
		
		// Increment index in direction of gradient
		do {
			offsetIndex[0] += grad[0]*ds;
			offsetIndex[1] += grad[1]*ds;
			offsetIndex[2] += grad[2]*ds;
			
			// Ensure new index is within image boundaries
			for (int j=0; j < 3; j++)
			{
				if (offsetIndex[j] < startIndex[j]) {	offsetIndex[j] = startIndex[j]; }
				if (offsetIndex[j] > endIndex[j])   {	offsetIndex[j] = endIndex[j]; }
			}

		} while (	(index[0] == round(offsetIndex[0])) &&
					(index[1] == round(offsetIndex[1])) &&
					(index[2] == round(offsetIndex[2]))	);	// Until a new discrete voxel has been found
		
		// Set new index and store it
		for (int k=0; k < 3; k++) { index[k] = round(offsetIndex[k]); }

		voxel = voxelEdge->GetPixel( index );
		
		indexVector.push_back(index);

		i++;
	}	
}

void FindVoxelsByGradient2(ImageType::Pointer voxelEdge, ImageType::IndexType index, ImageType::IndexType startIndex, ImageType::IndexType endIndex, CovariantVectorType grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector) {
	int i=0;
	ImageType::PixelType voxel = voxelEdge->GetPixel( index );

	while ( (i < numOfVoxels) && (voxel != 1) && (voxel != 3) ) // stop finding once you hit stool or tissue
	{
		index[0] += round( grad[0] );
		index[1] += round( grad[1] );
		index[2] += round( grad[2] );

		voxel = voxelEdge->GetPixel( index );
		
		indexVector.push_back(index);

		i++;
	}	
}

void WriteITK(ImageType::Pointer image, std::string ss, int count) {
	std::cout << "Writing " << ss << std::endl;
	
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	
	std::stringstream ss2;
	ss2 << count << "_" << ss;
	writer->SetFileName(ss2.str().c_str());
	writer->SetInput(image);
	writer->GlobalWarningDisplayOff();

	try  
	{
		writer->Update();
	} catch( itk::ExceptionObject & excp ) 
	{
		std::cerr << excp << std::endl;
	} 
}

void ReadITK(ImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	image = reader->GetOutput();
}

double round(float d)
{
  return floor(d + 0.5);
}