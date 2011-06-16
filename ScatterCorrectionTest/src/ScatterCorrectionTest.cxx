/*
Test scatter correction parameters
*/
const unsigned int Dimension = 3;



#include <itkImage.h> 						
#include <iostream> 			
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkGradientMagnitudeImageFilter.h>

#include <itkAddConstantToImageFilter.h>
PixelType BACKGROUND = -1025;
ImageType::RegionType REGION;
ImageType::RegionType OLDREGION;

#include <scattercorrection.cxx>



									
int main(int argc, char * argv[])				
{ 					
	// Get args
	if( argc < 6 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory slice1 slice2 filterP scaleP SCALE";
		system("pause");
		return EXIT_FAILURE;
	}

	// Load image
	ImageType::Pointer inputOriginal = ReadDicom <ImageType> (argv[1],atoi(argv[2]),atoi(argv[3]));
	ImageType::RegionType region = inputOriginal->GetLargestPossibleRegion();

	// Load parameters
	filterP = atof(argv[4]);
	scaleP = atoi(argv[5]);
	SCALE = atoi(argv[6]);
		
	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(inputOriginal);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	ImageType::Pointer input = CropByRegion <ImageType> (inputOriginal,colonRegion);
	Mask <ImageType,ByteImageType> (input,colon,-1025);

	OLDREGION = colonRegion;

	WriteITK <ImageType> (inputOriginal,"inputOriginal.nii");
	WriteITK <ImageType> (input,"input.nii");
	//WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ImageType::SizeType size = region.GetSize();

	REGION = region;

	// Scatter correction
	ImageType::Pointer scatterInput = ScatterCorrection(inputOriginal,colon);
	std::stringstream ss;
	ss << "scatterInput_filterP_" << filterP << "_scaleP_" << scaleP << "_SCALE_" << SCALE << ".nii";
	WriteITK <ImageType> (scatterInput,ss.str());

	// Get difference image
	typedef itk::SubtractImageFilter<ImageType> SubtracterType;
	SubtracterType::Pointer subtracter = SubtracterType::New();
	subtracter->SetInput1(input);
	subtracter->SetInput2(scatterInput);
	subtracter->Update();

	ss.str("");
	ss << "scatterDiff_filterP_" << filterP << "_scaleP_" << scaleP << "_SCALE_" << SCALE << ".nii";
	WriteITK <ImageType> (subtracter->GetOutput(),ss.str());

	return 0;
}