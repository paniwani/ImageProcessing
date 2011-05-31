#include "RemoveStool4.h"
#include "io.cxx"
#include "scattercorrection.cxx"
#include "QR.cxx"
#include "EM.cxx"
//#include "HessianFunctions.cxx"

int main(int argc, char * argv[])
{
	/*
	*	TODO
	*	- smooth gradient
	*	- colon air set to 0?
	*
	*/

	clock_t init,final;
	init = clock();


	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}

	ImageType::Pointer				inputOriginal			= ImageType::New();
	ImageType::Pointer				input					= ImageType::New();
	ByteImageType::Pointer			colon					= ByteImageType::New();
	FloatImageType::Pointer			gradientMagnitude		= FloatImageType::New();
	VoxelImageType::Pointer			vmap					= VoxelImageType::New();
	ArrayImageType::Pointer			partial					= ArrayImageType::New();


	debug.open("debug.txt");

	// Load images and segment colon
	Setup(argv[1],inputOriginal,input,colon,gradientMagnitude);

	// Initial segmentation, save tissue stool threshold
	PixelType tst = SingleMaterialClassification(input, gradientMagnitude, vmap, colon);

	// Apply scatter correction
	input = ScatterCorrection(inputOriginal,colon,vmap);

	// Update vmap with new scatter input
	ApplyThresholdRules(input,gradientMagnitude,vmap,colon,tst);

	Write(vmap,"scatter_vmap.nii");

	//DirectionalGradient(input, colon, vmap);

	//TextureAnalysis(input);

	//// Determine boundary types
	partial = QuadraticRegression(input,colon,vmap,gradientMagnitude,tst);

	//// EM
	EM(partial,colon,input);

	// end clock
	final = clock() - init;

	std::cout << (double) final / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;

	system("pause");
	return 0;
}

//void TextureAnalysis(ImageType::Pointer &input)
//{
//	// get gradient (no smoothing)
//	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
//	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
//	gmFilter->SetInput(input);
//
//	// threshold gradient
//	typedef itk::BinaryThresholdImageFilter<FloatImageType,ByteImageType> BinaryThresholdImageFilterType;
//	BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
//	thresholdFilter->SetInsideValue(1);
//	thresholdFilter->SetOutsideValue(0);
//	thresholdFilter->SetLowerThreshold(400);
//	thresholdFilter->SetInput(gmFilter->GetOutput());
//	thresholdFilter->Update();
//
//	// write mask
//	typedef itk::MultiplyByConstantImageFilter<ByteImageType,unsigned char,ByteImageType> MultiplyFilterType;
//	MultiplyFilterType::Pointer multiplier = MultiplyFilterType::New();
//	multiplier->SetInput(thresholdFilter->GetOutput());
//	multiplier->SetConstant(255);
//	multiplier->Update();
//	Write(multiplier->GetOutput(),"mask.nii");
//
//	// get image minimum and maximum
//	typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
//	MinimumMaximumImageCalculatorType::Pointer imageCalc = MinimumMaximumImageCalculatorType::New();
//	imageCalc->SetImage(input);
//	imageCalc->ComputeMinimum();
//	imageCalc->ComputeMaximum();
//
//	typedef otb::ScalarImageToTexturesFilter2<ImageType,FloatImageType,ByteImageType> TextureFilterType;
//	TextureFilterType::Pointer textureFilter = TextureFilterType::New();
//	textureFilter->SetInput(input);
//	textureFilter->SetMaskImage(thresholdFilter->GetOutput());
//	textureFilter->SetInputImageMinimum(imageCalc->GetMinimum());
//	textureFilter->SetInputImageMaximum(imageCalc->GetMaximum());
//	
//	ImageType::OffsetType offset;
//	offset[0] = 1;
//	offset[1] = 0;
//	offset[2] = 0;
//
//	textureFilter->SetOffset(offset);
//
//	ImageType::SizeType radius;
//	radius[0] = 1;
//	radius[1] = 1;
//	radius[2] = 0;
//
//	textureFilter->SetRadius(radius);
//	textureFilter->Update();
//
//	for (int i=0; i<8; i++)
//	{
//		std::stringstream ss;
//		ss << "texture" << i << ".nii";
//		Write(textureFilter->GetOutput(i),ss.str());
//	}
//
//}

void DirectionalGradient(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap)
{
	// create binary tissue mask
	ByteImageType::Pointer mask = ByteImageType::New();
	mask->SetSpacing(input->GetSpacing());
	mask->SetRegions(REGION);
	mask->Allocate();
	mask->FillBuffer(0);
	
	ByteIteratorType maskIt(mask,REGION);
	VoxelIteratorType vmapIt(vmap,REGION);

	for (maskIt.GoToBegin(), vmapIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++vmapIt)
	{
		if (vmapIt.Get() == Air)
			maskIt.Set(255);
	}

	Write(mask,"airMask.nii");

	// mask dg adjacent to air only

	typedef itk::ImageDuplicator<ByteImageType> ImageDuplicatorType;
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
	duplicator->SetInputImage(mask);
	duplicator->Update();
	ByteImageType::Pointer maskDilated = duplicator->GetOutput();

	Dilate(maskDilated,3);
	Write(maskDilated,"airMaskDilated.nii");

	typedef itk::SubtractImageFilter<ByteImageType,ByteImageType,ByteImageType> SubtractImageFilterType;
	SubtractImageFilterType::Pointer subtracter = SubtractImageFilterType::New();
	subtracter->SetInput1(maskDilated);
	subtracter->SetInput2(mask);
	subtracter->Update();
	ByteImageType::Pointer airBorder = subtracter->GetOutput();

	// only see dg in unclassified air border
	ByteIteratorType airBorderIt(airBorder,REGION);

	for (airBorderIt.GoToBegin(), vmapIt.GoToBegin(); !airBorderIt.IsAtEnd(); ++airBorderIt, ++vmapIt)
	{
		if (vmapIt.Get() != Unclassified)
			airBorderIt.Set(0);
	}

	Write(airBorder,"airBorder.nii");

	// directional gradient
	typedef itk::DirectionalGradientImageFilter2<ImageType,ByteImageType,ImageType> DirectionalGradientImageFilterType;
	DirectionalGradientImageFilterType::Pointer dgFilter = DirectionalGradientImageFilterType::New();
	dgFilter->SetInput(input);
	//dgFilter->SetSigma(input->GetSpacing()[0]);
	dgFilter->SetScale(-1);
	dgFilter->SetMaskImage(mask);
	dgFilter->SetOutsideValue(255);
	dgFilter->Update();
	ImageType::Pointer dg = dgFilter->GetOutput();
	IteratorType dgIt(dg,REGION);

	// set -ve to 0
	for (dgIt.GoToBegin(); !dgIt.IsAtEnd(); ++dgIt)
	{
		if (dgIt.Get() < 0)
			dgIt.Set(0);
	}

	Write(dg,"dg.nii");
	
	typedef itk::MaskImageFilter<ImageType,ByteImageType,ImageType> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(dg);
	masker->SetInput2(airBorder);
	masker->Update();
	Write(masker->GetOutput(),"dgMasked.nii");
	
}

void Dilate(ByteImageType::Pointer &img, unsigned int radius)
{
	StructuringElementType se;
	
	/*ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;*/

	se.SetRadius( radius );
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<ByteImageType, ByteImageType, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( img );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();

	img = dilater->GetOutput();
}


/*********************************************************
- Load image
- Shift input so that air is 0 HU
- Segment colon
- Mask input with colon
- Crop input and colon in XY plane
*********************************************************/
void Setup(std::string dataset, ImageType::Pointer  &inputOriginal, ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude)
{
	//----------------------------------------------
	// Load image
	//----------------------------------------------
	std::vector<std::string> datasetArray = explode( "\\", dataset );
	std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	if (truncateOn)
	{
		std::cout << "Truncating data" << std::endl;
		std::cout << "Slices: " << truncateArray[0] << " to " << truncateArray[1] << std::endl;
		
		std::stringstream ss;
		ss << dsname << "_" << truncateArray[0] << "_" << truncateArray[1];
		dsname = ss.str();
	}

	// Set writer prefix
	note = dsname;

	// Load dicom files
	if (!truncateOn)
	{
		inputOriginal = ReadDicom < ImageType > ( dataset );
	} else {
		inputOriginal = ReadDicom < ImageType > ( dataset, truncateArray[0], truncateArray[1] );
	}

	// Set global region
	REGION = inputOriginal->GetLargestPossibleRegion();	

	Write(inputOriginal,"inputOriginal.nii");

	// Get image minimum and set as global background
	typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minCalc = MinimumMaximumImageCalculatorType::New();
	minCalc->SetImage(inputOriginal);
	minCalc->ComputeMinimum();
	BACKGROUND = minCalc->GetMinimum();

	//----------------------------------------------
	// Segment colon
	//----------------------------------------------
	typedef itk::ColonSegmentationFilter< ImageType, ByteImageType > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( inputOriginal );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	//colonSegmenter->SetPrintImages(true);

	if ( truncateOn )
		colonSegmenter->SetRemoveBoneLung( false );

	colonSegmenter->Update();
	colon = colonSegmenter->GetOutput();

	//----------------------------------------------
	// Mask input with colon
	//----------------------------------------------
	typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( inputOriginal );
	masker->SetInput2( colon );
	masker->SetOutsideValue( BACKGROUND );
	masker->Update();
	input = masker->GetOutput();

	//----------------------------------------------
	// Calculate gradient magnitude
	//----------------------------------------------
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientMagnitudeFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientMagnitudeFilter->SetInput( inputOriginal );
	gradientMagnitudeFilter->SetSigma( inputOriginal->GetSpacing()[0] );
	gradientMagnitudeFilter->Update();
	gradientMagnitude = gradientMagnitudeFilter->GetOutput();

	//----------------------------------------------
	// Crop images in XY plane
	//----------------------------------------------

	ImageType::SizeType size = REGION.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0;
	short paddingXY = 5;
	
	IteratorType inputIt(input,REGION);
	ByteIteratorType colonIt(colon,REGION);

	for (inputIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt)
	{
		if ( colonIt.Get() != 0 )
		{
			ImageType::IndexType idx = inputIt.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];
		}
	}

	ImageType::IndexType edx;
	
	edx[0] = (minX-paddingXY) > 0 ? minX-paddingXY : 0;
	edx[1] = (minY-paddingXY) > 0 ? minY-paddingXY : 0;
	edx[2] = 0;

	ImageType::SizeType esize;
	esize[0] = maxX-minX+2*paddingXY+1 < size[0] ? maxX-minX+2*paddingXY+1 : size[0];
	esize[1] = maxY-minY+2*paddingXY+1 < size[1] ? maxY-minY+2*paddingXY+1 : size[1];
	esize[2] = size[2];

	ImageType::RegionType extractRegion;
	extractRegion.SetIndex( edx );
	extractRegion.SetSize( esize );

	OLDREGION = extractRegion;

	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( input );
	cropper->SetRegionOfInterest( extractRegion );
	cropper->Update();
	input = cropper->GetOutput();

	typedef itk::RegionOfInterestImageFilter<FloatImageType,FloatImageType> RegionOfInterestImageFilterFloatType;
	RegionOfInterestImageFilterFloatType::Pointer cropperFloat = RegionOfInterestImageFilterFloatType::New();
	cropperFloat->SetInput( gradientMagnitude );
	cropperFloat->SetRegionOfInterest( extractRegion );
	cropperFloat->Update();
	gradientMagnitude = cropperFloat->GetOutput();

	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();

	// Set cropped region globally
	REGION = colon->GetLargestPossibleRegion();

	Write(input,"input.nii");
	Write(colon,"colon.nii");
	Write(gradientMagnitude,"gradientMagnitude.nii");
}

/*********************************************************
- Allocate voxel map
- Compute otsu threshold separating tissue and stool
- Apply initial threshold rules
*********************************************************/
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon) 
{
	//----------------------------------------------
	// Allocate voxel map
	//----------------------------------------------
	vmap->SetRegions( input->GetLargestPossibleRegion() );
	vmap->SetSpacing( input->GetSpacing() );
	vmap->SetDirection( input->GetDirection() );
	vmap->CopyInformation( input );
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);

	//----------------------------------------------
	// Compute otsu threshold separating tissue and stool
	//----------------------------------------------

	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1500);
	otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	otsuCalculator->Compute();

	PixelType tissueStoolThreshold = otsuCalculator->GetThreshold();
	std::cout << "Tissue Stool Otsu Threshold: " << tissueStoolThreshold << std::endl;

	//----------------------------------------------
	// Apply initial threshold rules
	//----------------------------------------------
	ApplyThresholdRules( input, gradientMagnitude, vmap, colon, tissueStoolThreshold );

	Write(vmap,"vmap.nii");

	return tissueStoolThreshold;
}

void ApplyThresholdRules( ImageType::Pointer &input, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissueStoolThreshold )
{
	IteratorType inputIt(input,REGION);
	FloatIteratorType gradientMagnitudeIt(gradientMagnitude,REGION);
	VoxelIteratorType vmapIt(vmap,REGION);
	ByteIteratorType colonIt(colon,REGION);

	for ( inputIt.GoToBegin(), gradientMagnitudeIt.GoToBegin(), vmapIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++gradientMagnitudeIt, ++vmapIt, ++colonIt)
	{
		if (colonIt.Get() == 255)
		{
			VoxelType voxel;

			short I = inputIt.Get();
			float G = gradientMagnitudeIt.Get();

			if ( ( I >= tissueStoolThreshold && G < 0.8*I ) || I > 1000 )
			{
				voxel = Stool;
			} else if ( I <= -600 ) {
				voxel = Air;
			} else if ( I < tissueStoolThreshold && I > -300 && G <= 400 ) {
				voxel = Tissue;
			} else {
				voxel = Unclassified;
			}

			vmapIt.Set( voxel );
		}
	}
}