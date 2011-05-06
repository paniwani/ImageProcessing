#include "RemoveStool4.h"
#include "io.cxx"
#include "scattercorrection.cxx"
#include "QR.cxx"
//#include "HessianFunctions.cxx"

int main(int argc, char * argv[])
{
	/*
	*	TODO
	*	- smooth gradient
	*	- colon air set to 0?
	*
	*/


	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}

	ImageType::Pointer				input_original			= ImageType::New();
	ImageType::Pointer				input					= ImageType::New();
	ByteImageType::Pointer			colon					= ByteImageType::New();
	FloatImageType::Pointer			gradient_magnitude		= FloatImageType::New();
	VoxelImageType::Pointer			vmap					= VoxelImageType::New();
	ArrayImageType::Pointer			partial					= ArrayImageType::New();

	// Load images and segment colon
	Setup(argv[1],input_original,input,colon,gradient_magnitude);

	// Initial segmentation, save tissue stool threshold
	PixelType tst = SingleMaterialClassification(input, gradient_magnitude, vmap, colon);

	// Apply scatter correction
	input = ScatterCorrection(input_original,colon,vmap);

	// Update vmap with new scatter input
	ApplyThresholdRules(input,gradient_magnitude,vmap,colon,tst);

	Write(vmap,"scatter_vmap.nii");

	// Determine boundary types
	partial = QuadraticRegression(input,colon,vmap,gradient_magnitude);

	system("pause");
	return 0;
}

void Dilate(ByteImageType::Pointer &img, unsigned int radius)
{
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;

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
void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradient_magnitude)
{
	//----------------------------------------------
	// Load image
	//----------------------------------------------
	std::vector<std::string> dataset_ar = explode( "\\", dataset );
	std::string dsname = dataset_ar[ dataset_ar.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	if (truncate_on)
	{
		std::cout << "Truncating data" << std::endl;
		std::cout << "Slices: " << truncate_ar[0] << " to " << truncate_ar[1] << std::endl;
		
		std::stringstream ss;
		ss << dsname << "_" << truncate_ar[0] << "_" << truncate_ar[1];
		dsname = ss.str();
	}

	// Set writer prefix
	note = dsname;

	ShortImageType::Pointer input_s;

	// Load dicom files
	if (!truncate_on)
	{
		input_s = ReadDicom < ShortImageType > ( dataset );
	} else {
		input_s = ReadDicom < ShortImageType > ( dataset, truncate_ar[0], truncate_ar[1] );
	}

	//Write(input_s,"input_s.nii");

	ShortIteratorType input_s_iter( input_s, input_s->GetLargestPossibleRegion() );

	//----------------------------------------------
	// Shift input so that air is 0 HU
	//----------------------------------------------

	// Shift and cast to 16-bit unsigned short
	// Save original (input_original) and modified input with colon mask applied (input) 
	input_original->SetRegions( input_s->GetLargestPossibleRegion() );
	input_original->SetSpacing( input_s->GetSpacing() );
	input_original->SetDirection( input_s->GetDirection() );
	input_original->CopyInformation( input_s );
	input_original->Allocate();
	IteratorType input_original_iter(input_original, input_original->GetLargestPossibleRegion() );
	
	input->SetRegions( input_s->GetLargestPossibleRegion() );
	input->SetSpacing( input_s->GetSpacing() );
	input->SetDirection( input_s->GetDirection() );
	input->CopyInformation( input_s );
	input->Allocate();
	IteratorType input_iter(input, input->GetLargestPossibleRegion() );

	// Compute background pixel value and save globally
	typedef itk::MinimumMaximumImageCalculator< ShortImageType > MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer calc = MinimumMaximumImageCalculatorType::New();
	calc->SetImage( input_s );
	calc->ComputeMinimum();
	BACKGROUND = calc->GetMinimum();

	for (input_s_iter.GoToBegin(), input_original_iter.GoToBegin(); !input_s_iter.IsAtEnd(); ++input_s_iter, ++input_original_iter)
	{
		input_original_iter.Set( input_s_iter.Get() + abs( BACKGROUND ) );
	}

	input_s.~SmartPointer();

	// Set global region
	REGION = input_original->GetLargestPossibleRegion();

	//----------------------------------------------
	// Segment colon
	//----------------------------------------------
	typedef itk::ColonSegmentationFilter< ImageType, ByteImageType > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colon_segmenter = ColonSegmentationFilterType::New();
	colon_segmenter->SetInput( input_original );
	colon_segmenter->SetForegroundValue( 255 );
	colon_segmenter->SetBackgroundValue( 0 );
	//colon_segmenter->SetPrintImages(true);

	if ( truncate_on )
		colon_segmenter->SetRemoveBoneLung( false );

	colon_segmenter->Update();
	colon = colon_segmenter->GetOutput();

	//----------------------------------------------
	// Mask input with colon
	//----------------------------------------------
	typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( input_original );
	masker->SetInput2( colon );
	masker->Update();

	input = masker->GetOutput();

	Write(input_original,"input_original.nii");

	//----------------------------------------------
	// Calculate gradient magnitude and mask with colon
	//----------------------------------------------
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	
	GradientMagnitudeImageFilterType::Pointer gradient_magnitude_filter = GradientMagnitudeImageFilterType::New();
	gradient_magnitude_filter->SetInput( input_original );

	gradient_magnitude_filter->SetSigma( input_original->GetSpacing()[0] );

	gradient_magnitude_filter->Update();
	gradient_magnitude = gradient_magnitude_filter->GetOutput();

	//----------------------------------------------
	// Crop images in XY plane
	//----------------------------------------------
	
	input_iter = IteratorType(input,REGION);

	ImageType::SizeType size = REGION.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0;
	unsigned short paddingXY = 5;

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		if ( input_iter.Get() != 0 )
		{
			ImageType::IndexType idx = input_iter.GetIndex();

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
	cropperFloat->SetInput( gradient_magnitude );
	cropperFloat->SetRegionOfInterest( extractRegion );
	cropperFloat->Update();
	gradient_magnitude = cropperFloat->GetOutput();

	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();

	REGION = colon->GetLargestPossibleRegion();

	Write(input,"input.nii");
	Write(colon,"colon.nii");
	Write(gradient_magnitude,"gradient_magnitude.nii");
}

/*********************************************************
1. Calculate gradient magnitude
2. Allocate voxel map
3. Compute otsu threshold separating tissue and stool
4. Apply initial threshold rules
*********************************************************/
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon) 
{
	

	/*typedef itk::CastImageFilter<ImageType,FloatImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput(input);
	caster->Update();
	*/
	//// No smoothing
	//typedef itk::GradientMagnitudeImageFilter<FloatImageType,FloatImageType> GradientMagnitudeImageFilterNoSmoothingType;
	//GradientMagnitudeImageFilterNoSmoothingType::Pointer grad_no_smooth_filter = GradientMagnitudeImageFilterNoSmoothingType::New();
	//grad_no_smooth_filter->SetInput( caster->GetOutput() );
	//grad_no_smooth_filter->Update();
	//Write(grad_no_smooth_filter->GetOutput(), "grad_no_smooth.nii");

	//// Sobel
	//typedef itk::SobelEdgeDetectionImageFilter<FloatImageType,FloatImageType> SobelEdgeDetectionImageFilterType;
	//SobelEdgeDetectionImageFilterType::Pointer sobel_filter = SobelEdgeDetectionImageFilterType::New();
	//sobel_filter->SetInput( caster->GetOutput() );
	//sobel_filter->Update();
	//Write(sobel_filter->GetOutput(),"sobel.nii");

	//----------------------------------------------
	// 2. Allocate voxel map
	//----------------------------------------------
	vmap->SetRegions( input->GetLargestPossibleRegion() );
	vmap->SetSpacing( input->GetSpacing() );
	vmap->SetDirection( input->GetDirection() );
	vmap->CopyInformation( input );
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);
	
	IteratorType input_iter(input,REGION);
	FloatIteratorType gradient_magnitude_iter(gradient_magnitude,REGION);
	VoxelIteratorType vmap_iter(vmap,REGION);
	ByteIteratorType colon_iter(colon,REGION);

	//----------------------------------------------
	// 3. Compute otsu threshold separating tissue and stool
	//----------------------------------------------

	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(775);
	otsuCalculator->SetHistogramMax(2425);
	otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	otsuCalculator->Compute();

	PixelType tissue_stool_threshold = otsuCalculator->GetThreshold();

	std::cout << std::endl << std::endl;
	std::cout << "Tissue Stool Otsu Threshold: " << tissue_stool_threshold << std::endl;
	std::cout << std::endl << std::endl;

	//----------------------------------------------
	// 4. Apply initial threshold rules
	//----------------------------------------------
	ApplyThresholdRules( input, gradient_magnitude, vmap, colon, tissue_stool_threshold );

	Write(vmap,"vmap.nii");

	return tissue_stool_threshold;
}

void ApplyThresholdRules( ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissue_stool_threshold )
{
	IteratorType input_iter(input,REGION);
	FloatIteratorType gradient_magnitude_iter(gradient_magnitude,REGION);
	VoxelIteratorType vmap_iter(vmap,REGION);
	ByteIteratorType colon_iter(colon,REGION);

	tissue_stool_threshold += BACKGROUND;

	for ( input_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gradient_magnitude_iter, ++vmap_iter, ++colon_iter)
	{
		if (colon_iter.Get() == 255)
		{
			VoxelType voxel;

			short I = input_iter.Get() + BACKGROUND;
			float G = gradient_magnitude_iter.Get();

			if ( ( I >= tissue_stool_threshold && G < 0.8*I ) || I > 1000 )
			{
				voxel = Stool;
			} else if ( I <= -600 ) {
				voxel = Air;
			} else if ( I < tissue_stool_threshold && I > -300 && G <= 400 ) {
				voxel = Tissue;
			} else {
				voxel = Unclassified;
			}

			vmap_iter.Set( voxel );
		}
	}
}