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

	// Load images and segment colon
	Setup(argv[1],input_original,input,colon);

	// Initial segmentation, save tissue stool threshold
	PixelType tst = SingleMaterialClassification(input, gradient_magnitude, vmap, colon);

	// Apply scatter correction
	input = ScatterCorrection(input_original,colon,vmap);

	// Update vmap with new scatter input
	ApplyThresholdRules(input,gradient_magnitude,vmap,colon,tst);

	Write(vmap,"scatter_vmap.nii");

	// Determine boundary types
	QuadraticRegression(input,colon,vmap,gradient_magnitude);

	system("pause");
	return 0;
}

/*********************************************************
Using intensity vs gradient relationship, assign boundary
types to edge voxels

1. Determine egdes using canny
2. Get normalized gradient vector
3. Interpolate input and gradient magnitude
4. Compute local stool maximum
5. Run voxel edge classification

*********************************************************/
void QuadraticRegression(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, FloatImageType::Pointer &gradient_magnitude)
{
	//----------------------------------------------
	// 1. Determine egdes using canny
	//----------------------------------------------

	//typedef itk::CastImageFilter<ImageType, FloatImageType> CastImageFilterType;
	//CastImageFilterType::Pointer caster = CastImageFilterType::New();
	//caster->SetInput( input );

	//typedef itk::CannyEdgeDetectionImageFilterModified<FloatImageType,FloatImageType> CannyEdgeDetectionImageFilterType;
	//CannyEdgeDetectionImageFilterType::Pointer cannyFilter = CannyEdgeDetectionImageFilterType::New();
	//cannyFilter->SetInput( caster->GetOutput() );
	//cannyFilter->SetLowerThreshold(0); // minimum gradient threshold
	//cannyFilter->SetVariance(1);	
	//cannyFilter->Update();

	//Write(cannyFilter->GetOutput(), "canny.nii");

	/*typedef itk::CastImageFilter<FloatImageType, ByteImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput( cannyFilter->GetOutput() );
	caster->Update();
	
	ByteImageType::Pointer canny = caster->GetOutput();

	Dilate(canny, 1);*/

	//----------------------------------------------
	// 2. Get normalized gradient vector 
	//----------------------------------------------


	// output gradient vector to text file
	std::ofstream file;
	file.open("grad.txt");


	typedef itk::GradientImageFilter<ImageType> GradientImageFilterType;
	GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->Update();
	VectorImageType::Pointer gradient = gradientFilter->GetOutput();
	VectorIteratorType gradient_iter(gradient,REGION);

	ImageType::SpacingType spacing = input->GetSpacing();
	IteratorType input_iter(input,REGION);

	for (gradient_iter.GoToBegin(), input_iter.GoToBegin(); !gradient_iter.IsAtEnd(); ++gradient_iter, ++input_iter)
	{
		VectorType grad = gradient_iter.Get();

		float norm = grad.GetNorm();
		
		if (norm > 0)
		{
			grad[0] /= norm;
			grad[1] /= norm;
			grad[2] /= norm;
		}
		
		// Set tolerance
		for (int i=0; i<3; i++)
			grad[i] = abs(grad[i]) > 0.001 ? grad[i] : 0;

		gradient_iter.Set(grad);

		ImageType::IndexType idx = gradient_iter.GetIndex();
		idx[0] -= REGION.GetIndex()[0];
		idx[1] -= REGION.GetIndex()[1];
		idx[2] -= REGION.GetIndex()[2];

		if (input_iter.Get() == 721)
		{

			file << idx[0]*spacing[0] << "\t" << idx[1]*spacing[1] << "\t" << idx[2]*spacing[2] << "\t";
			file << input_iter.Get() << "\t";
			file << grad[0] << "\t" << grad[1] << "\t" << grad[2] << "\n";
		
		}
	}

	file.close();

	// write grad as .mhd
	typedef itk::ImageFileWriter<VectorImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(gradient);
	writer->SetFileName("gradient.mhd");
	writer->Update();

	Write(input,"input.mhd");


	//----------------------------------------------
	// 3. Interpolate input and gradient magnitude
	//----------------------------------------------
	
	// Convert input to short because B spline interp does not work well with unsigned data
	// Shift back to air at -1024
	ShortImageType::Pointer input_short = ShortImageType::New();
	input_short->SetRegions(REGION);
	input_short->SetSpacing(input->GetSpacing());
	input_short->SetDirection(input->GetDirection());
	input_short->CopyInformation(input);
	input_short->Allocate();
	ShortIteratorType input_short_iter(input_short,REGION);

	for (input_iter.GoToBegin(), input_short_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++input_short_iter)
	{
		input_short_iter.Set( input_iter.Get()+BACKGROUND );
	}

	Write(input_short,"input_short.nii");

	typedef itk::BSplineInterpolateImageFunction<ShortImageType> BSplineInterpolateImageFunctionImageType;
	BSplineInterpolateImageFunctionImageType::Pointer input_interp = BSplineInterpolateImageFunctionImageType::New();
	input_interp->SetSplineOrder(3);
	input_interp->SetInputImage(input_short);

	ShortImageType::IndexType index = input_interp->GetStartIndex();
	std::cout << index[0] << " " << index[1] << " " << index[2] << std::endl;

	typedef itk::BSplineInterpolateImageFunction<FloatImageType> BSplineInterpolateImageFunctionFloatImageType;
	BSplineInterpolateImageFunctionFloatImageType::Pointer gradient_magnitude_interp = BSplineInterpolateImageFunctionFloatImageType::New();
	gradient_magnitude_interp->SetSplineOrder(3);
	gradient_magnitude_interp->SetInputImage(gradient_magnitude);

	//----------------------------------------------
	// 4. Compute local stool maximum
	//----------------------------------------------
	ImageType::Pointer smax = ComputeNeighborhoodSmax(input,vmap,colon);
	IteratorType smax_iter(smax,REGION);

	//----------------------------------------------
	// 5. Run voxel edge classification
	//----------------------------------------------

	// Perform QR over several distance pairs
	float d1[3] = {1.0,0.5,0.3};
	float d2[3] = {1.5,1.0,0.6};

	VoxelIteratorType vmap_iter(vmap,REGION);
	ByteIteratorType colon_iter(colon,REGION);
	//IteratorType input_iter(input,REGION);

	ContinuousIndexType idx1 = input_interp->GetStartContinuousIndex();
	ContinuousIndexType idx2 = input_interp->GetEndContinuousIndex();

	std::cout << idx1[0] << " " << idx1[1] << " " << idx1[2] << std::endl;
	std::cout << idx2[0] << " " << idx2[1] << " " << idx2[2] << std::endl;

	int count=0;

	for (vmap_iter.GoToBegin(), smax_iter.GoToBegin(), gradient_iter.GoToBegin(), colon_iter.GoToBegin(), input_short_iter.GoToBegin(); !vmap_iter.IsAtEnd();
		++vmap_iter, ++smax_iter, ++gradient_iter, ++colon_iter, ++input_short_iter)
	{
		if (colon_iter.Get() == 255 && vmap_iter.Get() == Unclassified)
		{
			ImageType::IndexType index		=	input_short_iter.GetIndex();
			VectorType grad					=	gradient_iter.Get();

			VoxelType voxel = Unclassified;

			if (count > 1756)
			{
				std::cout << "stop here" << std::endl;

				ImageType::IndexType sx = REGION.GetIndex();
				std::cout << index[0]-sx[0] << " " << index[1]-sx[1] << " " << index[2]-sx[2] << std::endl;
			}
			
			if (smax_iter.Get() > 0)
			{

				//std::cout << grad[0] << "\t" << grad[1] << "\t" << grad[2] << std::endl;


				float distance = -1;

				for (int i=0; i<3; i++)
				{
					// Compute average orthogonal distance to each parabola
					
					// Create offset index
					BSplineInterpolateImageFunctionFloatImageType::ContinuousIndexType offset_index;
					offset_index[0]=index[0];
					offset_index[1]=index[1];
					offset_index[2]=index[2];
					
					// Store intensity and gradient magnitude for 5 locations
					std::vector< double > temp_intensity;
					std::vector< double > temp_gradient_magnitude;

					temp_intensity.push_back( input_interp->EvaluateAtContinuousIndex(offset_index) );
					temp_gradient_magnitude.push_back( gradient_magnitude_interp->EvaluateAtContinuousIndex(offset_index) );
					
					// Add values in positive direction,
					// +d1
					offset_index[0]=index[0]+grad[0]*d1[i];
					offset_index[1]=index[1]+grad[1]*d1[i];
					offset_index[2]=index[2]+grad[2]*d1[i];

					if ( checkBounds(REGION, offset_index) )
					{
						temp_intensity.push_back( input_interp->EvaluateAtContinuousIndex(offset_index) );
						temp_gradient_magnitude.push_back( gradient_magnitude_interp->EvaluateAtContinuousIndex(offset_index) );

						// +d2
						offset_index[0]=index[0]+grad[0]*d2[i];
						offset_index[1]=index[1]+grad[1]*d2[i];
						offset_index[2]=index[2]+grad[2]*d2[i];

						if ( checkBounds(REGION, offset_index) )
						{	
							temp_intensity.push_back( input_interp->EvaluateAtContinuousIndex(offset_index) );
							temp_gradient_magnitude.push_back( gradient_magnitude_interp->EvaluateAtContinuousIndex(offset_index) );
						}
					}

					// Add values in negative direction
					// -d1
					offset_index[0]=index[0]-grad[0]*d1[i];
					offset_index[1]=index[1]-grad[1]*d1[i];
					offset_index[2]=index[2]-grad[2]*d1[i];

					if ( checkBounds(REGION, offset_index) )
					{
						temp_intensity.push_back( input_interp->EvaluateAtContinuousIndex(offset_index) );
						temp_gradient_magnitude.push_back( gradient_magnitude_interp->EvaluateAtContinuousIndex(offset_index) );

						// -d2
						offset_index[0]=index[0]-grad[0]*d2[i];
						offset_index[1]=index[1]-grad[1]*d2[i];
						offset_index[2]=index[2]-grad[2]*d2[i];

						if ( checkBounds(REGION, offset_index) )
						{	
							//std::cout << offset_index[0] << " " << offset_index[1] << " " << offset_index[2] << std::endl;
							//std::cout << input_interp->EvaluateAtContinuousIndex(offset_index) << std::endl;

							temp_intensity.push_back( input_interp->EvaluateAtContinuousIndex(offset_index) );
							temp_gradient_magnitude.push_back( gradient_magnitude_interp->EvaluateAtContinuousIndex(offset_index) );
						}
					}

					/*for (int j=0; j < temp_intensity.size(); j++)
					{
						std::cout << temp_intensity[j] << " ";
					}
					std::cout << std::endl;*/

					// Calculate distances to each parabola
					float distanceTS=0;
					float distanceSA=0;
					float distanceTA=0;
					float shortest_distance=0;
					VoxelType shortest_voxel=Unclassified;

					distanceTS=AverageTissueStoolDist(smax_iter.Get(),temp_intensity,temp_gradient_magnitude);
					distanceSA=AverageStoolAirDist(smax_iter.Get(),temp_intensity,temp_gradient_magnitude); 
					distanceTA=AverageTissueAirDist(temp_intensity,temp_gradient_magnitude);

					// Keep shortest one and compare to previous iterations
					if (distanceSA<=distanceTS && distanceSA<=distanceTA) {
						shortest_distance=distanceSA;
						shortest_voxel=StoolAir;
					
					} else if (distanceTS<=distanceTA && distanceTS<=distanceSA) {
						shortest_distance=distanceTS;
						shortest_voxel=TissueStool;
					
					} else {
						shortest_distance=distanceTA;
						shortest_voxel=TissueAir;
					}

					if ( i==0 ) // initial iteration
					{
						
						distance = shortest_distance;
						voxel = shortest_voxel;
					
					} else {
						
						if ( shortest_distance < distance )
						{
							distance = shortest_distance;
							voxel = shortest_voxel;
						}
					
					}
				}

			} else { // If no stool nearby, assume boundary is Tissue-Air
				voxel = TissueAir;
			}

			vmap_iter.Set( voxel );

			std::cout << count++ << std::endl;
		}

	}

	Write(vmap,"qr.nii");
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
1. Load image
2. Shift input so that air is 0 HU
3. Segment colon
4. Mask input with colon
5. Crop input and colon in XY plane
*********************************************************/
void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon)
{
	//----------------------------------------------
	// 1. Load image
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
	// 2. Shift input so that air is 0 HU
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

	// Set global REGION
	REGION = input_original->GetLargestPossibleRegion();

	//----------------------------------------------
	// 3. Segment colon
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
	// 4. Mask input with colon
	//----------------------------------------------
	typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( input_original );
	masker->SetInput2( colon );
	masker->Update();

	input = masker->GetOutput();

	Write(input_original,"input_original.nii");

	//----------------------------------------------
	// 5. Crop input and colon in XY plane
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

	// Set global REGION as new colon REGION
	REGION = extractRegion;

	typedef itk::ExtractImageFilter<ImageType,ImageType> ExtractImageFilterType;
	ExtractImageFilterType::Pointer extracter = ExtractImageFilterType::New();
	extracter->SetInput( input );
	extracter->SetExtractionRegion( REGION );
	extracter->Update();
	input = extracter->GetOutput();

	typedef itk::ExtractImageFilter<ByteImageType,ByteImageType> ExtractImageFilterByteType;
	ExtractImageFilterByteType::Pointer extracterByte = ExtractImageFilterByteType::New();
	extracterByte->SetInput( colon );
	extracterByte->SetExtractionRegion( REGION );
	extracterByte->Update();
	colon = extracterByte->GetOutput();

	Write(input,"input.nii");
	Write(colon,"colon.nii");
}

/*********************************************************
1. Calculate gradient magnitude
2. Allocate voxel map
3. Compute otsu threshold separating tissue and stool
4. Apply initial threshold rules
*********************************************************/
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon) 
{
	//----------------------------------------------
	// 1. Calculate gradient magnitude
	//----------------------------------------------
	//typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
	
	GradientMagnitudeImageFilterType::Pointer gradient_magnitude_filter = GradientMagnitudeImageFilterType::New();
	gradient_magnitude_filter->SetInput( input );

	gradient_magnitude_filter->SetSigma( input->GetSpacing()[0] );

	gradient_magnitude_filter->Update();
	gradient_magnitude = gradient_magnitude_filter->GetOutput();

	typedef itk::CastImageFilter<ImageType,FloatImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput(input);
	caster->Update();
	
	// No smoothing
	typedef itk::GradientMagnitudeImageFilter<FloatImageType,FloatImageType> GradientMagnitudeImageFilterNoSmoothingType;
	GradientMagnitudeImageFilterNoSmoothingType::Pointer grad_no_smooth_filter = GradientMagnitudeImageFilterNoSmoothingType::New();
	grad_no_smooth_filter->SetInput( caster->GetOutput() );
	grad_no_smooth_filter->Update();
	Write(grad_no_smooth_filter->GetOutput(), "grad_no_smooth.nii");

	// Sobel
	typedef itk::SobelEdgeDetectionImageFilter<FloatImageType,FloatImageType> SobelEdgeDetectionImageFilterType;
	SobelEdgeDetectionImageFilterType::Pointer sobel_filter = SobelEdgeDetectionImageFilterType::New();
	sobel_filter->SetInput( caster->GetOutput() );
	sobel_filter->Update();
	Write(sobel_filter->GetOutput(),"sobel.nii");


		

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

	Write(gradient_magnitude,"gradient_magnitude.nii");
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