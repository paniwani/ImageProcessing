#include "RemoveStool4.h"
#include "io.cxx"
#include "scattercorrection.cxx"
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
	//QuadraticRegression(input,colon,vmap);

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
	typedef itk::CannyEdgeDetectionImageFilterModified<ImageType,FloatImageType> CannyEdgeDetectionImageFilterType;
	CannyEdgeDetectionImageFilterType::Pointer cannyFilter = CannyEdgeDetectionImageFilterType::New();
	cannyFilter->SetInput(input);
	cannyFilter->SetLowerThreshold(0); // minimum gradient threshold
	cannyFilter->SetVariance(1);	

	typedef itk::CastImageFilter<FloatImageType, ByteImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput( cannyFilter->GetOutput() );
	caster->Update();
	
	ByteImageType::Pointer canny = caster->GetOutput();

	Dilate(canny, 1);

	//----------------------------------------------
	// 2. Get normalized gradient vector 
	//----------------------------------------------
	typedef itk::GradientImageFilter<ImageType> GradientImageFilterType;
	GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->Update();
	VectorImageType::Pointer gradient = gradientFilter->GetOutput();
	VectorIteratorType gradient_iter(gradient,region);

	for (gradient_iter.GoToBegin(); !gradient_iter.IsAtEnd(); ++gradient_iter)
	{
		VectorType grad = gradient_iter.Get();

		float norm = grad.GetNorm();
		
		if (norm > 0)
		{
			grad[0] /= norm;
			grad[1] /= norm;
			grad[2] /= norm;
			gradient_iter.Set(grad);
		}
	}

	//----------------------------------------------
	// 3. Interpolate input and gradient magnitude
	//----------------------------------------------
	typedef itk::BSplineInterpolateImageFunction<ImageType> BSplineInterpolateImageFunctionImageType;
	BSplineInterpolateImageFunctionImageType::Pointer input_interp = BSplineInterpolateImageFunctionImageType::New();
	input_interp->SetSplineOrder(3);
	input_interp->SetInputImage(input);

	typedef itk::BSplineInterpolateImageFunction<FloatImageType> BSplineInterpolateImageFunctionFloatImageType;
	BSplineInterpolateImageFunctionFloatImageType::Pointer gradient_magnitude_interp = BSplineInterpolateImageFunctionFloatImageType::New();
	gradient_magnitude_interp->SetSplineOrder(3);
	gradient_magnitude_interp->SetInputImage(gradient_magnitude);

	//----------------------------------------------
	// 4. Compute local stool maximum
	//----------------------------------------------
	ImageType::Pointer smax = ComputeNeighborhoodSmax(input,vmap,colon);
	IteratorType smax_iter(smax,region);

	//----------------------------------------------
	// 5. Run voxel edge classification
	//----------------------------------------------
	int count = 0;

	for (vmap_iter.GoToBegin(), smax_iter.GoToBegin(), gradient_iter.GoToBegin(), colon_iter.GoToBegin(), input_iter.GoToBegin(); !vmap_iter.IsAtEnd();
		++vmap_iter, ++smax_iter, ++gradient_iter, ++colon_iter, ++input_iter)
	{
		if (colon_iter.Get() == 255 && vmap_iter.Get() == Unclassified)
		{
			ImageType::IndexType index		=	input_iter.GetIndex();
			VectorType grad					=	gradient_iter.Get();


			VoxelType voxel = Unclassified;
			
			if (smax_iter.Get() > 0)
			{

				float distance = -1;

				for (int i=0; i<3; i++)
				{
					BSplineInterpolateImageFunctionFloatImageType::ContinuousIndexType offset_index;
					offset_index[0]=index[0];
					offset_index[1]=index[1];
					offset_index[2]=index[2];
					
					// Store intensity and gradient magnitude for 5 locations
					PixelType temp_intensity[5];
					float temp_gradient_magnitude[5];

					// Set target voxel
					temp_intensity[2]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
					temp_gradient_magnitude[2]=gradient_magnitude_interpolator->EvaluateAtContinousIndex(offset_index);

					//+d1
					offset_index[0]=index[0]+grad[0]*d1;
					offset_index[1]=index[1]+grad[1]*d1;
					offset_index[2]=index[2]+grad[2]*d1;
					temp_intensity[3]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
					temp_gradient_magnitude[3]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
					
					//-d1
					offset_index[0]=index[0]-grad[0]*d1;
					offset_index[1]=index[1]-grad[1]*d1;
					offset_index[2]=index[2]-grad[2]*d1;
					temp_intensity[1]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
					temp_gradient_magnitude[1]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
					
					//+d2
					offset_index[0]=index[0]+grad[0]*d2;
					offset_index[1]=index[1]+grad[1]*d2;
					offset_index[2]=index[2]+grad[2]*d2;
					temp_intensity[4]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
					temp_gradient_magnitude[4]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
					
					//-d2
					offset_index[0]=index[0]-grad[0]*d2;
					offset_index[1]=index[1]-grad[1]*d2;
					offset_index[2]=index[2]-grad[2]*d2;
					temp_intensity[0]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
					temp_gradient_magnitude[0]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);

					float distanceTS=0;
					float distanceSA=0;
					float distanceTA=0;

					distanceTS=AverageTissueStoolDist(smax, temp_intensity,temp_gradient_magnitude);
					distanceSA=AverageStoolAirDist(smax, temp_intensity,temp_gradient_magnitude); 
					distanceTA=AverageTissueAirDist(temp_intensity,temp_gradient_magnitude);

					
				}

				
			} else { // If no stool nearby, assume boundary is Tissue-Air
				voxel = TissueAir;
			}



		}

	}

	


    for (voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), gradient_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), input_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(); 
		!voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !gradient_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++voxel_type_iter, ++input_smax_iter, ++gradient_iter, ++chamfer_colon_iter, ++input_iter, ++gradient_magnitude_iter) {
		
		if (chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get()==Unclassified) {
			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
			CovariantVectorType grad = gradient_iter.Get();	

			float temp_threshold=-1;
			VoxelType temp_type = Unclassified;

			if (!Modified)
			{
				VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);
				
				VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);

				VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);
			} else {

				if ( input_smax_iter.Get() > 0 )
				{
					VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);
					
					VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);

					VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);
				
				} else {
					temp_type = TissueAir;
				}

			}

			voxel_type_iter.Set(temp_type);

			vefile << voxel_index[0] << "\t" << 511-voxel_index[1] << "\t" << voxel_index[2] << "\t";
			vefile << grad[0] << "\t" << grad[1] << "\t" << grad[2] << "\t";
			vefile << input_iter.Get() << "\t" << gradient_magnitude_iter.Get() << "\t" << input_smax_iter.Get() << "\t" << VoxelTypeToString(temp_type) << "\n";
		}

		if (++count % 100000 == 0)
			std::cout << "Voxel edge count: " << count << std::endl;
    }









	

	


	
	


}

ImageType::Pointer ComputeNeighborhoodSmax(ImageType::Pointer &input, VoxelImageType::Pointer &v, ByteIteratorType &mask_iter)
{
	ImageType::Pointer smax = AllocateNewImage(input->GetLargestPossibleRegion());
	IteratorType smax_iter(smax,input->GetLargestPossibleRegion());
	IteratorType input_iter(input,input->GetLargestPossibleRegion());

	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorVoxelType;
	NeighborhoodIteratorVoxelType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 2;
	radius[1] = 2;

	NeighborhoodIteratorVoxelType nit(radius, v, v->GetLargestPossibleRegion());

	VoxelTypeImage::SizeType size = nit.GetSize();
	int n = size[0]*size[1]*size[2] - 1;

	nit.GoToBegin();
	mask_iter.GoToBegin();
	smax_iter.GoToBegin();
	input_iter.GoToBegin();

	PixelType max=0;

	while (!nit.IsAtEnd())
	{
		if (mask_iter.Get() == 255 && nit.GetCenterPixel()==Unclassified)
		{
			max = 0;

			for (int i=0; i<n; i++)
			{
			
				VoxelTypeImage::IndexType idx = nit.GetIndex(i);
				
				if ( region.IsInside( idx ) )
				{
					PixelType val = input->GetPixel(idx);
					
					if (nit.GetPixel(i) == Stool)
					{

						if ( val > max )
						{
							max = val;
						}	
					}				
				}
			}

			smax_iter.Set( max );
		
		} else {
			smax_iter.Set( 0 );
		}

		++nit;
		++mask_iter;
		++smax_iter;
		++input_iter;
	}

	return smax;
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

	typedef itk::MinimumMaximumImageCalculator< ShortImageType > MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer calc = MinimumMaximumImageCalculatorType::New();
	calc->SetImage( input_s );
	calc->ComputeMinimum();
	ShortImageType::PixelType min = calc->GetMinimum();

	for (input_s_iter.GoToBegin(), input_original_iter.GoToBegin(); !input_s_iter.IsAtEnd(); ++input_s_iter, ++input_original_iter)
	{
		input_original_iter.Set( input_s_iter.Get() + abs( min ) );
	}

	input_s.~SmartPointer();

	// Set global region
	region = input_original->GetLargestPossibleRegion();

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

	//----------------------------------------------
	// 2. Allocate voxel map
	//----------------------------------------------
	vmap->SetRegions( input->GetLargestPossibleRegion() );
	vmap->SetSpacing( input->GetSpacing() );
	vmap->SetDirection( input->GetDirection() );
	vmap->CopyInformation( input );
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);
	
	IteratorType input_iter(input,region);
	FloatIteratorType gradient_magnitude_iter(gradient_magnitude,region);
	VoxelIteratorType vmap_iter(vmap,region);
	ByteIteratorType colon_iter(colon,region);

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
	IteratorType input_iter(input,region);
	FloatIteratorType gradient_magnitude_iter(gradient_magnitude,region);
	VoxelIteratorType vmap_iter(vmap,region);
	ByteIteratorType colon_iter(colon,region);

	tissue_stool_threshold -= 1024;

	for ( input_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gradient_magnitude_iter, ++vmap_iter, ++colon_iter)
	{
		if (colon_iter.Get() == 255)
		{
			VoxelType voxel;

			short I = input_iter.Get() - 1024;
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