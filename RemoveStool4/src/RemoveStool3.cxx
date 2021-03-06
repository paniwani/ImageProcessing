#include "RemoveStool3.h"
#include "HessianFunctions.cxx"
#include "ScatterCorrection.cxx"

int main(int argc, char * argv[])
{
	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}

	ImageType::Pointer input;
	ByteImageType::Pointer colon;

	Setup(argv[1],input,colon);


}

/*********************************************************
1. Load image
2. Segment colon
3. Mask input with colon
4. Rescale input so that air is 0 HU
*********************************************************/
void Setup(std::string dataset, ImageType::Pointer input, ByteImageType::Pointer colon)
{
	//----------------------------------------------
	// 1. Load image
	//----------------------------------------------
	std::vector<std::string> dataset_ar = explode( "\\", dataset );
	std::string dsname = dataset_ar[ dataset_ar.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	if (truncateOn)
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
	if (!truncateOn)
	{
		input_s = ReadDicom( dataset );
	} else {
		input_s = ReadDicom( dataset, truncate_ar[0], truncate_ar[1] );
	}

	// Scale and cast to 16-bit unsigned short
	typedef itk::MinimumMaximumImageCalculator< ShortImageType > MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage( input_s );
	minMaxCalc->ComputeMinimum();
	ShortImageType::PixelType min = minMaxCalc->GetMinimum();

	IteratorShortType
	



	//----------------------------------------------
	// 2. Segment colon
	//----------------------------------------------
	ColonSegmenationFilterType::Pointer colon_segmenter = ColonSegmenationFilterType::New();
	colon_segmenter->SetInput( input );

	if (truncateOn)
		colon_segmenter->SetRemoveBoneLung( false );

	colon_segmenter->SetBackgroundValue( 0 );
	colon_segmenter->SetForegroundValue( 1 );
	colon_segmenter->Update();

	ByteImageType::Pointer chamfer_colon = colon_segmenter->GetOutput();

	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);




}









int main2(int argc, char * argv[])
{
	//-------------------------------------------BEGIN SETUP-------------------------------------------------------

	// Get arguments
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory Modified";
		system("pause");
		return EXIT_FAILURE;
	}
	
	std::cerr<<"Started"<<std::endl;

	// Get dataset name
	std::string dataset = argv[1];
	std::vector<std::string> datasetArr = explode( "\\", dataset );
	std::string dsname = datasetArr[ datasetArr.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	if (truncateOn)
	{
		std::cout << "Truncating data" << std::endl;
		std::cout << "Slices: " << truncate_ar[0] << " to " << truncate_ar[1] << std::endl;
		std::stringstream ss;
		ss << dsname << "_" << truncate_ar[0] << "_" << truncate_ar[1];
		dsname = ss.str();
	}

	note = dsname;

	// Use parameter changes
	if (atoi(argv[2]) == 1)
	{
		Modified = true;
	}

	std::cout << "Modified: " << Modified << std::endl;

	// Read and write dicom input
	ImageType::Pointer input = ReadDicom( dataset );

	WriteITK(input,"input.nii");

	/*

	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	input = ResampleImage(input, iso_spacing);
	WriteITK(input,"input_isotropic.nii");

	*/
	
	// Get region and bounds
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);

	IteratorTypeFloat4WithIndex input_iter(input,fullRegion);

	// Segment colon
	ColonSegmenationFilterType::Pointer colon_segmenter = ColonSegmenationFilterType::New();
	colon_segmenter->SetInput( input );

	if (truncateOn)
		colon_segmenter->SetRemoveBoneLung( false );

	colon_segmenter->SetBackgroundValue( 0 );
	colon_segmenter->SetForegroundValue( 1 );
	colon_segmenter->Update();

	ByteImageType::Pointer chamfer_colon = colon_segmenter->GetOutput();

	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);

	WriteITK(chamfer_colon, "colon.nii");

	// Get min, max, and mean intensity within colon
	float input_min = itk::NumericTraits<PixelType>::max();
	float input_max = itk::NumericTraits<PixelType>::NonpositiveMin();
	float input_mean = 0;
	unsigned int input_count = 0;

	for (input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++chamfer_colon_iter)
	{
		if (chamfer_colon_iter.Get()==1)
		{
			++input_count;

			PixelType I = input_iter.Get();

			if ( I < input_min )
				input_min = I;

			if ( I > input_max )
				input_max = I;

			input_mean += I;

		}
	}

	input_mean /= input_count;

	std::cout << std::endl << std::endl;

	std::cout << "Input minimum: " << input_min << std::endl;
	std::cout << "Input maximum: " << input_max << std::endl;
	std::cout << "Input mean: " << input_mean << std::endl;

	std::cout << std::endl << std::endl;


	RemoveStool5(input,chamfer_colon);

	/*

	// Scatter correction
	ImageType::Pointer input_scatter = ScatterCorrection(input, chamfer_colon);

	// Write change only image
	ImageType::Pointer scatter_change = AllocateNewImage(fullRegion);
	
	IteratorTypeFloat4WithIndex scatter_change_iter(scatter_change,fullRegion);
	IteratorTypeFloat4WithIndex input_scatter_iter(input_scatter,fullRegion);

	for (input_iter.GoToBegin(), input_scatter_iter.GoToBegin(), scatter_change_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++input_scatter_iter, ++scatter_change_iter)
	{
		scatter_change_iter.Set( input_iter.Get() - input_scatter_iter.Get() );		
	}	

	WriteITK(scatter_change, "scatter_change.nii");

	scatter_change.~SmartPointer();

	input = input_scatter;

	WriteITK(input, "input_scatter.nii");

	input_scatter.~SmartPointer();

	*/

	////////////////////////////////TESTING AREA
	//double gamma[6] = {0.25,.5,1,2,4,10};
	//for (int i=0; i<6; i++)

	//double sigma[2] = {0.56,1.4};
	//double gamma = 1;
	//
	//Imag*/eType::Pointer hessian = SatoModifiedResponse(input,chamfer_colon,sigma[i],gamma);


	//ImageType::Pointer test = HessianResponse(input,ALPHA,BETA,GAMMA,ETA);
	//Objectness(input);
	

	//GetHessian(input);
	//UnderstandHessian(input);
	
	/*
	GetHessian(input);

	double alpha = 1;
	double gamma = 0.5;

	;*/
	//ProcessHessian(hessian, voxel_type);



	//ImageType::Pointer sh = SatoHessianEdgeEnhancingDiffusion(input);
	//WriteITK(sh,"sh.nii");


	/*
	// Sato params
	double alpha = 0.25;
	double gamma[5] = {0.5,1,2,5,10};

	for (int i=0; i<5; i++)
	{
		ImageType::Pointer sato_hessian = SatoResponse(input, alpha, gamma[i]);
	}
	*/
	////////////////////////////////////////////////////////

	// Creates an image to store the type of voxel
	VoxelTypeImage::Pointer voxel_type = VoxelTypeImage::New();
	voxel_type->SetRegions(fullRegion);
	voxel_type->Allocate();
	voxel_type->FillBuffer(Unclassified);
	IteratorTypeVoxelType voxel_type_iter(voxel_type,fullRegion);

	// Get gradient vectors
    GradientFilterType::Pointer gradient_filter = GradientFilterType::New();
    gradient_filter->SetInput(input);
    gradient_filter->Update();
    ImageVectorType::Pointer gradient = gradient_filter->GetOutput();
	gradient_filter.~SmartPointer();
	IteratorImageVectorType gradient_iter(gradient,fullRegion);

	// Calculate gradient magnitude and normalize gradient
	ImageType::Pointer gradient_magnitude = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex gradient_magnitude_iter(gradient_magnitude,fullRegion);

	for (gradient_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin();
		 !gradient_iter.IsAtEnd() && !gradient_magnitude_iter.IsAtEnd();
		 ++gradient_iter, ++gradient_magnitude_iter)
	{
		CovariantVectorType grad = gradient_iter.Get();
		float norm = grad.GetNorm();
		
		if (norm > 0)
		{
			grad[0] /= norm;
			grad[1] /= norm;
			grad[2] /= norm;
			gradient_iter.Set(grad);
		}

		gradient_magnitude_iter.Set( norm );
	}
	
	// Write out the gradient magnitude for the input image
	WriteITK(gradient_magnitude, "gradient_1x.nii");



	// Make mask of G < 0.8*I

	ByteImageType::Pointer gmask = AllocateNewByteImage(fullRegion);
	IteratorTypeByteWithIndex gmask_iter(gmask,fullRegion);

	for (input_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), gmask_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gradient_magnitude_iter, ++chamfer_colon_iter, ++gmask_iter)
	{
		if (chamfer_colon_iter.Get() == 1)
		{
			if ( gradient_magnitude_iter.Get() < (0.8 * input_iter.Get()) )
			{
				gmask_iter.Set( 255 );
			} else {
				gmask_iter.Set( 0 );
			}
		}
	}

	WriteITK( gmask, "gmask.nii");

	// Create temporary image and iter
	ImageType::Pointer temp = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex temp_iter(temp,fullRegion);

	// Compute an initial hard threshold for the voxel type
	SingleMaterialClassification(input, gradient_magnitude, voxel_type, chamfer_colon);

	WriteITK(voxel_type, "voxel_type_start.nii");

	// Clean up stool voxels
	CleanIsolatedStool( voxel_type );

	WriteITK(voxel_type,"voxel_type_stool_clean.nii");

	/*
	// Hessian Analysis

	double alpha = 1;
	double gamma = 3;

	ImageType::Pointer hessian = SatoResponse(input, chamfer_colon, alpha, gamma);

	WriteITK(hessian,"hessian.nii");

	ProcessHessian( hessian, voxel_type );

	//double contrast[7] = {5,10,20,30,50,100,1000};
	/*double contrast = 15.0;
	double sigma[5] = {0.5,1,2,4,10};
	
	for (int i=0; i<5; i++)
		DiffusionTest( hessian , contrast, sigma[i]);
	
	//ProcessHessian(hessian,voxel_type);

	*/		
	
	if (Modified)
	{
		// Close any voxel components surrounded by stool
		HeuristicClosing(voxel_type,chamfer_colon);
		WriteITK(voxel_type,"voxel_type_stool_closed_heuristic.nii");
		
		//ImageType::Pointer hessian = SatoResponse(input, alpha, gamma);
	

		//// Run hessian and update voxel classification
		//ImageType::Pointer hessian = HessianResponse(input);
	
		//WriteITK(hessian,"hessian.nii");

		//// Rescale hessian to [0,1]
		//RescaleIntensityFilterType::Pointer rescaler = RescaleIntensityFilterType::New();
		//rescaler->SetOutputMaximum(1);
		//rescaler->SetOutputMinimum(0);
		//rescaler->SetInput(hessian);
		//rescaler->Update();
		//hessian = rescaler->GetOutput();

		//WriteITK(hessian,"hessian_rescaled.nii");


		//EnhanceVoxelType(input,hessian,voxel_type,chamfer_colon);

		//WriteITK(voxel_type,"voxel_type_hessian_enhanced.nii");
	}

	//-------------------------------------------END SETUP-------------------------------------------------------

	//-------------------------------------------BEGIN QR--------------------------------------------------------
	
    // Get cubic bspline image interpolation
    InterpolationType::Pointer input_interpolator = InterpolationType::New();
    input_interpolator->SetSplineOrder(3);
    input_interpolator->SetInputImage(input);

	for (gradient_magnitude_iter.GoToBegin(); !gradient_magnitude_iter.IsAtEnd(); ++gradient_magnitude_iter)
	{
		gradient_magnitude_iter.Set(2*gradient_magnitude_iter.Get());
	}

	WriteITK(gradient_magnitude, "gradient_magnitude_2x.nii");
	
	InterpolationType::Pointer gradient_magnitude_interpolator = InterpolationType::New();
	gradient_magnitude_interpolator->SetSplineOrder(3);
	gradient_magnitude_interpolator->SetInputImage(gradient_magnitude);

	// Storage for stool maximum
	ImageType::Pointer input_smax;

	if (!Modified)
	{
		input_smax = AllocateNewImage(fullRegion);
		input_smax->FillBuffer(0);
	} else {
		input_smax = ComputeNeighborhoodSmax(input, voxel_type, chamfer_colon_iter, startIndex, endIndex);
	}

	// Load smax
	//ReadITK(input_smax, "C:/ImageData/mr10_092_13p.i0344_100-105_smax2.nii");
	
	IteratorTypeFloat4WithIndex input_smax_iter(input_smax,fullRegion);

	std::cerr << "Running voxel edge classification..." << std::endl;

	// Write voxel edge to text file
	std::ofstream vefile;
	vefile.open("ve.txt");
	vefile << "Xindex\tYindex\tZindex\tXgrad\tYgrad\tZgrad\tIntensity\tGradient\tSmax\tType\n";

	//ReadITK(voxel_type, "Modified_6_voxel_edge_class.nii");
	voxel_type_iter = IteratorTypeVoxelType(voxel_type,fullRegion);

	int count = 0;
	//Runs edge classification with 3 different increments inside the colon only
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

	WriteITK(input_smax,"smax.nii");

	vefile.close();

	// Deletes the Interpolators
	input_interpolator.~SmartPointer();

	WriteITK(voxel_type,"voxel_edge_class.nii");

 	std::cerr<<"Done Edge"<<std::endl;

	// Remove TA and TS that are far from tissue
	// Create tissue mask
	ByteImageType::Pointer chamfer_tissue = AllocateNewByteImage(fullRegion);
	chamfer_tissue->FillBuffer(0);
	IteratorTypeByteWithIndex chamfer_tissue_iter(chamfer_tissue,fullRegion);

	for (chamfer_tissue_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !chamfer_tissue_iter.IsAtEnd(); ++chamfer_tissue_iter, ++voxel_type_iter)
	{
		if (voxel_type_iter.Get() == Tissue)
		{
			chamfer_tissue_iter.Set(1); // distance to tissue
		}
	}

	WriteITK(chamfer_tissue, "tissue_involvement.nii");

	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_tissue);
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_tissue=chamfer_filter->GetOutput();
	chamfer_tissue_iter = IteratorTypeByteWithIndex(chamfer_tissue,fullRegion);
	chamfer_filter.~SmartPointer();

	WriteITK(chamfer_tissue, "chamfer_tissue_distance.nii");

	for (chamfer_tissue_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !chamfer_tissue_iter.IsAtEnd(); ++chamfer_tissue_iter, ++voxel_type_iter)
	{
		if (voxel_type_iter.Get() == TissueAir || voxel_type_iter.Get() == TissueStool)
		{
			if (chamfer_tissue_iter.Get() > 5) 
			{
				voxel_type_iter.Set(StoolAir);
			}
		}
	}

	WriteITK(voxel_type, "voxel_type_tissue_fix_5.nii");

	//SubtractStool( input, voxel_type, chamfer_colon );


	//Optimize SA transitions using gradient information
	//std::cerr<<"Running voxel edge optimization"<<std::endl;
	
	//OptimizeVoxelEdge(input, voxel_type, gradient);
	
	//WriteITK(voxel_type, "voxel_edge_optimized_noz.nii");

	//// Find thin stool regions [Carston 2.3.6]
	//ByteImageType::Pointer chamfer_from_stool_involvement=AllocateNewByteImage(fullRegion);
	//IteratorTypeByteWithIndex chamfer_from_stool_involvement_iter(chamfer_from_stool_involvement,fullRegion);

	//for(voxel_type_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(); 
	//	!voxel_type_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(); 
	//	++voxel_type_iter, ++chamfer_from_stool_involvement_iter) 
 //   {
	//	if (voxel_type_iter.Get()==Stool || voxel_type_iter.Get()==TissueStool || voxel_type_iter.Get()==StoolAir) {
	//		chamfer_from_stool_involvement_iter.Set(0);     //stool
	//	} else {
	//		chamfer_from_stool_involvement_iter.Set(1);     //non stool is the object
	//	}
	//}

	////WriteITK(chamfer_from_stool_involvement,"non_stool_mask.nii");

	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_from_stool_involvement);
	//int weights[3]={3,4,5};	//3d distance weight recommended by julian
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(true);
	//chamfer_filter->Update();
	//chamfer_from_stool_involvement=chamfer_filter->GetOutput();
	//chamfer_from_stool_involvement_iter = IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	//chamfer_filter.~SmartPointer();

	////WriteITK(chamfer_from_stool_involvement, "non_stool_distance.nii");

	//// Find thin stool voxels < 8 and no neighbors > 8 (which were not already classified as stool)
	//int chamfercutoff = 8;
	//int thinStoolCount = 0;

	//for(voxel_type_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
	//	!voxel_type_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(); 
	//	++voxel_type_iter, ++chamfer_from_stool_involvement_iter) 
 //   {
	//	//if (voxel_type_iter.Get() != Stool && voxel_type_iter.Get() != TissueStool) // do not turn these into thin stool
	//	//{
	//		if (chamfer_from_stool_involvement_iter.Get() < chamfercutoff && chamfer_from_stool_involvement_iter.Get() > 0)
	//		{
	//			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();

	//			bool thinNeighbors = true;

	//			for(int i=-1;i<=1; i++) {
	//				if (thinNeighbors && index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
	//					for (int j=-1;j<=1 && thinNeighbors; j++) {
	//						if (thinNeighbors && index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
	//							for (int k=-1;k<=1; k++) {
	//								if (thinNeighbors && index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
	//									
	//									ImageType::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]+k};
	//									
	//									if (chamfer_from_stool_involvement->GetPixel(neighbor_index) > chamfercutoff) {
	//										thinNeighbors = false;
	//									}
	//								
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}

	//			if (thinNeighbors)
	//			{
	//				voxel_type_iter.Set(ThinStool);
	//				thinStoolCount++;
	//			}
	//		}
	//	//}
	//}

	//std::cout << "Number of thin stool voxels: " << thinStoolCount << std::endl;

	//WriteITK(voxel_type,"voxel_type_thinstool.nii");

	//// Apply voting scheme to voxel edge
	////VoteVoxels(voxel_type, chamfer_colon);
	//
	////voxel_type_iter=IteratorTypeVoxelType(voxel_type,fullRegion);
	//
	////WriteITK(voxel_type, "voxel_edge_voting.nii");

	//// Find thin stool distance to air
	//ByteImageType::Pointer chamfer_air=AllocateNewByteImage( fullRegion );
 //   IteratorTypeByteWithIndex chamfer_air_iter(chamfer_air,fullRegion);

 //   for(voxel_type_iter.GoToBegin(), chamfer_air_iter.GoToBegin();
 //       !voxel_type_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd();
 //       ++voxel_type_iter, ++chamfer_air_iter)
 //   {
	//	if (voxel_type_iter.Get()==Air) {
 //           chamfer_air_iter.Set(1);					
 //       } else {
 //           chamfer_air_iter.Set(0);	
	//	}
 //   }

	////compute chamfer for air
	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_air);
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(true);
	//chamfer_filter->Update();
	//chamfer_air=chamfer_filter->GetOutput();
	//chamfer_air_iter= IteratorTypeByteWithIndex(chamfer_air,fullRegion);
	//chamfer_filter.~SmartPointer();
	////WriteITK(chamfer_air,"chamfer_to_air.nii");
	//
	//// Convert stool and air chamfers to float and normalize
	//ImageType::Pointer ds = AllocateNewImage(fullRegion);
	//IteratorTypeFloat4WithIndex ds_iter(ds,fullRegion);

	//ImageType::Pointer da = AllocateNewImage(fullRegion);
	//IteratorTypeFloat4WithIndex da_iter(da,fullRegion);

	//for(chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), da_iter.GoToBegin(), ds_iter.GoToBegin();
	//	!chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(), !da_iter.IsAtEnd(), !ds_iter.IsAtEnd();
 //       ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter, ++da_iter, ++ds_iter)
 //   {
	//	if ( chamfer_from_stool_involvement_iter.Get() >= 8)
	//	{
	//		ds_iter.Set(1);
	//	} else {
	//		ds_iter.Set( (float) chamfer_from_stool_involvement_iter.Get() / 8 );
	//	}

	//	if ( chamfer_air_iter.Get() >= 8)
	//	{
	//		da_iter.Set(1);
	//	} else {
	//		da_iter.Set( (float) chamfer_air_iter.Get() / 8 );
	//	}
	//}

	//da_iter = IteratorTypeFloat4WithIndex(da,fullRegion);
	//ds_iter = IteratorTypeFloat4WithIndex(ds,fullRegion);

	////WriteITK(da,"normalized_air.nii");
	////WriteITK(ds,"normalized_stool.nii");

	//-------------------------------------------------COMPUTE PARTIAL VOLUMES---------------------------------------------------------------
	
	// Declares a structure to store the air/tissue/stool partial volumes for a voxel in the image
	ImageVectorType::Pointer partialVector = ImageVectorType::New();
    partialVector->SetRegions(fullRegion);
    partialVector->Allocate();
	partialVector->SetSpacing(input->GetSpacing());
    IteratorImageVectorType partialVector_iter(partialVector,fullRegion);

	int total_vertices=0;

	//std::ofstream file;
	//file.open("pv.txt");
	//file << "Index\tIntensity\tGradient\tClass\tPa\tPt\tPs\n";

	for(partialVector_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin()/*, da_iter.GoToBegin(), ds_iter.GoToBegin()*/;
		!partialVector_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !gradient_magnitude_iter.IsAtEnd()/* && !da_iter.IsAtEnd() && !ds_iter.IsAtEnd()*/ ; 
        ++partialVector_iter, ++input_iter, ++chamfer_colon_iter, ++voxel_type_iter, ++input_smax_iter, ++gradient_magnitude_iter/*, ++da_iter, ++ds_iter*/) 
	{
		if (chamfer_colon_iter.Get() == 1)
		{
			total_vertices++;

			ImageType::IndexType idx = input_iter.GetIndex();
			//file << "(" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t";
			//file << input_iter.Get() << "\t" << gradient_magnitude_iter.Get() << "\t";
			
			float value[3];	// air/tissue/stool
			
			switch (voxel_type_iter.Get()) {
				case Air:
					value[0]=1;					
					value[1]=0;
					value[2]=0;

					//file << "Air";
					break;
				case Tissue:
					value[0]=0;
					value[1]=1;					
					value[2]=0;

					//file << "Tissue";
					break;
				case Stool:
					value[0]=0;
					value[1]=0;
					value[2]=1;		

					//file << "Stool";
					break;
				case TissueAir:
					value[1]=1+(input_iter.Get()/1000);
					
					if (value[1] <= 0) { value[1] = 0; }
					if (value[1] >= 1) { value[1] = 1; }

					value[0]=1-value[1];
					value[2]=0;

					//file << "TissueAir";
					break;
				case TissueStool:
					value[0]=0;
					value[2]=input_iter.Get()/input_smax_iter.Get();
				
					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[1]=1-value[2];

					//file << "TissueStool";
					break; 
				case StoolAir:
					value[2]=(input_iter.Get()+1000)/(input_smax_iter.Get()+1000);
					
					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[0]=1-value[2];
					value[1]=0;

					//file << "StoolAir";
					break;
				case ThinStool:
					//value[0]=0.5*(1-vnl_erf((da_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));

					//if (Modified)
					//{
					//	value[0]=0.5*(1-vnl_erf((da_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2)))); //check eq
					//} else {
					//	value[0]=0.5*(1+vnl_erf((da_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
					//}
					//*/
					//
					//if (value[0] <= 0) { value[0] = 0; }
					//if (value[0] >= 1) { value[0] = 1; }

					//value[2]=0.5*(1+vnl_erf((ds_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));

					//if (value[2] <= 0) { value[2] = 0; }
					//if (value[2] >= 1) { value[2] = 1; }

					//value[1]=1-value[0]-value[2];

					//if (value[1] <= 0) { value[1] = 0; }
					//if (value[1] >= 1) { value[1] = 1; }

					//file << "ThinStool";
					break;
			}   

			//file << "\t" << value[0] << "\t" << value[1] << "\t" << value[2] << "\n";

			CovariantVectorType data(value);	//sets the partial to be stored in the actual structure
			partialVector_iter.Set(data);		//stores the partial in the structure
		}
    }

	std::cout << "Partials initialized" << std::endl;

	partialVector_iter = IteratorImageVectorType(partialVector,fullRegion);

	//file.close();

	WritePartialImages(partialVector, chamfer_colon, "QR");

	/*
	// Delete low pt
	for (partialVector_iter.GoToBegin(); !partialVector_iter.IsAtEnd(); ++partialVector_iter)
	{
		CovariantVectorType p = partialVector_iter.Get();
		if (p[1] < 0.3)
		{
			p[1] = 0;
			partialVector_iter.Set(p);
		}
	}

	WritePartialImages(partialVector, chamfer_colon, "Threshold");
	*/

	/*
	// Remove unconnected pt > 0 voxels
	ByteImageType::Pointer pt_mask = AllocateNewByteImage(fullRegion);
	IteratorTypeByteWithIndex pt_mask_iter(pt_mask,fullRegion);
	
	pt_mask_iter.GoToBegin();
	partialVector_iter.GoToBegin(); 
	chamfer_colon_iter.GoToBegin();

	while (!pt_mask_iter.IsAtEnd())
	{
		// Create mask of connected pt
		if ( chamfer_colon_iter.Get() == 1 && partialVector_iter.Get()[1] > 0 )
		{
			pt_mask_iter.Set( 1 );
		} else {
			pt_mask_iter.Set( 0 );
		}

		++pt_mask_iter;
		++partialVector_iter;
		++chamfer_colon_iter;
	}

	//WriteITK(pt_mask,"pt_mask.nii");

	// Run cc on pt mask
	ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput(pt_mask);
	ccFilter->Update();
	cc=ccFilter->GetOutput();

	//WriteITK(cc, "pt_mask_cc.nii");

	converter = LabelImageToShapeLabelMapFilterType::New();
	converter->SetInput(cc);
	converter->Update();
	labelMap = converter->GetOutput();
	
	converter.~SmartPointer();
	ccFilter.~SmartPointer();

	labelVector = labelMap->GetLabelObjects();
	sort(labelVector.begin(), labelVector.end(), compareSizeOnBorder);
	int tissueLabel = labelVector[0]->GetLabel();

	labelVector.clear();
	labelMap.~SmartPointer();

	cc_iter=IteratorTypeIntWithIndex(cc,fullRegion );
	
	// Set to air partial
	for (partialVector_iter.GoToBegin(), cc_iter.GoToBegin();
        !partialVector_iter.IsAtEnd() && !cc_iter.IsAtEnd();  
        ++partialVector_iter, ++cc_iter) 
	{
		if (cc_iter.Get() != tissueLabel)
		{
			CovariantVectorType p;
			p[0] = 1;
			p[1] = 0;
			p[2] = 0;
			partialVector_iter.Set( p );
		}
	}

	cc.~SmartPointer();

	WritePartialImages(partialVector, chamfer_colon, "QR_CC");
	*/

	// Smooth partial vector
	SmoothPartialVector(partialVector,chamfer_colon,startIndex,endIndex);

	//WritePartialImages(partialVector, chamfer_colon, "QR_Smooth");

	// Subtract stool via QR
	ImageType::Pointer output = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex output_iter(output,fullRegion);

	input_iter.GoToBegin();
	output_iter.GoToBegin();
	partialVector_iter.GoToBegin();
	chamfer_colon_iter.GoToBegin();
	voxel_type_iter.GoToBegin();

	while(!input_iter.IsAtEnd())
	{
		if (chamfer_colon_iter.Get() == 1)
		{

			float pt = partialVector_iter.Get()[1];

			if ( pt < 0.05 )
			{
				// Set to air partial
				CovariantVectorType p;
				p[0] = 1;
				p[1] = 0;
				p[2] = 0;
				partialVector_iter.Set( p );

				output_iter.Set(-1024);
			} else if ( pt >= 0.05 && pt <= 0.95) {
				output_iter.Set( pt*(input_iter.Get()+1000)-1000 ); 
			} else {
				output_iter.Set(input_iter.Get());
			}
		} else {
			output_iter.Set(input_iter.Get());
		}

		++input_iter;
		++output_iter;
		++partialVector_iter;
		++chamfer_colon_iter;
		++voxel_type_iter;
	}

	WriteITK(output,"outputQR.nii");

	/*

	// Remove stool based on tissue partials
	for(partialVector_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin();
		!partialVector_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !gradient_magnitude_iter.IsAtEnd();
        ++partialVector_iter, ++input_iter, ++chamfer_colon_iter, ++gradient_magnitude_iter)
	{
		if (chamfer_colon_iter.Get() == 1 ) 
		{
			CovariantVectorType partial = partialVector_iter.Get();

			
			// Air
			if (partial[0] > 0.9 && gradient_magnitude_iter.Get() < 250)
			{
				input_iter.Set(-1025);
			}

			// Stool
			if (partial[2] > 0.9 && gradient_magnitude_iter.Get() < 0.8*input_iter.Get())
			{
				input_iter.Set(-1025);
			}
			
			//input_iter.Set( partial[1]*(input_iter.Get()+1000)-1000 );
			
		}	
	}

	WriteITK(input_temp,"output.nii");
	*/

	//-------------------------------------------END QR-----------------------------------------------------

	//-------------------------------------------BEGIN EM--------------------------------------------------------

 //   int counter=0;							//used to total certain values (description given when set)
 //   double mean[3]={0, 0, 0};				//stores the mean of the air/tissue/stool classes
 //   double sum[3]={0, 0, 0};				//stores the total # of partials of the three classes
 //   double variance[3]={0, 0, 0};			//stores the variance of the air/tissue/stool classes
	//float weight[3]={0,0,0};				//stores the weights of the air/tissue/stool classes

	//// Computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
	//for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
	//	!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
	//	++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//{
	//	if (chamfer_colon_iter.Get()==1) {		
	//		CovariantVectorType partial = partialVector_iter.Get();
	//		for (int i=0; i<3; i++)
	//		{
	//			sum[i] += partial[i];
	//			mean[i] += partial[i]*input_iter.Get();
	//		}
	//	}
	//}
	//for (int i=0;i<3;i++) { mean[i]=mean[i]/sum[i]; } 

	//// Compute variance and weights
	//for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
	//	!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
	//	++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//{
	//	if (chamfer_colon_iter.Get()==1) {		
	//		CovariantVectorType partial = partialVector_iter.Get();
	//		for (int i=0; i<3; i++)
	//		{
	//			variance[i] += partial[i]*vnl_math_sqr(input_iter.Get()-mean[i]);
	//		}
	//	}
	//}

	//double sum_all = sum[0]+sum[1]+sum[2];

	//for (int i=0;i<3;i++) 
	//{ 
	//	variance[i]=variance[i]/sum[i];
	//	weight[i]=sum[i]/sum_all;
	//} 

	////outputs the initial mean and variance of each class
	//std::ofstream em;
	//em.open("em.txt");
	//em<<"EM\tMean0\tMean1\tMean2\tVar0\tVar1\tVar2\tWeight0\tWeight1\tWeight2\n";

	//std::ostringstream ss;
	//ss<<"0\t"<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<variance[0]<<"\t"<<variance[1]<<"\t"<<variance[2]<<"\t"<<weight[0]<<"\t"<<weight[1]<<"\t"<<weight[2]<<"\n";
	//em<<ss.str();

	//std::cerr<<std::endl;
	//std::cerr<<"EM0"<<std::endl;
	//std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	//std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	//std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;
	//
	//// Computes iterations of the Maximization Algorithm
 //   for (int emNum=0;emNum<20;emNum++) // used 20 in original test case
	//{
	//	//gives the temporary storage of the variables corresponding to sum, variance, and mean for each class on the i+1 iteration
	//	double sum_temp[3]={0,0,0};
 //       double variance_temp[3]={0,0,0};
	//	double mean_temp[3]={0,0,0};

	//	for (input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
	//		!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
	//		++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//	{
	//		if (chamfer_colon_iter.Get()==1) {		
	//			CovariantVectorType partial = partialVector_iter.Get();		//retrieves the partial informations
	//			float Z[3]={partial[0],partial[1],partial[2]};

	//			// Only update uncertain partials
	//			if ( !(partial[0] == 1 || partial[1] == 1 || partial[2] == 1) )
	//			{
	//				vnl_vector<float> Z_update=expectation(input_iter.Get(),mean, variance, weight, GetNeighbor(partialVector,input_iter.GetIndex()), Z);	//updates the partial values
	//				partial[0]=Z_update[0];										
	//				partial[1]=Z_update[1];
	//				partial[2]=Z_update[2];
	//				partialVector_iter.Set(partial);
	//			}

	//			//updates the new mean total partial sum for each class accordingly
	//			for (int i=0;i<3;i++) 
	//			{
	//				mean_temp[i]+=partial[i]*input_iter.Get();	
	//				sum_temp[i]+=partial[i];
	//			}

	//		}
 //       }

	//	partialVector_iter = IteratorImageVectorType(partialVector,fullRegion);

	//	for (int i=0;i<3;i++) { mean_temp[i]=mean_temp[i]/sum_temp[i]; } 

	//	// Compute variance and weights
	//	for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
	//		!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
	//		++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//	{
	//		if (chamfer_colon_iter.Get()==1) {		
	//			CovariantVectorType partial = partialVector_iter.Get();
	//			for (int i=0; i<3; i++)
	//			{
	//				variance_temp[i] += partial[i]*vnl_math_sqr(input_iter.Get()-mean_temp[i]);
	//			}
	//		}
	//	}

	//	double sum_all_temp = sum_temp[0]+sum_temp[1]+sum_temp[2];

	//	for (int i=0;i<3;i++) 
	//	{ 
	//		mean[i]=mean_temp[i];
	//		variance[i]=variance_temp[i]/sum_temp[i];
	//		weight[i]=sum_temp[i]/sum_all_temp;
	//	}

	//	std::cerr<<std::endl;
	//	std::cerr<<"EM"<<emNum+1<<std::endl;
	//	std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	//	std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	//	std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;

	//	std::stringstream ss;
	//	ss<<emNum+1;
	//	ss<<"\t"<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<variance[0]<<"\t"<<variance[1]<<"\t"<<variance[2]<<"\t"<<weight[0]<<"\t"<<weight[1]<<"\t"<<weight[2]<<"\n";
	//	em<<ss.str();

	//	std::stringstream ss2;
	//	ss2<<"EM"<<emNum+1;

	//	WritePartialImages(partialVector,chamfer_colon,ss2.str());
 //   }

	//em.close();

	////set modification images type.
	////Updates the existing voxel type classifications with knowledge found in the EM
	//EMClassification(input_iter, voxel_type_iter, partialVector_iter, chamfer_colon_iter, temp_iter);
	//WriteITK(temp,"Global_EM_classification.nii");

	//for(input_iter.GoToBegin(), em_output_image_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
	//	!input_iter.IsAtEnd() && !em_output_image_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
	//	++input_iter, ++em_output_image_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//{
	//	CovariantVectorType data = partialVector_iter.Get();		//retrieves the partial informations
	//	em_output_image_iter.Set(data[1]);		//sets the voxel intensity for the i-th iteration output image
	//}

	//WriteITK(em_output_image, "tissue_partial.nii");

	//std::cerr<<"Finished EM"<<std::endl;
	//em_output_image.~SmartPointer();
	//partialVector2.~SmartPointer();
	//std::cerr<<"partialVector2 has been destroyed: "<<partialVector2.IsNull()<<std::endl;

	//-------------------------------------------END EM-------------------------------------------------------

//-------------------------------------------BEGIN PARTIAL VOLUME COMP----------------------------------
	
	//std::cerr<<"Compute partial volumes "<<std::endl;
 //   //0: ps
 //   //1: pt
 //   //2: pa
 //   ImageType::Pointer partialVolume = AllocateNewImage(fullRegion);
 //   IteratorTypeFloat4WithIndex partialVolume_iter(partialVolume,fullRegion);
	////bool random=true;

	//ByteImageType::Pointer chamfer_from_stool_involvement=AllocateNewByteImage(fullRegion);
	//IteratorTypeByteWithIndex chamfer_from_stool_involvement_iter(chamfer_from_stool_involvement,fullRegion);
	//
	//std::cerr << "Partial volume computation..." << std::endl;
	////Sets the partial volumes
 //   for(voxel_type_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin(), input_smax_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
 //       !voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
 //       ++voxel_type_iter, ++input_iter, ++partialVolume_iter, ++input_smax_iter, ++chamfer_from_stool_involvement_iter) 
 //   {
	//	//prepares information for chamfer stool calculation
	//	if (voxel_type_iter.Get()==Stool || voxel_type_iter.Get()==TissueStool || voxel_type_iter.Get()==StoolAir) {
	//		chamfer_from_stool_involvement_iter.Set(1);     //stool
	//	} else {
	//		chamfer_from_stool_involvement_iter.Set(0);     //non stool is the object
	//	}
	//	partialVolume_iter.Set(0);
	//	float value=0;
 //       switch(voxel_type_iter.Get()) {
	//		case Tissue:
	//			partialVolume_iter.Set(1.0);
	//			break;
	//		case Stool:
	//		case Air:
	//			partialVolume_iter.Set(0.0);
	//			break;
 //           case TissueAir:
	//			//std::cerr<<"TA reached"<<std::endl;
	//			
	//			//partialVolume_iter.Set(1.0);	
	//			value = 1.0+(input_iter.Get()-1024)/1000;
	//			if (value<=1 && value>=0) {
	//				partialVolume_iter.Set(value);
	//				std::cout << "Set partial volume correctly!" << std::endl;
	//			} else if (value<0) {
	//				partialVolume_iter.Set(0);	
	//			} else {
	//				partialVolume_iter.Set(1);	
	//			}
 //               break;
 //           case TissueStool:
	//				value=1.0-(input_iter.Get())/(input_smax_iter.Get());
	//				if (value<=1 && value>=0) {
	//					partialVolume_iter.Set(value);
	//				} else if (value<0) {
	//					partialVolume_iter.Set(0);	
	//				} else {
	//					partialVolume_iter.Set(1);	
	//				}

	//				////special case when there is nothing really to fit
	//				//if (input_smax_iter.Get()-input_iter.Get()<50 && input_iter.Get()-1024<500) {
	//				//	partialVolume_iter.Set(1);	
	//				//}
 //               break;
 //           case StoolAir:
 //               partialVolume_iter.Set(0);
 //               break;
 //       }
	//	if (-partialVolume_iter.Get()==std::numeric_limits<float>::infinity()) {
	//		//std::cerr<<partialVolume_iter.GetIndex()[0]<<" "<<partialVolume_iter.GetIndex()[1]<<" "<<partialVolume_iter.GetIndex()[2]<<std::endl;
	//	}
 //   }

	//input_smax.~SmartPointer();

	////Write out partial volume
	//WriteITK(partialVolume,"tissue_partial_one.nii");

	////Write out unmodified image
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_input.nii");

	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_from_stool_involvement);
	//int weights[3]={3,4,5};	//3d distance weight recommended by julian
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(false);
	//chamfer_filter->Update();
	//chamfer_from_stool_involvement=chamfer_filter->GetOutput();
	//chamfer_from_stool_involvement_iter = IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	//chamfer_filter.~SmartPointer();


	////Write out chamfer thickness of stool
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool.nii");

	////Updates chamfer by setting all stool attached to an 8+ chamfer stool to also 8 chamfer units
	//int chamfercutoff=8;
	//for (int i=0;i<2;i++) {
	//	//Forward
	//	for(chamfer_from_stool_involvement_iter.GoToBegin(); !chamfer_from_stool_involvement_iter.IsAtEnd(); ++chamfer_from_stool_involvement_iter){
	//		ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
	//		float temp_value=chamfer_from_stool_involvement_iter.Get();
	//		if (temp_value>=chamfercutoff) {
	//			for(int i=-1;i<=1;i++) {
	//				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
	//					for (int j=-1;j<=1;j++) {
	//						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
	//							for (int k=-1;k<=1;k++) {
	//								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
	//									ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
	//									if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
	//										chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
	//										//std::cout << "Exploding Eights Forward!" << std::endl;
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}

	//	//Write the modified chamfer map
	//	WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_forward.nii");
	//	
	//	//Backward
	//	for(chamfer_from_stool_involvement_iter.GoToReverseBegin(); !chamfer_from_stool_involvement_iter.IsAtReverseEnd(); --chamfer_from_stool_involvement_iter){
	//		ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
	//		float temp_value=chamfer_from_stool_involvement_iter.Get();
	//		if (temp_value>=chamfercutoff) {
	//			for(int i=-1;i<=1;i++) {
	//				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
	//					for (int j=-1;j<=1;j++) {
	//						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
	//							for (int k=-1;k<=1;k++) {
	//								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
	//									ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
	//									if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
	//										chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
	//										//std::cout << "Exploding Eights Backward!" << std::endl;
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//	//write the modified chamfer map
	//	WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_back.nii");
	//}

	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool.nii");



 //   ByteImageType::Pointer chamfer_air=AllocateNewByteImage( fullRegion );
 //   IteratorTypeByteWithIndex chamfer_air_iter(chamfer_air,fullRegion);

 //   for(voxel_type_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin();
 //       !voxel_type_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd();
 //       ++voxel_type_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter, ++input_iter, ++partialVolume_iter)
 //   {
	//	//prepare for computation of chamfer to air
	//	if (voxel_type_iter.Get()==Air) {
 //           chamfer_air_iter.Set(0);					//distance to AIR
 //       } else {
 //           chamfer_air_iter.Set(1);					//non important locations
 //       }
	//	float temp_value=chamfer_from_stool_involvement_iter.Get();
	//	//the temp_value>3 is uncessary
	//	if ((temp_value>=chamfercutoff && temp_value>3 && voxel_type_iter.Get()!=TissueStool) || voxel_type_iter.Get()==Stool) {
	//		partialVolume_iter.Set(0);
	//		//std::cout << "Set partial tissue volume to 0!" << std::endl;
	//	} else if (temp_value<chamfercutoff && temp_value>3) {
	//		voxel_type_iter.Set(ThinStool);		//if the chamfer < 8 then it is set to thin stool
	//		//std::cout << "Set ThinStool!" << std::endl;
	//	}
 //   }

	////compute chamfer for air
	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_air);
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(false);
	//chamfer_filter->Update();
	//chamfer_air=chamfer_filter->GetOutput();
	//chamfer_air_iter= IteratorTypeByteWithIndex(chamfer_air,fullRegion);
	//chamfer_filter.~SmartPointer();
	////Writes out the chamfer to air
	//WriteITK(chamfer_air,"chamfer_from_air.nii");
	////normalize the data

 //   for(chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
 //       !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
 //       ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter)
 //   {
 //       chamfer_from_stool_involvement_iter.Set(chamfer_from_stool_involvement_iter.Get()/chamfercutoff);
 //       if (chamfer_air_iter.Get()>=chamfercutoff) {
 //           chamfer_air_iter.Set(1);
 //       } else {
 //           chamfer_air_iter.Set(chamfer_air_iter.Get()/8);
 //       }
 //   }

 //   //Calculating partial air
 //   for(voxel_type_iter.GoToBegin(), partialVolume_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
 //       !voxel_type_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
	//	++voxel_type_iter, ++partialVolume_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter) 
	//{
 //       if(voxel_type_iter.Get()==ThinStool) {
	//		float value=0;
	//		float ps=1/2*(1+vnl_erf((chamfer_from_stool_involvement_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
	//		float pa=1/2*(1+vnl_erf((chamfer_air_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
	//		value=1-ps-pa;
	//		partialVolume_iter.Set(value);
 //       }
	//}

	//for(partialVolume_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(); 
	//	!partialVolume_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd(); 
	//	++partialVolume_iter, ++chamfer_colon_iter, ++voxel_type_iter) 
	//{
	//	if (chamfer_colon_iter.Get()==1) { 
	//		if (partialVolume_iter.Get()<.3) {
	//			partialVolume_iter.Set(0);
	//		}
	//		if (voxel_type_iter.Get() == TissueAir) {
	//			partialVolume_iter.Set(1);
	//		}
	//	}
	//}

	////Write out the partial tissue with thinstool modification done
	//WriteITK(partialVolume,"tissue_partial_two.nii");

	////Prepares the output image
	//ImageType::Pointer output = AllocateNewImage(fullRegion);
	//output->SetSpacing(input->GetSpacing());
	//IteratorTypeFloat4WithIndex output_iter(output,fullRegion);
 //   //Modifies to intensity for output
 //   for(output_iter.GoToBegin(), input_iter.GoToBegin(), temp_iter.GoToBegin(), voxel_type_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), partialVolume_iter.GoToBegin();
	//	!output_iter.IsAtEnd() && !input_iter.IsAtEnd() && !temp_iter.IsAtEnd() &&  !voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd();
	//	++output_iter, ++input_iter, ++temp_iter, ++voxel_type_iter, ++chamfer_colon_iter, ++partialVolume_iter) 
 //   {	
	//	if (voxel_type_iter.Get()==Air || voxel_type_iter.Get()==StoolAir) {
	//		//output_iter.Set(mean[0]);
	//		output_iter.Set(-1024);
	//	} else if (chamfer_colon_iter.Get()==1 || voxel_type_iter.Get()==TissueStool  || voxel_type_iter.Get()==Stool) { 
	// //       ImageType::IndexType index =input_iter.GetIndex();
	//		output_iter.Set(partialVolume_iter.Get()*input_iter.Get());
	//	} else {
	//		output_iter.Set(input_iter.Get());
	//	}
	//}
	////voxel_type.~SmartPointer();
	//chamfer_air.~SmartPointer();
	////Writes out the output preclosing
	//WriteITK(output,"output_one.nii");
	//for(output_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), voxel_type_iter.GoToBegin();
 //       !output_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd();
	//	++output_iter, ++chamfer_from_stool_involvement_iter, ++voxel_type_iter) 
	//{
	//	if (output_iter.Get()>300 || voxel_type_iter.Get()==Tissue || voxel_type_iter.Get()==TissueAir) {
	//		chamfer_from_stool_involvement_iter.Set(1);
	//	} else {
	//		chamfer_from_stool_involvement_iter.Set(0);
	//	}
	//}

	////Now chamfer_from_stool_involvement is a map to known tissue

	////Write out unmodified image
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_tissue_input.nii");
	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_from_stool_involvement);
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(false);
	//chamfer_filter->Update();

	//chamfer_from_stool_involvement = chamfer_filter->GetOutput();
	//chamfer_from_stool_involvement_iter=IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	//chamfer_filter.~SmartPointer();

	////Write out chamfer thickness of stool
	//WriteITK(chamfer_from_stool_involvement,"chamfer_tissue_stool.nii");

	////-------------------------------------------BEGIN PARTIAL VOLUME COMP----------------------------------
	

	//Write out the final output	
	std::cerr<<"Ended"<<std::endl;
	//WriteITK(output, "output.nii");
	//system("pause");
	return 0;
}

double Probability(double Y, double mean, double variance,  double current_partial, double local_variance, double local_mean) {  
	if (variance>0.01) {
		if (local_variance>0.01) { //?
		//return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance))*exp(-beta*total_neighbor);
			//return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance))*exp(-vnl_math_sqr(current_partial-local_mean)/(2*local_variance))/(sqrt(2*PI*local_variance));
			return exp(-vnl_math_sqr(current_partial-local_mean)/(2*local_variance))/(sqrt(2*PI*local_variance));
		} else {
			//return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance));
			return exp(-vnl_math_sqr(current_partial-local_mean)/(2*local_variance))/(sqrt(2*PI*local_variance));
		}
	} else {
		return 0;
	}
}

vnl_vector<float> expectation(double Y, double mean[], double variance[], float weight[],  vnl_matrix<float> neighbor, float current_partial[]) {
	int nKernal=3;
	float pFK[3]={0,0,0};	//probabilities of each class: air/tissue/stool per voxel
	float pKF[4]={0,0,0};   //""	"", [3] = sum of probabilities of air/tissue/stool
	float sumPFK=0.0;
	for (int i=0;i<3;i++) {
		double total_neighbor=0;
		
		/*
		if (current_partial[0]<.6 && current_partial[1]<.6 && current_partial[2]<.6) {
			beta=.7;
		} else {
			beta=0;
		}*/
		
		// Local neighborhood probability inputs
		float local_mean=0;
		float sum=0;
		for (int j=0;j<6;j++) {
			local_mean+=neighbor.get(i,j)*neighbor_weight[j % 3];
			sum+=neighbor_weight[j % 3];
		}
		local_mean=local_mean/sum;
		
		float local_variance=0;
		for (int j=0;j<6;j++) {
			local_variance+=vnl_math_sqr(local_mean-neighbor.get(i,j))*neighbor_weight[j % 3];
		}

		local_variance=local_variance/sum;
			//if (gradient>100) {
			
				//switch (i) {
				//	case 0:
				//		total_neighbor+=sqrt(vnl_math_sqr(current_partial[1]-neighbor.get(1,j))*vnl_math_sqr(current_partial[2]-neighbor.get(2,j)))*neighbor_weight[j % 3];
				//		break;
				//	case 1:
				//		total_neighbor+=sqrt(vnl_math_sqr(current_partial[0]-neighbor.get(0,j))*vnl_math_sqr(current_partial[2]-neighbor.get(2,j)))*neighbor_weight[j % 3];
				//		break;
				//	case 2:
				//		total_neighbor+=sqrt(vnl_math_sqr(current_partial[1]-neighbor.get(1,j))*vnl_math_sqr(current_partial[0]-neighbor.get(0,j)))*neighbor_weight[j % 3];
				//		break;
				//}
			//}
		//std::cerr<<sum<<" "<<local_variance<<" "<<local_mean<<std::endl;
		
		// Calculate probability
		pFK[i]=Probability(Y,mean[i],variance[i], current_partial[i], local_variance, local_mean)*weight[i];
		
		// Sum probabilities across classes
		sumPFK+=pFK[i];
	}
	
	// Normalize probabilities across classes so that prob(air/tissue/stool)=1
	for(int i=0;i<3;i++) {
		pKF[i]=pFK[i]/sumPFK;
	}
	pKF[3]=sumPFK;
	vnl_vector<float> return_value(4, 4, pKF);	//Return prob(air/tissue/stool)
	return return_value;
}



vnl_vector<float> expectation(double Y, float mean[], float variance[], float weight[], float current_partial[]) {
	float pFK[2]={0,0};
	float pKF[3]={0,0,0};
	float sumPFK=0.0;

	for (int i=0;i<2;i++) {
		pFK[i]=Probability(Y,mean[i],variance[i], 0, 0, 0)*weight[i];
		sumPFK+=pFK[i];
	}

	for(int i=0;i<2;i++) {
		pKF[i]=pFK[i]/sumPFK;
	}
	pKF[2]=sumPFK;
	vnl_vector<float> return_value(3, 3, pKF);
	return return_value;
}

vnl_matrix<float> GetNeighbor(ImageVectorType::Pointer partialVector, ImageType::IndexType index) {
    ImageType::RegionType fullregion = partialVector->GetLargestPossibleRegion();
    vnl_matrix<float> temp_return(3,6);
    float filler[3]={0,0,0};
    for(int i=0;i<3;i++) {
        ImageType::IndexType temp_index_1(index);
        temp_index_1[i]+=1;
        CovariantVectorType data_1(filler);
        if (fullregion.IsInside(temp_index_1)) {
            data_1 = partialVector->GetPixel(temp_index_1);    //account for image issues
        }

        ImageType::IndexType temp_index_2(index);
        temp_index_2[i]-=1;
        CovariantVectorType data_2(filler);
        if (fullregion.IsInside(temp_index_2)) {
            data_2 = partialVector->GetPixel(temp_index_2);
        } 

        for (int j=0;j<3;j++) {
            temp_return.put(j,i,data_1[j]);
			temp_return.put(j,i+3,data_2[j]);
        }
    }
    return temp_return;
}




void EMClassification(IteratorTypeFloat4WithIndex input_iter, IteratorTypeVoxelType voxel_type_iter, IteratorImageVectorType partialVector_iter,
					  IteratorTypeByteWithIndex chamfer_colon_iter, IteratorTypeFloat4WithIndex temp_iter) {
	
	// count number of changes made by EM
	int count = 0;

	for(input_iter.GoToBegin(), voxel_type_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), temp_iter.GoToBegin(); 
		!input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd(); 
		++input_iter, ++voxel_type_iter, ++partialVector_iter, ++chamfer_colon_iter, ++temp_iter) 
	{
		if (chamfer_colon_iter.Get()==1) { 
			CovariantVectorType data = partialVector_iter.Get();	//gets the partial information
			switch (voxel_type_iter.Get()) {
				case Unclassified:
					//if an unclassified voxel has a high partial (probability) of being a solid 
					//class and doesn't have a high gradient then reclassify it as the solid class
					//
					//it seems to mis classify a lot of stool and air, so for now just reclassifying the solid tissues
					
					if (data[0]>.9 && temp_iter.Get()<250) {
						voxel_type_iter.Set(Air);        //sets the voxel type
						std::cout << "Changed Air!" << std::endl;
						count++;
					} else if (data[1]>.8 && temp_iter.Get()<300*data[1] && (input_iter.Get()<600)) {
						voxel_type_iter.Set(Tissue);        //sets the voxel type
						std::cout << "Changed Tissue!" << std::endl;
						count++;
					} else if (data[2]>.9 && temp_iter.Get()<0.8*input_iter.Get()) {
						voxel_type_iter.Set(Stool);        //sets the voxel type
						std::cout << "Changed Stool!" << std::endl;
						count++;
					} else {
						voxel_type_iter.Set(Unclassified);        //sets the voxel type
					}
					break;
					
				case Tissue:
					if (data[1]<.1 /*&& temp_iter.Get()>300*/) {
						voxel_type_iter.Set(Unclassified);
						std::cout << "Changed Tissue to Unclassified!" << std::endl;
						count++;
					}
					break;
			}
		}

		switch(voxel_type_iter.Get()) {
			case Stool:
                temp_iter.Set(2);
                break;
            case Tissue:
                temp_iter.Set(1);
				break;
            case Unclassified:
                temp_iter.Set(-1);
                break;
            case Air:
                temp_iter.Set(0);
                break;
        }
	}

	std::cerr << "Number of voxels changed by EM: " << count << std::endl;

}

void VoxelEdgeClassification(float * threshold, VoxelType * previous, double d2, double d1,

                                      InterpolationType::Pointer &input_interpolator, 

									  InterpolationType::Pointer &gradient_magnitude_interpolator,

                                      IteratorTypeFloat4WithIndex &input_smax,

                                      ImageType::IndexType &index,
									
									  CovariantVectorType &gradient) 

{
	ImageType::RegionType fullRegion = input_interpolator->GetInputImage()->GetLargestPossibleRegion();
    InterpolationType::ContinuousIndexType offset_index;
    offset_index[0]=index[0];
    offset_index[1]=index[1];
    offset_index[2]=index[2];
	ImageType::IndexType endIndex = input_interpolator->GetEndIndex();
	ImageType::IndexType startIndex = input_interpolator->GetStartIndex();

    //stores the intensity
    float temp_intensity[5];
    float temp_gradient_magnitude[5];
	//stores the target voxel
    temp_intensity[2]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[2]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
    //+d1
    offset_index[0]=index[0]+gradient[0]*d1;
    offset_index[1]=index[1]+gradient[1]*d1;
    offset_index[2]=index[2]+gradient[2]*d1;
    temp_intensity[3]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[3]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
    //-d1
    offset_index[0]=index[0]-gradient[0]*d1;
    offset_index[1]=index[1]-gradient[1]*d1;
    offset_index[2]=index[2]-gradient[2]*d1;
    temp_intensity[1]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[1]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
    //+d2
    offset_index[0]=index[0]+gradient[0]*d2;
    offset_index[1]=index[1]+gradient[1]*d2;
    offset_index[2]=index[2]+gradient[2]*d2;
    temp_intensity[4]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[4]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
    //-d2
    offset_index[0]=index[0]-gradient[0]*d2;
    offset_index[1]=index[1]-gradient[1]*d2;
    offset_index[2]=index[2]-gradient[2]*d2;
    temp_intensity[0]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[0]=gradient_magnitude_interpolator->EvaluateAtContinuousIndex(offset_index);
	//std::cerr<<"vector created"<<std::endl;
	float min_distance=*threshold;
	float stool_tissue_Smax=0;
	float stool_air_Smax=0;
	float smax=0;
	float distance=0;
	VoxelType voxel_type=Unclassified;

	float distanceTS=0;
	float distanceSA=0;
	float distanceTA=0;



	if (!Modified)
	{
		stool_tissue_Smax=ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
		stool_air_Smax=Stool_Air_ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
		smax = 0;

		distanceTS=AverageTissueStoolDist(stool_tissue_Smax, temp_intensity,temp_gradient_magnitude);
		distanceSA=AverageStoolAirDist(stool_air_Smax, temp_intensity,temp_gradient_magnitude);
	} else {
		smax = input_smax.Get();

		distanceTS=AverageTissueStoolDist(smax, temp_intensity,temp_gradient_magnitude);
		distanceSA=AverageStoolAirDist(smax, temp_intensity,temp_gradient_magnitude); 
	}	
	
	distanceTA=AverageTissueAirDist(temp_intensity,temp_gradient_magnitude);

	if (distanceSA<=distanceTS && distanceSA<=distanceTA) {
		distance=distanceSA;
		voxel_type=StoolAir;

		if (!Modified)
			smax=stool_air_Smax;

	} else if (distanceTS<=distanceTA && distanceTS<=distanceSA) {
		distance=distanceTS;
		voxel_type=TissueStool;

		if (!Modified)
			smax=stool_tissue_Smax;

	} else {
		distance=distanceTA;
		voxel_type=TissueAir;

		if (!Modified)
			smax=0;
	}

    if (min_distance>=distance || min_distance==-1 || *previous==Unclassified) {
		*threshold=distance;
        *previous=voxel_type;
		if (!Modified)
			input_smax.Set(smax);
    }
}

float Stool_Air_ComputeSmax(float intensity[], float gradient_magnitude[], int size) {
	//std::cerr<<"Value ";
	float intensity2[5];
	for (int i=0;i<size;i++) {
		intensity2[i]=intensity[i]+1000;
	//	std::cerr<<intensity2[i]<<" ";
	}
	//std::cerr<<std::endl;
	return ComputeSmax(intensity2,gradient_magnitude, size)-1000;
}
float ComputeSmax(float intensity[], float gradient_magnitude[], int size) {
	double Smax =0;
	double subVariable=0;
	double square =0;
	for (int i=0;i<size;i++) {
		square = intensity[i]*intensity[i];
		Smax+=square*gradient_magnitude[i];
		subVariable+=square*intensity[i];
	}
	Smax-=4*subVariable;
	subVariable=0;
	for (int i=0;i<size;i++) {
		square = intensity[i]*intensity[i];
		square *= square;
		subVariable+=square;
	}
	return -4*subVariable/Smax;
}
float ComputeSmaxFit(float intensity[], float gradient_magnitude[], float Smax) {
	float R=0;
	for (int i=0;i<5;i++) {
		float difference=gradient_magnitude[i]-(-4/Smax*intensity[i]*intensity[i]+4*intensity[i]);
		R+=	difference*difference;
	}
	return R;
}
float AverageTissueAirDist(float intensity[], float gradient_magnitude[]) {
//	std::cerr<<"TissueAir"<<std::endl;
	double coefficients[3]={-4.0/1000,-4.0,0};

	vnl_real_polynomial poly(coefficients,3);
    double average=0;
    for(int i=0;i<5;i++) {
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }
    return average/5;
}
float AverageTissueStoolDist(float Smax, float intensity[], float gradient_magnitude[]) {
//	std::cerr<<"TissueStool"<<std::endl;
    double coefficients[3] ={-4.0/Smax,4.0,0};
    vnl_real_polynomial poly(coefficients,3);
    double average=0;
    for(int i=0;i<5;i++) {
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }
    return average/5;
}

float AverageStoolAirDist(float Smax, float intensity[], float gradient_magnitude[]) {
//	std::cerr<<"StoolAir"<<std::endl;
    double coefficients[3] ={-4/(Smax+1000),-4*(1000-Smax)/(Smax+1000),4000*Smax/(Smax+1000)};
    vnl_real_polynomial poly(coefficients,3);
    double average=0;
    for(int i=0;i<5;i++) {
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }
    return average/5;
}

float PolyMinDist(vnl_real_polynomial poly, float x, float y) {
    double coefficients[4] ={
		vnl_math_sqr(poly[0])*2,
		poly[0]*3*poly[1],
		(2*poly[0]*poly[2]-2*poly[0]*y+vnl_math_sqr(poly[1])+1),
        poly[1]*poly[2]-poly[1]*y-x
    };
	
    vnl_real_polynomial polyRoot(coefficients,4);
    vnl_rpoly_roots rooter(polyRoot);
    //rooter.compute();
    vnl_vector<double> rRoots = rooter.realroots(0.0005);
    int size = rRoots.size();
    float distance = PolyDist(rRoots(0),poly,x,y);
    for(int i=1;i<size;i++) {
        float temp_dist=PolyDist(rRoots(i),poly,x,y);
        if (distance>temp_dist) {
            distance=temp_dist;
        }
    }
	//std::cerr<<distance<<std::endl;
    return distance;
}

float PolyDist(float X, vnl_real_polynomial poly, float x, float y) {
    return sqrtf(vnl_math_sqr(X-x)+vnl_math_sqr(poly.evaluate(X)-y));
}

void SingleMaterialClassification(ImageType::Pointer input, ImageType::Pointer gradient, VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon) 
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex gradient_iter(gradient,region);
	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);

	// Determine tissue stool intensity threshold using otsu threshold on histogram

	// Mask input with colon segmentation
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( input );
	masker->SetInput2( chamfer_colon );
	masker->SetOutsideValue( -1500 ); // arbitrary
	masker->Update();
    input = masker->GetOutput();
	input_iter = IteratorTypeFloat4WithIndex(input,region);

	// Compute otsu threshold for tissue and stool intensity
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-250);
	otsuCalculator->SetHistogramMax(1400);
	otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	otsuCalculator->Compute();

	float tissue_stool_threshold = otsuCalculator->GetThreshold();

	std::cout << std::endl << std::endl;
	std::cout << "Tissue Stool Otsu Threshold: " << tissue_stool_threshold << std::endl;
	std::cout << std::endl << std::endl;

	// Apply thresholds and save voxel type

	for ( input_iter.GoToBegin(), gradient_iter.GoToBegin(), voxel_type_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++gradient_iter, ++voxel_type_iter, ++chamfer_colon_iter)
	{
		if (chamfer_colon_iter.Get() == 1)
		{
			VoxelType voxel;

			float I = input_iter.Get();
			float G = gradient_iter.Get();

			if (Modified) 
			{

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

			} else {
				
				if ( I >= 180 && G < 0.8*I )
				{
					voxel = Stool;
				} else if ( I <= -800 && G <= 250) {
					voxel = Air;
				} else if ( I < 150 && I > -250 && G <= 300 ) {
					voxel = Tissue;
				} else {
					voxel = Unclassified;
				}
			}

			voxel_type_iter.Set( voxel );
		}
	}
}


ByteImageType::Pointer AllocateNewByteImage(ImageType::RegionType fullRegion) {
    ByteImageType::Pointer newImage = ByteImageType::New();
    newImage->SetLargestPossibleRegion( fullRegion );
    newImage->SetBufferedRegion( fullRegion );
    newImage->SetRequestedRegion( fullRegion );
    newImage->Allocate();
    return newImage;
}

ImageType::Pointer AllocateNewImage(ImageType::RegionType fullRegion) {
    ImageType::Pointer newImage = ImageType::New();
    newImage->SetLargestPossibleRegion( fullRegion );
    newImage->SetBufferedRegion( fullRegion );
    newImage->SetRequestedRegion( fullRegion );
	newImage->Allocate();
    return newImage;
}

ImageType::Pointer ComputePseudoGradient(ImageType::Pointer input) {
	ImageType::RegionType fullRegion=input->GetLargestPossibleRegion();
	ImageType::Pointer temp = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex input_iter(input,fullRegion);
	IteratorTypeFloat4WithIndex temp_iter(temp,fullRegion);
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);
	for(input_iter.GoToBegin(), temp_iter.GoToBegin(); !input_iter.IsAtEnd() && !temp_iter.IsAtEnd(); ++input_iter, ++temp_iter){
		float gradient=100000;
		ImageType::IndexType index = input_iter.GetIndex();
		if (index[0]+1<=endIndex[0] && index[0]-1>=startIndex[0] && 
			index[1]+1<=endIndex[1] && index[1]-1>=startIndex[1] &&
			index[2]+1<=endIndex[2] && index[2]-1>=startIndex[2]) {
			for(int i=-1;i<=1;i++) {
				for (int j=-1;j<=1;j++) {
					for (int k=-1;k<=1;k++) {
						if (i!=0 || j!=0 || k!=0) {
							ImageType::IndexType temp_index_1={index[0]+i,index[1]+j,index[2]+k};
							ImageType::IndexType temp_index_2={index[0]-i,index[1]-j,index[2]-k};
							float temp_gradient=(input->GetPixel(temp_index_1)-input->GetPixel(temp_index_2));
							temp_gradient=sqrtf(temp_gradient*temp_gradient/(vnl_math_sqr(i*input->GetSpacing()[0])+vnl_math_sqr(j*input->GetSpacing()[1])+vnl_math_sqr(k*input->GetSpacing()[2])));
							//std::cerr<<temp_gradient<<std::endl;
							if (temp_gradient<gradient) {
								gradient=temp_gradient;
							}
						}
					}
				}
			}
		} else {
			gradient=0;
		}
		if (gradient>75) {
		//	std::cerr<<gradient<<std::endl;
		}
		temp_iter.Set(gradient);
	}
	return temp;
}

double round(float d)
{
  return floor(d + 0.5);
}

void WriteITK(ImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";

	if ( writeNum )
		ss << writeCount++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void WriteITK(FuzzySceneType::Pointer image, std::string name) {
	// Use ITK Writer
	typedef itk::ImageFileWriter< FuzzySceneType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";

	if ( writeNum )
		ss << writeCount++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);

	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";

	if ( writeNum )
		ss << writeCount++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();
}

void WriteITK(IntImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< IntImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);

	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";

	if ( writeNum )
		ss << writeCount++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();
}

void WriteITK(UintImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< UintImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);

	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";

	if ( writeNum )
		ss << writeCount++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();
}

void WriteITK(VoxelTypeImage::Pointer vimage, std::string name) {
	// Get voxel type iter
	IteratorTypeVoxelType vit( vimage, vimage->GetLargestPossibleRegion() );
	
	// Create helper image
	ImageType::Pointer temp = ImageType::New();
	temp->SetRegions( vimage->GetLargestPossibleRegion() );
	temp->Allocate();
	IteratorTypeFloat4WithIndex tit( temp, temp->GetLargestPossibleRegion() );

	for (vit.GoToBegin(), tit.GoToBegin(); !vit.IsAtEnd() && !tit.IsAtEnd(); ++vit, ++tit)
	{
		 switch(vit.Get()) {
            case Stool:
                tit.Set(1);
                break;
			case Air:
				tit.Set(2);
				break;
			case Tissue:
				tit.Set(3);
				break;
			case Unclassified:
				tit.Set(4);
				break;
			case StoolAir:
				tit.Set(5);
				break;
			case TissueAir:
				tit.Set(6);
				break;
			case TissueStool:
				tit.Set(7);
				break;
			case ThinStool:
				tit.Set(8);
				break;
			default:
                tit.Set(0);
                break;
        }
	}

	WriteITK(temp, name);

}

void ReadITK(ImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	image = reader->GetOutput();
}

void ReadITK(ByteImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ByteImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	image = reader->GetOutput();
}

void ReadITK(VoxelTypeImage::Pointer &vimage, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ByteImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fileName);
	reader->Update();
	ByteImageType::Pointer vim = reader->GetOutput();
	IteratorTypeByteWithIndex vim_iter(vim, vimage->GetLargestPossibleRegion() );

	IteratorTypeVoxelType voxel_type_iter(vimage,vimage->GetLargestPossibleRegion());

	for (vim_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !vim_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd(); ++vim_iter, ++voxel_type_iter)
	{
		 switch(vim_iter.Get()) {
            case 1:
                voxel_type_iter.Set(Stool);
                break;
			case 2:
				voxel_type_iter.Set(Air);
				break;
			case 3:
				voxel_type_iter.Set(Tissue);
				break;
			case 4:
				voxel_type_iter.Set(Unclassified);
				break;
			case 5:
				voxel_type_iter.Set(StoolAir);
				break;
			case 6:
				voxel_type_iter.Set(TissueAir);
				break;
			case 7:
				voxel_type_iter.Set(TissueStool);
				break;
			case 8:
				voxel_type_iter.Set(ThinStool);
				break;
        }
	}
}

bool compareSizeOnBorder(LabelObjectType::Pointer a, LabelObjectType::Pointer b)
{
//Comparison function object that, taking two values of the same type than those contained in the range, returns true if the first 
//argument goes before the second argument in the specific strict weak ordering it defines, and false otherwise.

	return (a->GetSizeOnBorder() > b->GetSizeOnBorder());
}

bool compareSize(LabelObjectType::Pointer a, LabelObjectType::Pointer b)
{
//Comparison function object that, taking two values of the same type than those contained in the range, returns true if the first 
//argument goes before the second argument in the specific strict weak ordering it defines, and false otherwise.

	return (a->GetSize() > b->GetSize());
}

int VoxelTypeToNum(VoxelType type)
{
	switch(type) {
		case Stool:
			return 1;
		case Air:
			return 2;
		case Tissue:
			return 3;
		case Unclassified:
			return 4;
		case StoolAir:
			return 5;
		case TissueAir:
			return 6;
		case TissueStool:
			return 7;
		case ThinStool:
			return 8;
	}
	return 0;
}

std::string VoxelTypeToString(VoxelType type)
{
	switch(type) {
		case Stool:
			return "Stool";
		case Air:
			return "Air";
		case Tissue:
			return "Tissue";
		case Unclassified:
			return "Unclassified";
		case StoolAir:
			return "StoolAir";
		case TissueAir:
			return "TissueAir";
		case TissueStool:
			return "TissueStool";
		case ThinStool:
			return "ThinStool";
	}
	return "Error";
}

VoxelType NumToVoxelType(int num)
{
	switch(num)
	{
		case 1:
			return Stool;
		case 2:
			return Air;
		case 3:
			return Tissue;
		case 4:
			return Unclassified;
		case 5:
			return StoolAir;
		case 6:
			return TissueAir;
		case 7:
			return TissueStool;
		case 8:
			return ThinStool;
		default:
			return Unclassified;
	}
	return Unclassified;
}

template <typename T>
typename T::Pointer ReadDicom( std::string path, int slice1=0, int slice2=-1)
{	
	// Create reader
	itk::ImageSeriesReader<T>::Pointer reader = itk::ImageSeriesReader<T>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	/*
	
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	
	nameGenerator->SetDirectory( path );
	
	const std::vector< std::string > seriesUID = nameGenerator->GetSeriesUIDs();
	std::string seriesIdentifier = seriesUID.begin()->c_str();
	
	std::vector< std::string > names = nameGenerator->GetFileNames( seriesIdentifier );

	*/
	
	
	
	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	if (slice2 > 0 && slice2 > slice1)
	{
		names.erase( names.begin(), names.begin()+slice1);
		names.erase( names.begin()+slice2-slice1, names.end() );
	}
	

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error reading dicom: " << err << std::endl;
		return 0;
	}
	
	T::Pointer output = reader->GetOutput();
	
	/*
    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<T,T>::Pointer orienter = itk::OrientImageFilter<T,T>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput( output );
    orienter->Update();
	output = orienter->GetOutput();
	*/
	
	
    return output;
}

void SmoothPartialVector(ImageVectorType::Pointer pv, ByteImageType::Pointer chamfer_colon, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex)
{
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,pv->GetLargestPossibleRegion());

	// Find edge of chamfer colon
	//ByteImageType::Pointer chamfer_colon_edge = FindBinaryEdge(chamfer_colon,startIndex,endIndex);
	//IteratorTypeByteWithIndex chamfer_colon_edge_iter(chamfer_colon_edge,pv->GetLargestPossibleRegion());
	//WriteITK(chamfer_colon_edge,"chamfer_colon_edge.nii");

	IteratorImageVectorType pv_iter(pv,pv->GetLargestPossibleRegion());
	
	ImageType::Pointer part = AllocateNewImage(pv->GetLargestPossibleRegion());
	part->SetSpacing(pv->GetSpacing());
	IteratorTypeFloat4WithIndex part_iter(part,part->GetLargestPossibleRegion());

	
	
	for (int i=0; i<3; i++)
	{
		// Ensure iters are set properly
		pv_iter=IteratorImageVectorType(pv,pv->GetLargestPossibleRegion());
		part_iter=IteratorTypeFloat4WithIndex(part,part->GetLargestPossibleRegion());

		// Extract each dimension of pv vector image
		for (pv_iter.GoToBegin(), part_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
			 !pv_iter.IsAtEnd() && !part_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
			 ++pv_iter, ++part_iter, ++chamfer_colon_iter)
		{
			if (chamfer_colon_iter.Get() == 1)
			{
				CovariantVectorType p = pv_iter.Get();
				part_iter.Set(p[i]);
			} else {
				part_iter.Set(0);
			}
		}

		//WriteITK(part,"part.nii");

		// Apply smoothing filter to one dimension
		DiscreteGaussianFilterType::Pointer dgFilter = DiscreteGaussianFilterType::New();
		dgFilter->SetInput(part);

		std::cout << "Smoothing partial by " << pv->GetSpacing()[0] << "mm" << std::endl;

		dgFilter->SetVariance(pv->GetSpacing()[0]*pv->GetSpacing()[0]);
		dgFilter->Update();
		part=dgFilter->GetOutput();
		part_iter=IteratorTypeFloat4WithIndex(part,part->GetLargestPossibleRegion());

		//WriteITK(part,"part_smooth.nii");

		// Set smoothed value back to partial vector image
		for (pv_iter.GoToBegin(), part_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();//chamfer_colon_edge_iter.GoToBegin();
			!pv_iter.IsAtEnd() && !part_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();// && !chamfer_colon_edge_iter.IsAtEnd();
			 ++pv_iter, ++part_iter, ++chamfer_colon_iter)//, ++chamfer_colon_edge_iter)
		{
			if (chamfer_colon_iter.Get() == 1 /*&& chamfer_colon_edge_iter.Get()==0*/) // in mask but not on edge
			{
				CovariantVectorType p = pv_iter.Get();

				if ( p[i] > 0 && p[i] < 1 )	 // only smooth uncertain partials
				{
					p[i] = part_iter.Get();
					pv_iter.Set(p);
				}
			}
		}
	}
}



ImageType::Pointer ResampleImage(ImageType::Pointer input, ImageType::SpacingType output_spacing)
{
	ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
	
	//typedef itk::AffineTransform< double, 3> TransformType;
	//TransformType::Pointer transform = TransformType::New();
	//resampleFilter->SetTransform(transform);

	//typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  InterpolatorType;
	//InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//resampleFilter->SetInterpolator( interpolator );

	// Use sinc interpolation
	//typedef itk::WindowedSincInterpolateImageFunction<ImageType,3> SincInterpolatorType;
	//SincInterpolatorType::Pointer interpolator = SincInterpolatorType::New();
	//resampleFilter->SetInterpolator( interpolator );
	
	// Uuse B-spline interpolation
	//InterpolationType::Pointer interpolator = InterpolationType::New();
	//resampleFilter->SetInterpolator( interpolator );


	resampleFilter->SetDefaultPixelValue( -1024 );

	ImageType::SpacingType input_spacing = input->GetSpacing();
	ImageType::SizeType input_size = input->GetLargestPossibleRegion().GetSize();

	// Adjusto 
	ImageType::SizeType output_size;

	output_size[0] = input_size[0] * input_spacing[0] / output_spacing[0];
	output_size[1] = input_size[1] * input_spacing[1] / output_spacing[1];
	output_size[2] = input_size[2] * input_spacing[2] / output_spacing[2];

	resampleFilter->SetOutputSpacing( output_spacing );
	resampleFilter->SetOutputOrigin( input->GetOrigin() );
	resampleFilter->SetOutputDirection( input->GetDirection() );

	resampleFilter->SetSize( output_size );
	resampleFilter->SetInput( input );
	resampleFilter->Update();

	return resampleFilter->GetOutput();
}

bool MatchVoxels(std::string s) 
{
	if ( (s.compare(0,1,"2")==0) && (s.compare(s.size()-1,1,"1")==0) ) // start: air, end: stool
	{
		for (int i=1; i < s.size()-1; i++)
		{
			if ( !( s.compare(i,1,"5")==0 || s.compare(i,1,"6")==0 || s.compare(i,1,"7")==0 || s.compare(i,1,"1")==0 || s.compare(i,1,"2")==0 ) ) // SA, TA, TS
			{
				return false;
			}
		}

		return true;
	}

	return false;

	// "^2*[5-7]+1+$"
	// (vs.compare(i,1,"6") == 0 || vs.compare(i,1,"7") == 0)
	/*
	const int ovecount = 6;
	pcre *re;
	const char *error;
	int erroffset;
	int ovector[ovecount];
	int rc;

	char *data;
	data = new char[inputString.length()+1];
	strcpy(data, inputString.c_str());

	re = pcre_compile(
	regex, // the pattern
	0, // default options
	&error, // for error message
	&erroffset, // for error offset
	NULL); // use default character table
	if (! re)
	{
		fprintf(stderr,"PCRE compilation failed at expression offset%d: %s\n", erroffset, error);
		return 1;
	}

	rc = pcre_exec(
	re, // the compiled pattern 
	NULL, // no extra data - we didn't study the pattern
	data, // the subject string
	strlen(data), // the length of the subject
	0, // start at offset 0 in the subject
	0, // default options 
	ovector, // output vector for substring information
	ovecount); //number of elements in the output vector

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
	*/
}

void VoteVoxels(VoxelTypeImage::Pointer v, ByteImageType::Pointer mask)
{

	/*
	typedef itk::ConstantBoundaryCondition< VoxelTypeImage > ConstantBoundaryConditionType;
	ConstantBoundaryConditionType bc;
	bc.SetConstant( Air );
	*/

	ImageType::RegionType fullRegion = v->GetLargestPossibleRegion();
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);
	
	typedef itk::NeighborhoodIterator< VoxelTypeImage > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	NeighborhoodIteratorType v_iter(radius, v, v->GetLargestPossibleRegion() );
	//v_iter.SetBoundaryCondition(bc);

	IteratorTypeByteWithIndex mask_iter(mask,mask->GetLargestPossibleRegion());
	
	// Create helper image
	VoxelTypeImage::Pointer v2 = VoxelTypeImage::New();
	v2->SetRegions(v->GetLargestPossibleRegion());
	v2->SetSpacing(v->GetSpacing());
	v2->Allocate();
	IteratorTypeVoxelType v2_iter(v2, v2->GetLargestPossibleRegion() );

	
	// Copy image
	for (v_iter.GoToBegin(), v2_iter.GoToBegin(); !v_iter.IsAtEnd() && !v2_iter.IsAtEnd(); ++v_iter, ++v2_iter)
	{
		v2_iter.Set( v_iter.GetCenterPixel() );
	}

	std::ofstream file;
	file.open("vote.txt");

	// Set neighbors based on majority vote
	for (v_iter.GoToBegin(), v2_iter.GoToBegin(), mask_iter.GoToBegin();
		!v_iter.IsAtEnd() && !v2_iter.IsAtEnd() && !mask_iter.IsAtEnd();
		++v_iter, ++v2_iter, ++mask_iter)
	{

		if (mask_iter.Get() == 1) // in colon mask
		{
			VoxelType voxel = v_iter.GetCenterPixel();
			VoxelTypeImage::IndexType index = v_iter.GetIndex();
			if (index[0]==389 && index[1]==511-326 && index[2]==1)
			{
				std::cout<<"test"<<std::endl;
			}

			if ( voxel == StoolAir || voxel == TissueAir || voxel == TissueStool ) // only on transition voxels
			{
				int countType[8] = {0,0,0,0,0,0,0,0};

				for (int i=0; i<26; i++)
				{
					VoxelType neighbor = v_iter.GetPixel(i);

					int j = VoxelTypeToNum(neighbor);
					++countType[j-1];

					if ( neighbor != voxel && countType[j-1] >= 13 ) // if neighbors are different and majority
					{
						VoxelTypeImage::IndexType idx = v_iter.GetIndex();

						if (idx[0] >= startIndex[0] && idx[0] <= endIndex[0] && idx[1] >= startIndex[1] && idx[1] <= endIndex[1] && idx[2] >= startIndex[2] && idx[2] <= endIndex[2])
						{
							v2->SetPixel( idx , NumToVoxelType(j) );
							file << idx[0]<<"\t"<<511-idx[1]<<"\t"<<idx[2]<<"\n";
							break;
						} else {
							std::cout << "out of bounds" << std::endl;
						}

						
					}
				}
			}
		}
	}

	file.close();

	// Copy image
	for (v_iter.GoToBegin(), v2_iter.GoToBegin(), mask_iter.GoToBegin();
		!v_iter.IsAtEnd() && !v2_iter.IsAtEnd() && !mask_iter.IsAtEnd();
		++v_iter, ++v2_iter, ++mask_iter)
	{
		v_iter.SetCenterPixel( v2_iter.Get() );
	}
}

ImageType::Pointer ComputeNeighborhoodSmax(ImageType::Pointer input, VoxelTypeImage::Pointer v, IteratorTypeByteWithIndex &mask_iter, ImageType::IndexType &sidx, ImageType::IndexType &eidx)
{
	ImageType::Pointer smax = AllocateNewImage(input->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex smax_iter(smax,input->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_iter(input,input->GetLargestPossibleRegion());

	NeighborhoodIteratorVoxelType::RadiusType radius;

	radius.Fill(2);

	NeighborhoodIteratorVoxelType nit(radius, v, v->GetLargestPossibleRegion());

	VoxelTypeImage::SizeType size = nit.GetSize();
	int n = size[0]*size[1]*size[2] - 1;

	nit.GoToBegin();
	mask_iter.GoToBegin();
	smax_iter.GoToBegin();
	input_iter.GoToBegin();

	float max=0;

	while (!nit.IsAtEnd())
	{
		if (mask_iter.Get() == 1 && nit.GetCenterPixel()==Unclassified)
		{
			max=0;

			for (int i=0; i<n;i++)
			{
				//if (i != n/2 )
				//{
					VoxelTypeImage::IndexType idx = nit.GetIndex(i);
					
					if (idx[0] >= sidx[0] && idx[0] <= eidx[0] && idx[1] >= sidx[1] && idx[1] <= eidx[1] && idx[2] >= sidx[2] && idx[2] <= eidx[2])
					{
						float val = input->GetPixel(idx);
						
						if (nit.GetPixel(i) == Stool /*|| val > N_SVAL*/)
						{

							if ( val > max )
							{
								max = val;
							}	
						}
						
					//}

					
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

void WritePartialImages(ImageVectorType::Pointer partialVector, ByteImageType::Pointer chamfer_colon, std::string name)
{
	ImageType::RegionType fullRegion = partialVector->GetLargestPossibleRegion();

	ImageType::Pointer partial_image = AllocateNewImage(fullRegion);

	IteratorImageVectorType partialVector_iter(partialVector,fullRegion);
	IteratorTypeFloat4WithIndex partial_image_iter(partial_image,fullRegion);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);

	for (int i=0; i < 3; i++)
	{
		for (partialVector_iter.GoToBegin(), partial_image_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
			!partialVector_iter.IsAtEnd() && !partial_image_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
			++partialVector_iter, ++partial_image_iter, ++chamfer_colon_iter)
		{
			if (chamfer_colon_iter.Get() == 1)
			{
				CovariantVectorType partial = partialVector_iter.Get();
				partial_image_iter.Set( partial[i] );
			} else {
				partial_image_iter.Set(0);
			}
		}
		std::stringstream ss;

		ss << "partial_" << i << "_" << name << ".nii";
		
		if (i==1)
		{
			WriteITK(partial_image, ss.str() );
		}
	}
}

vnl_matrix_fixed< float, 3, 3> EvaluateAtNeighborhood (const vnl_matrix_fixed<float,3,3> Next, const vnl_matrix_fixed<float,3,3> Previous, float weights[3])
{
	unsigned i, j;
	vnl_matrix_fixed<float,3,3> J;

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			J[i][j] = weights[i]* 0.5 * (Next[i][j] - Previous[i][j]);
		}
	}

	return J;
}

ByteImageType::Pointer FindBinaryEdge(ByteImageType::Pointer im, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex)
{
	IteratorTypeByteWithIndex im_iter(im,im->GetLargestPossibleRegion());

	ByteImageType::Pointer edge = AllocateNewByteImage(im->GetLargestPossibleRegion());
	edge->FillBuffer(0);

	IteratorTypeByteWithIndex edge_iter(edge,im->GetLargestPossibleRegion());

	int radius = 2;

	for (im_iter.GoToBegin(); !im_iter.IsAtEnd(); ++im_iter)
	{

		if ( im_iter.Get() == 1)
		{
			bool isEdge = false;

			ImageType::IndexType index = im_iter.GetIndex();
			for(int i=-radius;i<=radius;i++) {
				if (!isEdge && index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
					for (int j=-radius;j<=radius;j++) {
						if (!isEdge && index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
							for (int k=-radius;k<=radius;k++) {
								if (!isEdge && index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
									
									ImageType::IndexType neighborIndex={index[0]+i,index[1]+j,index[2]+k};
									if ( im->GetPixel( neighborIndex ) == 0 )
									{
										edge->SetPixel( index, 1);
										isEdge = true;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return edge;
}


void RunFuzzy(ImageType::Pointer input, ByteImageType::Pointer chamfer_colon, VoxelTypeImage::Pointer voxel_type)
{
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,chamfer_colon->GetLargestPossibleRegion());
	IteratorTypeVoxelType voxel_type_iter(voxel_type,voxel_type->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_iter(input,input->GetLargestPossibleRegion());

	// Create smaller input region
	ImageType::Pointer input_sub = AllocateNewImage(input->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_sub_iter(input_sub,input->GetLargestPossibleRegion());

	for (input_sub_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
		!input_sub_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
		++input_sub_iter, ++input_iter, ++chamfer_colon_iter)
	{
		if ( chamfer_colon_iter.Get() == 1) {
			input_sub_iter.Set(input_iter.Get());
		} else {
			input_sub_iter.Set(-1025);
		}
	}

	// Calculate mean and variance of tissue voxels
	float mean = 0;

	// Choose variance such that max tissue near stool (TNS) is within +/- 3 std dev

	/*
	float sum = 0;
	int count = 0;

	for (voxel_type_iter.GoToBegin(), input_iter.GoToBegin(); !voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd(); ++voxel_type_iter, ++input_iter)
	{
		if (voxel_type_iter.Get() == Tissue)
		{
			sum += input_iter.Get();
			count++;
		}
	}

	float mean = sum/count;

	float variance = 0;

	for (voxel_type_iter.GoToBegin(), input_iter.GoToBegin(); !voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd(); ++voxel_type_iter, ++input_iter)
	{
		if (voxel_type_iter.Get() == Tissue)
		{
			variance += vnl_math_sqr(input_iter.Get()-mean);
		}
	}
	
	variance = variance/count;
	*/

	float Tmax[2] = {800,900};
	float variance[2];

	for (int i=0; i<2; i++)
	{
		variance[i] = Tmax[i]*Tmax[i]/9;

		std::cout << "Fuzzy variance: " << variance[i] << std::endl;

		FuzzyFilterType::Pointer fuzzyFilter = FuzzyFilterType::New();
		fuzzyFilter->SetInput(input_sub);
		fuzzyFilter->SetOutsideValue(0);
		fuzzyFilter->SetInsideValue(1);
		fuzzyFilter->SetThreshold(0.1);

		fuzzyFilter->SetMean(mean);
		fuzzyFilter->SetVariance(variance[i]);

		input_iter.GoToBegin();
		chamfer_colon_iter.GoToBegin();
		voxel_type_iter.GoToBegin();

		ImageType::IndexType idx;

		while (!input_iter.IsAtEnd())
		{

			if ( chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get() == Tissue )
			{
				idx = input_iter.GetIndex();
				fuzzyFilter->SetObjectSeed(idx);
				break;
			}

			++input_iter;
			++chamfer_colon_iter;
			++voxel_type_iter;
		}

		fuzzyFilter->Update();

		typedef itk::CastImageFilter< FuzzySceneType, ImageType > CastFilterType;
		CastFilterType::Pointer caster = CastFilterType::New();
		caster->SetInput(fuzzyFilter->GetFuzzyScene());
		
		typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescaleFilterType;
		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
		rescaler->SetInput(caster->GetOutput());
		rescaler->SetOutputMinimum(0);
		rescaler->SetOutputMaximum(1);
		rescaler->Update();

		std::stringstream ss;
		ss << "fuzzy_scene_mean" << mean << "_variance_" << variance[i] << "_Tmax_" << Tmax[i] << ".nii";

		WriteITK(rescaler->GetOutput(), ss.str());
	}

	/*

	typedef itk::ImageFileWriter< FuzzySceneType > FuzzyWriterType;
	FuzzyWriterType::Pointer fwriter = FuzzyWriterType::New();
	fwriter->SetFileName("fuzzy_scene.nii");
	fwriter->SetInput(fuzzyFilter->GetFuzzyScene());
	fwriter->Update();
	*/
}

void RunConnectedThresholdGrowing(ImageType::Pointer input, ByteImageType::Pointer chamfer_colon, VoxelTypeImage::Pointer voxel_type)
{
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,chamfer_colon->GetLargestPossibleRegion());
	IteratorTypeVoxelType voxel_type_iter(voxel_type,voxel_type->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_iter(input,input->GetLargestPossibleRegion());

	// Create smaller input region
	ImageType::Pointer input_sub = AllocateNewImage(input->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_sub_iter(input_sub,input->GetLargestPossibleRegion());

	for (input_sub_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
		!input_sub_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
		++input_sub_iter, ++input_iter, ++chamfer_colon_iter)
	{
		if ( chamfer_colon_iter.Get() == 1) {
			input_sub_iter.Set(input_iter.Get());
		} else {
			input_sub_iter.Set(-1025);
		}
	}

	// Smooth image
	typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> GradientAnisotropicDiffusionFilterType;
	GradientAnisotropicDiffusionFilterType::Pointer smoother = GradientAnisotropicDiffusionFilterType::New();
	smoother->SetNumberOfIterations(5);
	smoother->SetTimeStep(0.125);
	smoother->SetInput(input_sub);
	smoother->Update();
	WriteITK(smoother->GetOutput(),"input_sub_smoothed.nii");

	float lower = -300;
	float upper = 600;

	typedef itk::ConnectedThresholdImageFilter<ImageType,ImageType> ConnectedThresholdFilterType;
	ConnectedThresholdFilterType::Pointer connecter = ConnectedThresholdFilterType::New();
	connecter->SetInput(smoother->GetOutput());
	connecter->SetLower(lower);
	connecter->SetUpper(upper);
	connecter->SetReplaceValue(1);

	

	input_iter.GoToBegin();
	chamfer_colon_iter.GoToBegin();
	voxel_type_iter.GoToBegin();

	int count = 0;

	while (!input_iter.IsAtEnd())
	{

		if ( chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get() == Tissue )
		{
			ImageType::IndexType idx = input_iter.GetIndex();
			
			if (count == 0)
			{
				connecter->SetSeed(idx);
			} else {
				connecter->AddSeed(idx);
			}

			count++;
		}

		++input_iter;
		++chamfer_colon_iter;
		++voxel_type_iter;
	}

	connecter->Update();

	std::stringstream ss;
	ss << "connected_threshold_lower" << lower << "_upper_" << upper << ".nii";
	WriteITK(connecter->GetOutput(),ss.str());

}

/********************************************************************

Remove all components fully surrounded by stool. 
Assumes that all non-stool, non-holes are connected.

*********************************************************************/
void HeuristicClosing(VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon)
{
	ImageType::RegionType region = voxel_type->GetLargestPossibleRegion();

	ByteImageType::Pointer mask = AllocateNewByteImage(region);

	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);
	IteratorTypeByteWithIndex mask_iter(mask,region);

	// Get non stool mask
	for (voxel_type_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), mask_iter.GoToBegin(); !mask_iter.IsAtEnd(); ++voxel_type_iter, ++chamfer_colon_iter, ++mask_iter)
	{
		if (chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get() == Stool)
		{
			mask_iter.Set(0);
		} else {
			mask_iter.Set(1);
		}
	}

	WriteITK(mask,"non_stool_mask.nii");

	// Apply closing slice by slice
	
	typedef itk::Image< unsigned char, 2> ByteImageType2D;
	typedef itk::Image< int, 2>			  IntImageType2D;

	typedef itk::ConnectedComponentImageFilter< IntImageType2D, IntImageType2D > ConnectedComponentFilterType2D;
	typedef itk::RelabelComponentImageFilter< IntImageType2D, IntImageType2D> RelabelFilterType2D;
	typedef itk::BinaryThresholdImageFilter< IntImageType2D, ByteImageType2D> BinaryThresholdType2D;
	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, ConnectedComponentFilterType2D, BinaryThresholdType2D> SliceBySliceFilterType;

	SliceBySliceFilterType::Pointer sliceFilter = SliceBySliceFilterType::New();

	ConnectedComponentFilterType2D::Pointer ccFilter = ConnectedComponentFilterType2D::New();
	ccFilter->FullyConnectedOn();

	RelabelFilterType2D::Pointer relabelFilter = RelabelFilterType2D::New();
	relabelFilter->SetInput( ccFilter->GetOutput() );

	BinaryThresholdType2D::Pointer binaryFilter = BinaryThresholdType2D::New();
	binaryFilter->SetInput ( relabelFilter->GetOutput() );
	binaryFilter->SetLowerThreshold(1);
	binaryFilter->SetUpperThreshold(1);
	binaryFilter->SetInsideValue(0);
	binaryFilter->SetOutsideValue(1);

	sliceFilter->SetInput( mask );
	sliceFilter->SetInputFilter( ccFilter );
	sliceFilter->SetOutputFilter( binaryFilter );
	sliceFilter->Update();

	ByteImageType::Pointer mask_relabel = sliceFilter->GetOutput();

	WriteITK(mask_relabel,"mask_relabel.nii");

	IteratorTypeByteWithIndex mask_relabel_iter(mask_relabel,region);

	for (voxel_type_iter.GoToBegin(), mask_relabel_iter.GoToBegin(); !voxel_type_iter.IsAtEnd(); ++voxel_type_iter, ++mask_relabel_iter)
	{
		if ( mask_relabel_iter.Get() == 1 )
		{
			voxel_type_iter.Set( Stool );
		}
	}
	
	/*

	// Run connected component
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->FullyConnectedOn(); // use face-edge-vertex
	ccFilter->SetInput(mask);
	ccFilter->Update();
	IntImageType::Pointer cc = ccFilter->GetOutput();
	IteratorTypeIntWithIndex cc_iter(cc,region);
	
	ccFilter.~SmartPointer();
	WriteITK(cc,"mask_cc.nii");

	// Convert to Label map
	LabelImageToShapeLabelMapFilterType::Pointer converter = LabelImageToShapeLabelMapFilterType::New();
	converter->SetInput(cc);
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();
	converter.~SmartPointer();

	// Find largest non-mask component
	std::vector<LabelObjectType::Pointer> labelVector = labelMap->GetLabelObjects();
	sort(labelVector.begin(), labelVector.end(), compareSize);
	
	for (voxel_type_iter.GoToBegin(), cc_iter.GoToBegin(); !voxel_type_iter.IsAtEnd(); ++voxel_type_iter, ++cc_iter)
	{
		if (cc_iter.Get() != 0 && cc_iter.Get() != labelVector[0]->GetLabel())
		{
			voxel_type_iter.Set(Stool);
		}
	}

	*/
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
}

ByteImageType::Pointer SegmentColon(ImageType::Pointer input)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	IteratorTypeFloat4WithIndex input_iter(input,region);

	// Find air lumen
	ByteImageType::Pointer air_mask = AllocateNewByteImage(region);
	IteratorTypeByteWithIndex air_mask_iter(air_mask,region);

	for (input_iter.GoToBegin(), air_mask_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_mask_iter)
	{
		if (input_iter.Get() < -600)
		{
			air_mask_iter.Set(1);
		} else {
			air_mask_iter.Set(0);
		}
	}

	WriteITK(air_mask, "air_lumen.nii" );

	// Apply median filter to remove noise
	typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType > BinaryMedianFilterType;
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput( air_mask );
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 1 );

	ByteImageType::SizeType radius;
	radius[0] = 4;
	radius[1] = 4;
	radius[2] = 0;

	medianFilter->SetRadius( radius );
	medianFilter->Update();
	air_mask = medianFilter->GetOutput();
	air_mask_iter = IteratorTypeByteWithIndex(air_mask,region);

	WriteITK(air_mask, "air_mask_median.nii");

	// Remove air component with most pixels on image border
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput( air_mask );
	ccFilter->Update();
	IntImageType::Pointer cc = ccFilter->GetOutput();
	ccFilter.~SmartPointer();

	WriteITK(cc, "air_cc.nii");

	// Relabel so that background comp is label 1
	typedef itk::ShapeRelabelImageFilter< IntImageType > ShapeRelabelImageFilterType;
	ShapeRelabelImageFilterType::Pointer relabeler = ShapeRelabelImageFilterType::New();
	relabeler->SetInput( cc );
	relabeler->SetAttribute("SizeOnBorder");
	relabeler->SetBackgroundValue(0);
	//relabeler->SetReverseOrdering(true);
	relabeler->Update();
	cc = relabeler->GetOutput();
	
	WriteITK(cc,"air_cc_relabel_sizeOnBorder.nii");

	IteratorTypeIntWithIndex cc_iter(cc,region);

	for(air_mask_iter.GoToBegin(), cc_iter.GoToBegin(); !air_mask_iter.IsAtEnd(); ++air_mask_iter, ++cc_iter)
	{
		if (cc_iter.Get() < 2)
		{
			air_mask_iter.Set(0);
		}
	}

	relabeler.~SmartPointer();
	cc.~SmartPointer();

	WriteITK(air_mask,"air_lumen_noBkg.nii");


	/*typedef itk::LabelShapeOpeningImageFilter< IntImageType > LabelShapeOpeningImageFilterType;
	LabelShapeOpeningImageFilterType::Pointer labelFilter = LabelShapeOpeningImageFilterType::New();

	labelFilter->SetInput( cc );
	labelFilter->SetAttribute("SizeOnBorder");
	labelFilter->SetBackgroundValue(0);
	labelFilter->SetLambda(512);
	labelFilter->SetReverseOrdering(true);
	labelFilter->Update();
	cc = labelFilter->GetOutput();

	WriteITK(cc,"air_cc_no_bkg_component.nii");
	*/



	// Dilate air mask

	// Create dilate filter
	typedef itk::BinaryDilateImageFilter< ByteImageType, ByteImageType, StructuringElementType> DilateFilterType;

	DilateFilterType::Pointer dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(1);
	dilater->SetInput( air_mask );

	// Create structuring element
	StructuringElementType  structuringElement;

	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;

    structuringElement.SetRadius( radius );
    structuringElement.CreateStructuringElement();

	dilater->SetKernel( structuringElement );
	dilater->Update();
	ByteImageType::Pointer air_mask_dilated = dilater->GetOutput();
	IteratorTypeByteWithIndex air_mask_dilated_iter(air_mask_dilated,region);

	WriteITK(air_mask_dilated, "air_dilated.nii");

	// Subtract original air_mask area from dilation
	for (air_mask_iter.GoToBegin(), air_mask_dilated_iter.GoToBegin(); !air_mask_iter.IsAtEnd(); ++air_mask_iter, ++air_mask_dilated_iter)
	{
		air_mask_dilated_iter.Set( air_mask_dilated_iter.Get() - air_mask_iter.Get() );
	}

	WriteITK(air_mask_dilated,"air_dilated_only.nii");

	// Setup region growing filter
	typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType> ConnectedThresholdImageFilterType;
	ConnectedThresholdImageFilterType::Pointer grower = ConnectedThresholdImageFilterType::New();
	grower->SetInput( input );
	grower->SetLower(200);
	grower->SetUpper(2000);
	grower->SetReplaceValue(1);

	// Set tagged regions within dilated area as seeds of region growing
	for (input_iter.GoToBegin(), air_mask_dilated_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_mask_dilated_iter)
	{
		if (air_mask_dilated_iter.Get() == 1 && input_iter.Get() >= 200)
		{
			grower->AddSeed( input_iter.GetIndex() );
		}
	}

	grower->Update();
	ByteImageType::Pointer tagged = grower->GetOutput();
	IteratorTypeByteWithIndex tagged_iter(tagged,region);

	WriteITK(tagged,"tagged.nii");

	// Add air lumen to tagged region
	for (tagged_iter.GoToBegin(), air_mask_iter.GoToBegin(); !tagged_iter.IsAtEnd(); ++tagged_iter, ++air_mask_iter)
	{
		if ( air_mask_iter.Get() == 1 )
			tagged_iter.Set(1);
	}

	WriteITK(tagged,"tagged_air.nii");

	// Dilate
	dilater = DilateFilterType::New();
	dilater->SetBackgroundValue(0);
	dilater->SetForegroundValue(1);
	dilater->SetInput( tagged );
	
	radius[0] = 8;
	radius[1] = 8;
	radius[2] = 0;

	structuringElement.SetRadius( radius );
	structuringElement.CreateStructuringElement();

	dilater->SetKernel( structuringElement );
	dilater->Update();

	WriteITK(dilater->GetOutput(),"colon.nii");

	return dilater->GetOutput();
}

/********************************************************************

Remove stool components with small size that are far away from air

*********************************************************************/
void CleanIsolatedStool(VoxelTypeImage::Pointer voxel_type)
{	
	// Calculuate distance to air
	ImageType::RegionType region = voxel_type->GetLargestPossibleRegion();
	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);

	ByteImageType::Pointer air_mask = AllocateNewByteImage(region);
	IteratorTypeByteWithIndex air_mask_iter(air_mask,region);

	// Find small stool objects
	ByteImageType::Pointer stool_mask = AllocateNewByteImage(region);
	IteratorTypeByteWithIndex stool_mask_iter(stool_mask,region);

	for (air_mask_iter.GoToBegin(), voxel_type_iter.GoToBegin(), stool_mask_iter.GoToBegin(); !air_mask_iter.IsAtEnd(); ++air_mask_iter, ++voxel_type_iter, ++stool_mask_iter)
	{
		if (voxel_type_iter.Get() == Air)
		{
			air_mask_iter.Set(1);
		} else {
			air_mask_iter.Set(0);
		}

		if (voxel_type_iter.Get() == Stool)
		{
			stool_mask_iter.Set(1);
		} else {
			stool_mask_iter.Set(0);
		}
	}

	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput( air_mask );
	int chamfer_weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();

	air_mask = chamfer_filter->GetOutput();
	air_mask_iter = IteratorTypeByteWithIndex(air_mask,region);
	chamfer_filter.~SmartPointer();

	//WriteITK(air_mask,"distance_to_air.nii");

	// Isolate small stool objects
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput( stool_mask );
	ccFilter->Update();
	IntImageType::Pointer cc = ccFilter->GetOutput();
	ccFilter.~SmartPointer();

	//WriteITK(cc,"stool_cc.nii");

	typedef itk::LabelShapeOpeningImageFilter< IntImageType > LabelShapeOpeningImageFilterType;
	LabelShapeOpeningImageFilterType::Pointer labelFilter = LabelShapeOpeningImageFilterType::New();
	labelFilter->SetInput( cc );
	labelFilter->SetAttribute("Size");
	labelFilter->SetBackgroundValue(0);
	labelFilter->SetLambda(20);
	labelFilter->SetReverseOrdering(true);
	labelFilter->Update();
	cc = labelFilter->GetOutput();
	IteratorTypeIntWithIndex cc_iter(cc,region);

	//WriteITK(cc,"stool_cc_small_only.nii");

	unsigned int count = 0;

	for (cc_iter.GoToBegin(), air_mask_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !cc_iter.IsAtEnd(); ++cc_iter, ++air_mask_iter, ++voxel_type_iter)
	{
		if (cc_iter.Get() != 0 && air_mask_iter.Get() >= 15)
		{
			voxel_type_iter.Set(Tissue);
			count++;
		}
	}

	std::cout << "Number of stool voxels cleaned: " << count << std::endl;
}

ImageType::Pointer Sharpen(ImageType::Pointer input)
{
	// Copy input to output
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	ImageType::Pointer output = AllocateNewImage(region);
	output->SetRegions(region);
	output->SetSpacing(input->GetSpacing());
	output->CopyInformation(input);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);

	for (input_iter.GoToBegin(), output_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter)
	{
		output_iter.Set( input_iter.Get() );
	}

	NeighborhoodIteratorType::RadiusType radius;

	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;

	NeighborhoodIteratorType nit(radius, input, input->GetLargestPossibleRegion());

	ImageType::SizeType size = nit.GetSize();
	int n = size[0]*size[1]*size[2] - 1;

	nit.GoToBegin();
	output_iter.GoToBegin();

	while (!nit.IsAtEnd())
	{
		float sum = 0;

		for (int i=0; i<9; i++)
		{
			if ( i != 4) //center pixel
			{
				sum += (-1 * nit.GetPixel(i));
			}
		}

		float val = ( 12*nit.GetCenterPixel() + sum ) / 4;

		output_iter.Set( val );

		++nit;
		++output_iter;
	}

	return output;
}

void DiffusionTest( ImageType::Pointer input, double ContrastParameter, double sigma)
{
	typedef itk::Image< double, 3 > DoubleImageType;

	// Cast image to type double
	typedef itk::CastImageFilter<ImageType, DoubleImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput( input );
	caster->Update();
	
	DoubleImageType::Pointer input_double = caster->GetOutput();

	// Run diffusion filter
	/*typedef itk::AnisotropicHybridDiffusionImageFilter< DoubleImageType, DoubleImageType> HybridFilterType;
	HybridFilterType::Pointer hybridFilter = HybridFilterType::New();
	hybridFilter->SetInput ( input_double );*/

	typedef itk::AnisotropicEdgeEnhancementDiffusionImageFilter<DoubleImageType, DoubleImageType> EdgeDiffusionFilterType;

	typedef itk::AnisotropicCoherenceEnhancingDiffusionImageFilter<DoubleImageType, DoubleImageType> CoherenceDiffusionFilterType;
	
	CoherenceDiffusionFilterType::Pointer diffusionFilter = CoherenceDiffusionFilterType::New();
	diffusionFilter->SetInput( input_double );

	diffusionFilter->SetContrastParameterLambdaC( ContrastParameter );
	diffusionFilter->SetSigma( sigma );

	// Cast back to float
	typedef itk::CastImageFilter<DoubleImageType, ImageType> CastImageFilterType2;
	CastImageFilterType2::Pointer caster2 = CastImageFilterType2::New();
	caster2->SetInput( diffusionFilter->GetOutput() );
	caster2->Update();

	std::stringstream ss;
	ss << "diffusion_coherence_contrast_" << ContrastParameter << "_sigma_" << sigma << ".nii";
	WriteITK(caster2->GetOutput(), ss.str());

}

void SubtractStool( ImageType::Pointer input, VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Allocate subtracted output image
	ImageType::Pointer output = AllocateNewImage(region);
	output->CopyInformation(input);
	output->SetSpacing(input->GetSpacing());

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);
	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);

	for (input_iter.GoToBegin(), output_iter.GoToBegin(), voxel_type_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter, ++voxel_type_iter, ++chamfer_colon_iter)
	{
		output_iter.Set( input_iter.Get() );

		if ( chamfer_colon_iter.Get() == 1 )
		{
			if ( voxel_type_iter.Get() == Stool || voxel_type_iter.Get() == Air || voxel_type_iter.Get() == TissueStool || voxel_type_iter.Get() == StoolAir || input_iter.Get() > 150)
			{
				output_iter.Set(-1025);
			}
		}
	}

	WriteITK(output,"basic_subtraction.nii");

	/*

	// Threshold to detect non-air within mask
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing( input->GetSpacing() );
	air->Allocate();

	IteratorTypeByteWithIndex air_iter(air,region);

	for (input_iter.GoToBegin(), air_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++air_iter, ++chamfer_colon_iter)
	{
		if (chamfer_colon_iter.Get()==1 && input_iter.Get() < -300) // Chosen to reconstruct entire mucosa (including regions with no stool) for uniformity
		{
			air_iter.Set( 0 );
		} else {
			air_iter.Set( 1 );
		}
	}

	// Remove non-air components that are not connected to the largest body of tissue
	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType > BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryFilter->SetInput( air );
	binaryFilter->SetBackgroundValue(0);
	binaryFilter->SetForegroundValue(1);
	binaryFilter->SetAttribute("Size");
	binaryFilter->SetNumberOfObjects(1);
	binaryFilter->Update();
	air = binaryFilter->GetOutput();
	air_iter = IteratorTypeByteWithIndex(air,region);

	for (output_iter.GoToBegin(), air_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter, ++air_iter, ++chamfer_colon_iter)
	{
		if (chamfer_colon_iter.Get()==1 && air_iter.Get() == 0)
		{
			output_iter.Set(-1025);
		}
	}

	WriteITK(output,"full_subtraction_connected.nii");

	// Apply mask to output
	typedef itk::MaskImageFilter< ImageType, ByteImageType> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( output );
	masker->SetInput2( chamfer_colon );
	masker->Update();
	output = masker->GetOutput();

	WriteITK(output,"full_subtraction_post_mask.nii");

	*/

	// Reconstruct mucosa slice by slice
	typedef itk::Image< PixelType, 2 > ImageType2D;

	typedef itk::MucosalReconstructionFilter< ImageType2D > MucosalReconstructionFilterType2D;
	typedef itk::SliceBySliceImageFilter< ImageType, ImageType, MucosalReconstructionFilterType2D, MucosalReconstructionFilterType2D > SliceBySliceFilterType;

	MucosalReconstructionFilterType2D::Pointer reconstructionFilter = MucosalReconstructionFilterType2D::New();
	SliceBySliceFilterType::Pointer sliceFilter = SliceBySliceFilterType::New();
	sliceFilter->SetInput( output );
	sliceFilter->SetFilter( reconstructionFilter );
	sliceFilter->Update();
	
	WriteITK(sliceFilter->GetOutput(),"output_reconstruction.nii");
}

void RemoveStool4(ImageType::Pointer input, ByteImageType::Pointer colon)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();	

	VoxelTypeImage::Pointer vmap = VoxelTypeImage::New();
	vmap->SetRegions(region);
	vmap->SetSpacing(input->GetSpacing());
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeByteWithIndex colon_iter(colon,region);
	IteratorTypeVoxelType vmap_iter(vmap,region);

	// Find optimal tissue stool threshold

	// Mask input with colon segmentation
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( input );
	masker->SetInput2( colon );
	masker->SetOutsideValue( -1500 ); // arbitrary
	masker->Update();
	ImageType::Pointer input_mask = masker->GetOutput();

	// Compute otsu threshold for tissue and stool intensity
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input_mask );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-250);
	otsuCalculator->SetHistogramMax(1400);
	otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	otsuCalculator->Compute();

	float tissue_stool_threshold = otsuCalculator->GetThreshold();

	otsuCalculator.~SmartPointer();
	input_mask.~SmartPointer();
	masker.~SmartPointer();
	
	std::cout << std::endl << std::endl;
	std::cout << "Tissue Stool Otsu Threshold: " << tissue_stool_threshold << std::endl;
	std::cout << std::endl << std::endl;

	// Find air and stool

	for (input_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++vmap_iter, ++colon_iter)
	{
		if (colon_iter.Get() == 1)
		{
			if (input_iter.Get() < -600 )
			{

				vmap_iter.Set( Air );

			} else if ( input_iter.Get() > tissue_stool_threshold ) {
				
				vmap_iter.Set( Stool );
			}
		}
	}

	WriteITK(vmap,"vmap.nii");

}

void RemoveStool5(ImageType::Pointer input, ByteImageType::Pointer colon)
{
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();

	IteratorTypeFloat4WithIndex input_iter(input,fullRegion);
	IteratorTypeByteWithIndex colon_iter(colon,fullRegion);

	
	WriteITK(colon, "colon.nii");

	// Shift input so that air is 0 HU
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage(input);
	minMaxCalc->ComputeMinimum();
	PixelType min = minMaxCalc->GetMinimum();

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		input_iter.Set( input_iter.Get() + abs(min) );
	}
	minMaxCalc.~SmartPointer();


	// Mask with colon
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(input);
	masker->SetInput2(colon);
	masker->SetOutsideValue(0); // arbitrary
	masker->Update();
	input=masker->GetOutput();
	masker.~SmartPointer();

	WriteITK(input, "input_masked.nii");

	// Otsu threshold intensity
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage( input );
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(750);
	otsuCalculator->SetHistogramMax(2400);
	otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	otsuCalculator->Compute();
	PixelType tissue_stool_threshold = otsuCalculator->GetThreshold();
	otsuCalculator.~SmartPointer();
	
	std::cout << std::endl << std::endl;
	std::cout << "Tissue Stool Otsu Threshold: " << tissue_stool_threshold << std::endl;
	std::cout << std::endl << std::endl;

	// Get gradient using canny, supress non-max
	CannyEdgeDetectionImageFilterType::Pointer canny = CannyEdgeDetectionImageFilterType::New();
	canny->SetInput(input);
	canny->SetLowerThreshold(100); // minimum gradient threshold
	canny->SetVariance(1);
	canny->Update();
	ImageType::Pointer grad = canny->GetOutput();

	WriteITK(grad,"canny.nii");

	// Mask gradient with colon (after 1px erosion) to remove outside region
	BinaryErodeImageFilterType::Pointer eroder = BinaryErodeImageFilterType::New();

	StructuringElementType se;
	ImageType::SizeType radius;
	radius.Fill(0);
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;
	se.SetRadius(radius);
	se.CreateStructuringElement();

	eroder->SetKernel(se);
	eroder->SetInput(colon);
	eroder->SetForegroundValue(1);
	eroder->SetBackgroundValue(0);

	eroder->Update();
	ByteImageType::Pointer colonErode = eroder->GetOutput();

	WriteITK(colonErode,"colonErode.nii");

	masker = MaskImageFilterType::New();
	masker->SetInput1(grad);
	masker->SetInput2(colonErode);
	masker->SetOutsideValue(0);
	masker->Update();
	grad=masker->GetOutput();
	
	eroder.~SmartPointer();
	colonErode.~SmartPointer();
	masker.~SmartPointer();

	WriteITK(grad,"canny_erode.nii");

	//// Cast float to 16-bit int
	//typedef itk::CastImageFilter<ImageType,UshortImageType> CastImageFilterType;
	//CastImageFilterType::Pointer caster = CastImageFilterType::New();
	//caster->SetInput(grad);
	//caster->Update();
	//UshortImageType::Pointer grad2 = caster->GetOutput();

	//WriteITK(grad2,"canny_erode_cast.nii");

	// Get grad max
	minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage(grad);
	minMaxCalc->ComputeMaximum();
	PixelType gradMax = minMaxCalc->GetMaximum();

	// Otsu gradient
	otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(grad);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(0);
	otsuCalculator->SetHistogramMax(gradMax);
	otsuCalculator->SetPrintHistogram(note+"_gradient.csv");
	otsuCalculator->Compute();

	PixelType grad_threshold = otsuCalculator->GetThreshold();

	otsuCalculator.~SmartPointer();
	
	std::cout << std::endl << std::endl;
	std::cout << "Tissue Stool Otsu Threshold: " << grad_threshold << std::endl;
	std::cout << std::endl << std::endl;

	// Regular gradient
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,ImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradMag = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradMag->SetInput(input);
	gradMag->SetSigma(1);
	gradMag->Update();
	WriteITK(gradMag->GetOutput(),"gradMag.nii");

}