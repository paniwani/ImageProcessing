#include "RemoveStool3.h"

float Modified=false;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
double beta=.7;
double weight_sum=2.5;
int writeCount=1;
std::string note = "";

const double ALPHA = 0.5;
const double BETA = 0.3;
const double GAMMA = 0.3;
const double ETA = 0.2;

//boost::xpressive::sregex rex = boost::xpressive::sregex::compile( "^2*[5-7]+1+$" );

int main(int argc, char * argv[])
{
	//-------------------------------------------BEGIN SETUP-------------------------------------------------------

	std::cerr<<"Started"<<std::endl;

	// Use parameter changes
	Modified = true;

	// Read input
	//ImageType::Pointer input = ReadDicom("C:/Documents and Settings/panjwanin/Desktop/Matt_to_Neil/nn1_008_10p.i0373/dcm");	
	ImageType::Pointer input = ImageType::New();
	//ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_86-94.hdr");
	//ReadITK(input, "C:/ImageData/mr10_092_13p.i0344.hdr");
	ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_100-105.hdr");
	//ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_75-125.hdr");
	
	// Set region
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);

	// declares a temporary storage of the image with the selected region
	ImageType::Pointer input_temp = AllocateNewImage(fullRegion);
	input_temp->CopyInformation(input);
	//input_temp->SetSpacing(input->GetSpacing());

	// declares two iterators, one for the original image, one for the sub image
	IteratorTypeFloat4WithIndex input_temp_iter(input,fullRegion);
	IteratorTypeFloat4WithIndex input_iter(input_temp,fullRegion);

	// Copy Image the sub region into the temporary container
	for(input_iter.GoToBegin() , input_temp_iter.GoToBegin() ; !input_iter.IsAtEnd() && !input_temp_iter.IsAtEnd(); ++input_iter, ++input_temp_iter) {
		input_iter.Set(input_temp_iter.Get());
	}

	// Print out the sub region of the image, un modified
	WriteITK(input_temp,"input_temp.hdr");

	//ComputeFrangiHessian(input);
	//ComputeHessianResponse3(input);
	//ComputeSatoHessian(input);

	//HessianMeasure(input);

	// Test resampling
	//ImageType::Pointer input_resample = AllocateNewImage(fullRegion);
	//input_resample = ResampleImage(input);
	//WriteITK(input_resample, "input_resample.hdr");

	// Test hessian response

	//ComputeHessianResponse(input);

	//ImageType::Pointer Himage = AllocateNewImage(fullRegion);
	//Himage = ComputeHessianResponse(input_temp);
	//WriteITK(Himage,"H.hdr");

	// Creates an image to store the type of voxel for each voxel in the sub image
	VoxelTypeImage::Pointer voxel_type = VoxelTypeImage::New();
	voxel_type->SetRegions(fullRegion);
	voxel_type->Allocate();
	IteratorTypeVoxelType voxel_type_iter(voxel_type,fullRegion);

	// Get gradient vectors
    GradientFilterType::Pointer gradient_filter = GradientFilterType::New();
    gradient_filter->SetInput(input_temp);
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

		gradient_magnitude_iter.Set( 2*norm );
	}
	
	//writes out the gradient magnitude for the input image
	WriteITK(gradient_magnitude, "gradient.hdr");

	// Create temporary image and iter
	ImageType::Pointer temp = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex temp_iter(temp,fullRegion);

	//stores the largest gradient magnitude, this is purely for reference on the validity of the image
	float gradient_max=0;

	//Computes an initial hard threshold for the voxel type and sets initial partials
	for(input_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !input_iter.IsAtEnd() && !gradient_magnitude_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd();
        ++input_iter, ++gradient_magnitude_iter, ++voxel_type_iter) 
	{
		//updates for maximum gradient
		if (gradient_max<gradient_magnitude_iter.Get()) { gradient_max=gradient_magnitude_iter.Get(); }

		//computes the hard threshold classification for a voxel in the image: Air/Tissue/Stool/Unclassified
        VoxelType type = SingleMaterialClassification(input_iter.Get(), gradient_magnitude_iter.Get());
		voxel_type_iter.Set(type);		//sets the type in the voxel_type structure
    }

	WriteITK(voxel_type, "voxel_type_start_new.hdr");

	// outputs the gradient magnitude purely for reference
	std::cerr<<"Gradient max: "<<gradient_max<<std::endl;


	// Fill in Stool Holes
	for (voxel_type_iter.GoToBegin(), temp_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++temp_iter)
    {
        switch (voxel_type_iter.Get()) {
            case Tissue:
                temp_iter.Set(0);
                break;
            case Stool:
                temp_iter.Set(1);
                break;
            case Unclassified:
                temp_iter.Set(0);
                break;
            case Air:
                temp_iter.Set(0);
                break;
        }
    }

	WriteITK(temp, "temp_stool.hdr");

	/*// Use voting filter to fill in holes
	HoleFillingFilterType::Pointer holeFilter = HoleFillingFilterType::New();
	holeFilter->SetInput(temp);
	ByteImageType::SizeType hole_radius;
	hole_radius[0] = 1;
	hole_radius[1] = 1;
	hole_radius[2] = 1;
	holeFilter->SetRadius(hole_radius);
	holeFilter->SetBackgroundValue(0);
	holeFilter->SetForegroundValue(1);
	holeFilter->SetMajorityThreshold( 11 );
	holeFilter->SetMaximumNumberOfIterations(1);
	holeFilter->Update();
	temp=holeFilter->GetOutput();
	WriteITK(temp,"stool_voting_hole_filled_11.hdr");
	*/

    StructuringElementType  structuringElement;
    structuringElement.SetRadius( 1 );  // .5
    structuringElement.CreateStructuringElement();
	ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
	closingFilter->SetKernel(structuringElement);	 
	closingFilter->SetInput(temp);
	closingFilter->Update();
	temp=closingFilter->GetOutput();
	temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	closingFilter.~SmartPointer();

	WriteITK(temp, "temp_stool_closed.hdr");
	

	// Apply stool closing and prepare for tissue closing
	for (temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++temp_iter)
    {
        switch(voxel_type_iter.Get()) {
            case Tissue:         
				if (temp_iter.Get()==1) {
					voxel_type_iter.Set(Stool);			//Modify tissue into stool holes
                }
                temp_iter.Set(1);
				break;
			case Stool:
                temp_iter.Set(-1);
                break;
            case Unclassified:
				if (temp_iter.Get()==1) {
					//voxel_type_iter.Set(Stool);			//Modify unclassified into stool holes
                }
				temp_iter.Set(0);
                break;
            case Air:
                temp_iter.Set(-1);
                break;
        }
    }

	WriteITK(voxel_type, "voxel_type_stool_closed.hdr");

	//WriteITK(temp, "temp_tissue.hdr");

	//closingFilter = ClosingFilterType::New();
	//closingFilter->SetKernel(structuringElement);
	//closingFilter->SetInput(temp);
	//closingFilter->Update();
	//temp=closingFilter->GetOutput();
	//temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	//closingFilter.~SmartPointer();

	////WriteITK(temp, "temp_tissue_closed.hdr");

 //   // Updates Information after morphological operations
 //   for (chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(), temp_iter.GoToBegin();
 //       !voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
 //       ++voxel_type_iter, ++chamfer_colon_iter, ++temp_iter) 
 //   {
 //       switch(voxel_type_iter.Get()) {
 //           case Unclassified:
 //               if (temp_iter.Get()==1) {					// unclassified into tissue
	//				voxel_type_iter.Set(Tissue);		
	//				chamfer_colon_iter.Set(0);
	//			} else {
	//				chamfer_colon_iter.Set(1);
	//			}
 //               break;
	//		case Tissue:
	//			chamfer_colon_iter.Set(0);
	//			break;
	//		default:
	//			chamfer_colon_iter.Set(1);
	//			break;
 //       }
 //   }

	// Chamfer colon to tissue
	ByteImageType::Pointer chamfer_colon=AllocateNewByteImage( fullRegion );
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);

	for (chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin();
		!voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();  
		++voxel_type_iter, ++chamfer_colon_iter) 
	{
		switch(voxel_type_iter.Get()) {
			case Tissue:
				chamfer_colon_iter.Set(0);
				break;
			default:
				chamfer_colon_iter.Set(1);
				break;
		}
	}

	//WriteITK(voxel_type, "voxel_type_closed.hdr");
	
	//WriteITK(chamfer_colon, "chamfer_colon_air_pre.hdr");

	// Remove external air, find component with largest pixels on border
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput(chamfer_colon);
	ccFilter->Update();
	IntImageType::Pointer cc = ccFilter->GetOutput();
	ccFilter.~SmartPointer();

	LabelImageToShapeLabelMapFilterType::Pointer converter = LabelImageToShapeLabelMapFilterType::New();
	converter->SetInput( cc );
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();
	converter.~SmartPointer();

	std::vector<LabelObjectType::Pointer> labelVector = labelMap->GetLabelObjects();
	sort(labelVector.begin(), labelVector.end(), compareSizeOnBorder);
	int externalAirLabel = labelVector[0]->GetLabel();
	
	//labelVector.erase( labelVector.begin() );
	//sort(labelVector.begin(), labelVector.end(), compareSize);
	//int colonLabel = labelVector[0]->GetLabel();

	labelVector.clear();
	labelMap.~SmartPointer();

	IteratorTypeIntWithIndex cc_iter( cc, fullRegion );
	chamfer_colon_iter = IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	
	for (chamfer_colon_iter.GoToBegin(), cc_iter.GoToBegin();
        !chamfer_colon_iter.IsAtEnd() && !cc_iter.IsAtEnd();  
        ++chamfer_colon_iter, ++cc_iter) 
	{
		if (cc_iter.Get() == externalAirLabel)
		{
			chamfer_colon_iter.Set(0);
		}
	}

	cc.~SmartPointer();

	//WriteITK(chamfer_colon,"chamfer_colon_no_bkg.hdr");

	// Smooth colon mask using median filter
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput(chamfer_colon);
	medianFilter->SetForegroundValue(1);
	medianFilter->Update();
	chamfer_colon=medianFilter->GetOutput();
	medianFilter.~SmartPointer();
	chamfer_colon_iter=IteratorTypeByteWithIndex(chamfer_colon,fullRegion);

	//WriteITK(chamfer_colon,"chamfer_colon_median.hdr");

	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_colon);
	int chamfer_weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_colon=chamfer_filter->GetOutput();
	chamfer_colon_iter = IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	chamfer_filter.~SmartPointer();

	// Distance map to air
	//WriteITK(chamfer_colon, "distance_map_to_air.hdr");

	chamfer_colon_iter.GoToBegin();
	chamfer_colon_iter.Set(-1);	//first pixel == -1

	for (chamfer_colon_iter.GoToBegin(), input_iter.GoToBegin();
        !chamfer_colon_iter.IsAtEnd() && ! input_iter.IsAtEnd();  
        ++chamfer_colon_iter, ++input_iter) 
    {
		if (chamfer_colon_iter.Get()<20 && chamfer_colon_iter.Get()>-1) {
			chamfer_colon_iter.Set(1);
		} else {
			chamfer_colon_iter.Set(0);
		}
    }

	WriteITK(chamfer_colon, "chamfer_colon_air_mask.hdr");

	// Test fuzzy

	//ImageType::Pointer RunFuzzy(ImageType::Pointer input, ByteImageType::Pointer chamfer_colon, VoxelTypeImage::Pointer voxel_type)

	//RunFuzzy(input_temp,chamfer_colon,voxel_type);

	// Test region growing
	//RunConnectedThresholdGrowing(input, chamfer_colon, voxel_type);


	// Read chamfer colon
	//ReadITK(chamfer_colon, "mr10_092_13p.i0344_100-105_chamfer_colon_air_mask.hdr");
	//chamfer_colon_iter=IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	//WriteITK(chamfer_colon,"chamfer_colon.hdr");

	//-------------------------------------------END SETUP-------------------------------------------------------

	//-------------------------------------------BEGIN QR--------------------------------------------------------
	
    // Get cubic bspline image interpolation
    InterpolationType::Pointer input_interpolator = InterpolationType::New();
    input_interpolator->SetSplineOrder(3);
    input_interpolator->SetInputImage(input);

	// Interpolate gradient magnitude image
	InterpolationType::Pointer gradient_magnitude_interpolator = InterpolationType::New();
	gradient_magnitude_interpolator->SetSplineOrder(3);
	gradient_magnitude_interpolator->SetInputImage(gradient_magnitude);

	// Storage for Smax
	//ImageType::Pointer input_smax = AllocateNewImage(fullRegion);
	ImageType::Pointer input_smax = ComputeNeighborhoodSmax(input_temp, voxel_type, chamfer_colon_iter, startIndex, endIndex);
	WriteITK(input_smax,"smax.hdr");

	// Load smax
	//ReadITK(input_smax, "C:/ImageData/mr10_092_13p.i0344_100-105_smax2.hdr");
	
	IteratorTypeFloat4WithIndex input_smax_iter(input_smax,fullRegion);

	std::cerr << "Running voxel edge classification..." << std::endl;

	//ReadITK(voxel_type, "Modified_6_voxel_edge_class.hdr");
	voxel_type_iter = IteratorTypeVoxelType(voxel_type,fullRegion);

	int count = 0;
	//Runs edge classification with 3 different increments inside the colon only
    for (voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), gradient_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
		!voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !gradient_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++voxel_type_iter, ++input_smax_iter, ++gradient_iter, ++chamfer_colon_iter) {
		
		if (chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get()==Unclassified) {
			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
			CovariantVectorType grad = gradient_iter.Get();	//NOTE: using 2x image gradient, (See Carston paper)

			float temp_threshold=-1;
			VoxelType temp_type = Unclassified;

			if (input_smax_iter.Get() > 0)
			{

				VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);
				
				VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);

				VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,input_interpolator,gradient_magnitude_interpolator,input_smax_iter,voxel_index,grad);

			} else {
				temp_type = TissueAir;
			}

			voxel_type_iter.Set(temp_type);
		}

		if (++count % 100000 == 0)
			std::cout << "Voxel edge count: " << count << std::endl;
    }

	// Deletes the Interpolators
	input_interpolator.~SmartPointer();

	WriteITK(voxel_type,"voxel_edge_class.hdr");

 	std::cerr<<"Done Edge"<<std::endl;

	// Optimize SA transitions using gradient information
	//std::cerr<<"Running voxel edge optimization"<<std::endl;
	
	//OptimizeVoxelEdge(input, voxel_type, gradient);
	
	//WriteITK(voxel_type, "voxel_edge_optimized_noz.hdr");

	// Find thin stool regions [Carston 2.3.6]
	ByteImageType::Pointer chamfer_from_stool_involvement=AllocateNewByteImage(fullRegion);
	IteratorTypeByteWithIndex chamfer_from_stool_involvement_iter(chamfer_from_stool_involvement,fullRegion);

	for(voxel_type_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(); 
		!voxel_type_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(); 
		++voxel_type_iter, ++chamfer_from_stool_involvement_iter) 
    {
		if (voxel_type_iter.Get()==Stool || voxel_type_iter.Get()==TissueStool || voxel_type_iter.Get()==StoolAir) {
			chamfer_from_stool_involvement_iter.Set(0);     //stool
		} else {
			chamfer_from_stool_involvement_iter.Set(1);     //non stool is the object
		}
	}

	//WriteITK(chamfer_from_stool_involvement,"non_stool_mask.hdr");

	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_from_stool_involvement);
	int weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_from_stool_involvement=chamfer_filter->GetOutput();
	chamfer_from_stool_involvement_iter = IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	chamfer_filter.~SmartPointer();

	//WriteITK(chamfer_from_stool_involvement, "non_stool_distance.hdr");

	// Find thin stool voxels < 8 and no neighbors > 8 (which were not already classified as stool)
	int chamfercutoff = 8;
	int thinStoolCount = 0;

	for(voxel_type_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
		!voxel_type_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(); 
		++voxel_type_iter, ++chamfer_from_stool_involvement_iter) 
    {
		if (voxel_type_iter.Get() != Stool && voxel_type_iter.Get() != TissueStool) // do not turn these into thin stool
		{
			if (chamfer_from_stool_involvement_iter.Get() < chamfercutoff && chamfer_from_stool_involvement_iter.Get() > 0)
			{
				ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();

				bool thinNeighbors = true;

				for(int i=-1;i<=1; i++) {
					if (thinNeighbors && index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1 && thinNeighbors; j++) {
							if (thinNeighbors && index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1; k++) {
									if (thinNeighbors && index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										
										ImageType::IndexType neighbor_index={index[0]+i,index[1]+j,index[2]+k};
										
										if (chamfer_from_stool_involvement->GetPixel(neighbor_index) > chamfercutoff) {
											thinNeighbors = false;
										}
									
									}
								}
							}
						}
					}
				}

				if (thinNeighbors)
				{
					voxel_type_iter.Set(ThinStool);
					thinStoolCount++;
				}
			}
		}
	}

	std::cout << "Number of thin stool voxels: " << thinStoolCount << std::endl;

	WriteITK(voxel_type,"voxel_type_thinstool.hdr");

	// Apply voting scheme to voxel edge
	//VoteVoxels(voxel_type, chamfer_colon);
	
	//voxel_type_iter=IteratorTypeVoxelType(voxel_type,fullRegion);
	
	//WriteITK(voxel_type, "voxel_edge_voting.hdr");

	// Find thin stool distance to air
	ByteImageType::Pointer chamfer_air=AllocateNewByteImage( fullRegion );
    IteratorTypeByteWithIndex chamfer_air_iter(chamfer_air,fullRegion);

    for(voxel_type_iter.GoToBegin(), chamfer_air_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd();
        ++voxel_type_iter, ++chamfer_air_iter)
    {
		if (voxel_type_iter.Get()==Air) {
            chamfer_air_iter.Set(1);					
        } else {
            chamfer_air_iter.Set(0);	
		}
    }

	//compute chamfer for air
	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_air);
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_air=chamfer_filter->GetOutput();
	chamfer_air_iter= IteratorTypeByteWithIndex(chamfer_air,fullRegion);
	chamfer_filter.~SmartPointer();
	//WriteITK(chamfer_air,"chamfer_to_air.hdr");
	
	// Convert stool and air chamfers to float and normalize
	ImageType::Pointer ds = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex ds_iter(ds,fullRegion);

	ImageType::Pointer da = AllocateNewImage(fullRegion);
	IteratorTypeFloat4WithIndex da_iter(da,fullRegion);

	for(chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), da_iter.GoToBegin(), ds_iter.GoToBegin();
		!chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(), !da_iter.IsAtEnd(), !ds_iter.IsAtEnd();
        ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter, ++da_iter, ++ds_iter)
    {
		if ( chamfer_from_stool_involvement_iter.Get() >= 8)
		{
			ds_iter.Set(1);
		} else {
			ds_iter.Set( (float) chamfer_from_stool_involvement_iter.Get() / 8 );
		}

		if ( chamfer_air_iter.Get() >= 8)
		{
			da_iter.Set(1);
		} else {
			da_iter.Set( (float) chamfer_air_iter.Get() / 8 );
		}
	}

	da_iter = IteratorTypeFloat4WithIndex(da,fullRegion);
	ds_iter = IteratorTypeFloat4WithIndex(ds,fullRegion);

	//WriteITK(da,"normalized_air.hdr");
	//WriteITK(ds,"normalized_stool.hdr");

	//-------------------------------------------------COMPUTE PARTIAL VOLUMES---------------------------------------------------------------
	
	// Declares a structure to store the air/tissue/stool partial volumes for a voxel in the image
	ImageVectorType::Pointer partialVector = ImageVectorType::New();
    partialVector->SetRegions(fullRegion);
    partialVector->Allocate();
	partialVector->SetSpacing(input->GetSpacing());
    IteratorImageVectorType partialVector_iter(partialVector,fullRegion);

	int total_vertices=0;

	std::ofstream file;
	file.open("pv.txt");
	file << "Index\tIntensity\tGradient\tClass\tPa\tPt\tPs\n";

	for(partialVector_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), gradient_magnitude_iter.GoToBegin(), da_iter.GoToBegin(), ds_iter.GoToBegin();
		!partialVector_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !gradient_magnitude_iter.IsAtEnd() && !da_iter.IsAtEnd() && !ds_iter.IsAtEnd() ; 
        ++partialVector_iter, ++input_iter, ++chamfer_colon_iter, ++voxel_type_iter, ++input_smax_iter, ++gradient_magnitude_iter, ++da_iter, ++ds_iter) 
	{
		if (chamfer_colon_iter.Get() == 1)
		{
			total_vertices++;

			ImageType::IndexType idx = input_iter.GetIndex();
			file << "(" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t";
			file << input_iter.Get() << "\t" << gradient_magnitude_iter.Get() << "\t";
			
			float value[3];	// air/tissue/stool
			
			switch (voxel_type_iter.Get()) {
				case Air:
					value[0]=1;					
					value[1]=0;
					value[2]=0;

					file << "Air";
					break;
				case Tissue:
					value[0]=0;
					value[1]=1;					
					value[2]=0;

					file << "Tissue";
					break;
				case Stool:
					value[0]=0;
					value[1]=0;
					value[2]=1;		

					file << "Stool";
					break;
				case TissueAir:
					value[1]=1+(input_iter.Get()/1000);
					
					if (value[1] <= 0) { value[1] = 0; }
					if (value[1] >= 1) { value[1] = 1; }

					value[0]=1-value[1];
					value[2]=0;

					file << "TissueAir";
					break;
				case TissueStool:
					value[0]=0;
					value[2]=input_iter.Get()/input_smax_iter.Get();
				
					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[1]=1-value[2];

					file << "TissueStool";
					break; 
				case StoolAir:
					value[2]=(input_iter.Get()+1000)/(input_smax_iter.Get()+1000);
					
					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[0]=1-value[2];
					value[1]=0;

					file << "StoolAir";
					break;
				case ThinStool:
					value[0]=0.5*(1-vnl_erf((da_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2)))); //check eq
					
					if (value[0] <= 0) { value[0] = 0; }
					if (value[0] >= 1) { value[0] = 1; }

					value[2]=0.5*(1+vnl_erf((ds_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));

					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[1]=1-value[0]-value[2];

					if (value[1] <= 0) { value[1] = 0; }
					if (value[1] >= 1) { value[1] = 1; }

					file << "ThinStool";
					break;
			}   

			file << "\t" << value[0] << "\t" << value[1] << "\t" << value[2] << "\n";

			CovariantVectorType data(value);	//sets the partial to be stored in the actual structure
			partialVector_iter.Set(data);		//stores the partial in the structure
		}
    }

	std::cout << "Partials initialized" << std::endl;

	partialVector_iter = IteratorImageVectorType(partialVector,fullRegion);

	file.close();

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

	//WriteITK(pt_mask,"pt_mask.hdr");

	// Run cc on pt mask
	ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput(pt_mask);
	ccFilter->Update();
	cc=ccFilter->GetOutput();

	//WriteITK(cc, "pt_mask_cc.hdr");

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

	// Smooth partial vector
	SmoothPartialVector(partialVector,chamfer_colon,startIndex,endIndex);

	WritePartialImages(partialVector, chamfer_colon, "QR_Smooth");

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

	WriteITK(output,"outputQR.hdr");

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

	WriteITK(input_temp,"output.hdr");
	*/

	//-------------------------------------------END QR-----------------------------------------------------

	//-------------------------------------------BEGIN EM--------------------------------------------------------

    int counter=0;							//used to total certain values (description given when set)
    double mean[3]={0, 0, 0};				//stores the mean of the air/tissue/stool classes
    double sum[3]={0, 0, 0};				//stores the total # of partials of the three classes
    double variance[3]={0, 0, 0};			//stores the variance of the air/tissue/stool classes
	float weight[3]={0,0,0};				//stores the weights of the air/tissue/stool classes

	// Computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
	for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
		!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
		++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	{
		if (chamfer_colon_iter.Get()==1) {		
			CovariantVectorType partial = partialVector_iter.Get();
			for (int i=0; i<3; i++)
			{
				sum[i] += partial[i];
				mean[i] += partial[i]*input_iter.Get();
			}
		}
	}
	for (int i=0;i<3;i++) { mean[i]=mean[i]/sum[i]; } 

	// Compute variance and weights
	for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
		!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
		++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	{
		if (chamfer_colon_iter.Get()==1) {		
			CovariantVectorType partial = partialVector_iter.Get();
			for (int i=0; i<3; i++)
			{
				variance[i] += partial[i]*vnl_math_sqr(input_iter.Get()-mean[i]);
			}
		}
	}

	double sum_all = sum[0]+sum[1]+sum[2];

	for (int i=0;i<3;i++) 
	{ 
		variance[i]=variance[i]/sum[i];
		weight[i]=sum[i]/sum_all;
	} 

	//outputs the initial mean and variance of each class
	std::ofstream em;
	em.open("em.txt");
	em<<"EM\tMean0\tMean1\tMean2\tVar0\tVar1\tVar2\tWeight0\tWeight1\tWeight2\n";

	std::ostringstream ss;
	ss<<"0\t"<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<variance[0]<<"\t"<<variance[1]<<"\t"<<variance[2]<<"\t"<<weight[0]<<"\t"<<weight[1]<<"\t"<<weight[2]<<"\n";
	em<<ss.str();

	std::cerr<<std::endl;
	std::cerr<<"EM0"<<std::endl;
	std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;
	
	// Computes iterations of the Maximization Algorithm
    for (int emNum=0;emNum<20;emNum++) // used 20 in original test case
	{
		//gives the temporary storage of the variables corresponding to sum, variance, and mean for each class on the i+1 iteration
		double sum_temp[3]={0,0,0};
        double variance_temp[3]={0,0,0};
		double mean_temp[3]={0,0,0};

		for (input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
			!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
			++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
		{
			if (chamfer_colon_iter.Get()==1) {		
				CovariantVectorType partial = partialVector_iter.Get();		//retrieves the partial informations
				float Z[3]={partial[0],partial[1],partial[2]};

				// Only update uncertain partials
				if ( !(partial[0] == 1 || partial[1] == 1 || partial[2] == 1) )
				{
					vnl_vector<float> Z_update=expectation(input_iter.Get(),mean, variance, weight, GetNeighbor(partialVector,input_iter.GetIndex()), Z);	//updates the partial values
					partial[0]=Z_update[0];										
					partial[1]=Z_update[1];
					partial[2]=Z_update[2];
					partialVector_iter.Set(partial);
				}

				//updates the new mean total partial sum for each class accordingly
				for (int i=0;i<3;i++) 
				{
					mean_temp[i]+=partial[i]*input_iter.Get();	
					sum_temp[i]+=partial[i];
				}

			}
        }

		partialVector_iter = IteratorImageVectorType(partialVector,fullRegion);

		for (int i=0;i<3;i++) { mean_temp[i]=mean_temp[i]/sum_temp[i]; } 

		// Compute variance and weights
		for(input_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin();
			!input_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();
			++input_iter, ++partialVector_iter, ++chamfer_colon_iter) 
		{
			if (chamfer_colon_iter.Get()==1) {		
				CovariantVectorType partial = partialVector_iter.Get();
				for (int i=0; i<3; i++)
				{
					variance_temp[i] += partial[i]*vnl_math_sqr(input_iter.Get()-mean_temp[i]);
				}
			}
		}

		double sum_all_temp = sum_temp[0]+sum_temp[1]+sum_temp[2];

		for (int i=0;i<3;i++) 
		{ 
			mean[i]=mean_temp[i];
			variance[i]=variance_temp[i]/sum_temp[i];
			weight[i]=sum_temp[i]/sum_all_temp;
		}

		std::cerr<<std::endl;
		std::cerr<<"EM"<<emNum+1<<std::endl;
		std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
		std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
		std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;

		std::stringstream ss;
		ss<<emNum+1;
		ss<<"\t"<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<variance[0]<<"\t"<<variance[1]<<"\t"<<variance[2]<<"\t"<<weight[0]<<"\t"<<weight[1]<<"\t"<<weight[2]<<"\n";
		em<<ss.str();

		std::stringstream ss2;
		ss2<<"EM"<<emNum+1;

		WritePartialImages(partialVector,chamfer_colon,ss2.str());
    }

	em.close();

	////set modification images type.
	////Updates the existing voxel type classifications with knowledge found in the EM
	//EMClassification(input_iter, voxel_type_iter, partialVector_iter, chamfer_colon_iter, temp_iter);
	//WriteITK(temp,"Global_EM_classification.hdr");

	//for(input_iter.GoToBegin(), em_output_image_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
	//	!input_iter.IsAtEnd() && !em_output_image_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
	//	++input_iter, ++em_output_image_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//{
	//	CovariantVectorType data = partialVector_iter.Get();		//retrieves the partial informations
	//	em_output_image_iter.Set(data[1]);		//sets the voxel intensity for the i-th iteration output image
	//}

	//WriteITK(em_output_image, "tissue_partial.hdr");

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
	//WriteITK(partialVolume,"tissue_partial_one.hdr");

	////Write out unmodified image
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_input.hdr");

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
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool.hdr");

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
	//	WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_forward.hdr");
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
	//	WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool_back.hdr");
	//}

	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_stool.hdr");



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
	//WriteITK(chamfer_air,"chamfer_from_air.hdr");
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
	//WriteITK(partialVolume,"tissue_partial_two.hdr");

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
	//WriteITK(output,"output_one.hdr");
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
	//WriteITK(chamfer_from_stool_involvement,"chamfer_from_tissue_input.hdr");
	//chamfer_filter = ChamferDistanceFilterType::New();
	//chamfer_filter->SetInput(chamfer_from_stool_involvement);
	//chamfer_filter->SetWeights(weights, weights+3);
	//chamfer_filter->SetDistanceFromObject(false);
	//chamfer_filter->Update();

	//chamfer_from_stool_involvement = chamfer_filter->GetOutput();
	//chamfer_from_stool_involvement_iter=IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	//chamfer_filter.~SmartPointer();

	////Write out chamfer thickness of stool
	//WriteITK(chamfer_from_stool_involvement,"chamfer_tissue_stool.hdr");

	////-------------------------------------------BEGIN PARTIAL VOLUME COMP----------------------------------
	

	//Write out the final output	
	std::cerr<<"Ended"<<std::endl;
	//WriteITK(output, "output.hdr");
	system("pause");
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
	//float stool_tissue_Smax;
	//float stool_air_Smax;
	
	//Tissue Stool

	//stool_tissue_Smax=ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
	//stool_air_Smax=Stool_Air_ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
	
	float smax = input_smax.Get();
	//float smax=0;

	float distance=0;
	VoxelType voxel_type=Unclassified;
	
	//float distanceTS=AverageTissueStoolDist(stool_tissue_Smax, temp_intensity,temp_gradient_magnitude);
	//float distanceSA=AverageStoolAirDist(stool_air_Smax, temp_intensity,temp_gradient_magnitude); // bias stool air distance to preserve tissue

	float distanceTS=AverageTissueStoolDist(smax, temp_intensity,temp_gradient_magnitude);
	float distanceSA=AverageStoolAirDist(smax, temp_intensity,temp_gradient_magnitude); // bias stool air distance to preserve tissue
	float distanceTA=AverageTissueAirDist(temp_intensity,temp_gradient_magnitude);

	if (distanceSA<=distanceTS && distanceSA<=distanceTA){
		distance=distanceSA;
		voxel_type=StoolAir;
	} else if (distanceTS<=distanceTA && distanceTS<=distanceSA ) {
		distance=distanceTS;
		voxel_type=TissueStool;
	} else {
		distance=distanceTA;
		voxel_type=TissueAir;
	}

    if (min_distance>=distance || min_distance==-1 || *previous==Unclassified) {
		*threshold=distance;
        *previous=voxel_type;
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
	
	/*if (Modified) {
		coefficients[0]= -6.0/1000;
		coefficients[1]= -6.0;
		coefficients[2]=	0;
	}*/

	if (Modified) {
		coefficients[0]= -4.0/1000;
		coefficients[1]= -4.0;
		coefficients[2]=	0;
	}

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
VoxelType SingleMaterialClassification(ImageType::PixelType input_pixel, ImageType::PixelType input_gradient_pixel) {
	if (Modified) {

		/*
		if ((input_pixel >=800 && input_gradient_pixel<=0.8*input_pixel) || input_pixel >= 1000) {
			return Stool;
		} else if (input_pixel<=-800) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
			return Tissue;
		} else {
			return Unclassified;
		}
		*/

		/*

		if (input_pixel >=650) {
			return Stool;
		} else if (input_pixel<=-600) {
			return Air;
		} else if (input_pixel<=300  && input_pixel>=-300 && input_gradient_pixel<=700) {
			return Tissue;
		} else {
			return Unclassified;
		}

		*/

		if (input_pixel >=650) {
			return Stool;
		} else if (input_pixel<=-600) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=400) {
			return Tissue;
		} else {
			return Unclassified;
		}

	} else {
		if ((input_pixel >=180 && input_gradient_pixel<=0.8*input_pixel)) {
			return Stool;
		} else if (input_pixel<=-800 && input_gradient_pixel<=250) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
			return Tissue;
		} else {
			return Unclassified;
		}
	}
    return Stool;
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

void OptimizeVoxelEdge(ImageType::Pointer input, VoxelTypeImage::Pointer voxelEdge, ImageVectorType::Pointer gradient ) {

	// Get info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	
	// Create iterators
	IteratorTypeFloat4WithIndex inputIt(input,region);
	IteratorImageVectorType gradientIt(gradient,region);
	IteratorTypeVoxelType voxelEdgeIt(voxelEdge,region);

	// Pattern matching
	char *regex = "^2*[5-7]+1+$";

	int numOfVoxels = 10; // to move in direction of gradient

	// Setup text file
	std::ofstream myfile;
	myfile.open ("OptimizeVoxelEdgeOutput.txt");

	myfile << "Number of voxels to move in direction of gradient: " << numOfVoxels << "\n";

	// Setup boundaries
	ImageType::IndexType endIndex = region.GetIndex();
	ImageType::IndexType startIndex = region.GetIndex();	
	endIndex[0]+=(region.GetSize()[0]-1);
	endIndex[1]+=(region.GetSize()[1]-1);
	endIndex[2]+=(region.GetSize()[2]-1);

	// Copy voxel edge into new image
	VoxelTypeImage::Pointer voxelEdgeUpdate = VoxelTypeImage::New();
	voxelEdgeUpdate->SetRegions( region );
	voxelEdgeUpdate->Allocate();
	voxelEdgeUpdate->SetSpacing( spacing );

	IteratorTypeVoxelType voxelEdgeUpdateIt ( voxelEdgeUpdate, region );

	for(voxelEdgeIt.GoToBegin(), voxelEdgeUpdateIt.GoToBegin();
		!voxelEdgeIt.IsAtEnd() && !voxelEdgeUpdateIt.IsAtEnd(); 
		++voxelEdgeIt, ++voxelEdgeUpdateIt) 
	{
		voxelEdgeUpdateIt.Set(voxelEdgeIt.Get());
	}	
	
	// Make binary air mask
	std::cout << "Making air mask" << std::endl;
	
	ImageType::Pointer airMask = ImageType::New();
	airMask->SetRegions( region );
	airMask->Allocate();
	airMask->SetSpacing( spacing );
	
	IteratorTypeFloat4WithIndex airMaskIt( airMask, region );

	for (	voxelEdgeIt.GoToBegin(), airMaskIt.GoToBegin();
			!voxelEdgeIt.IsAtEnd() && !airMaskIt.IsAtEnd();
			++voxelEdgeIt, ++airMaskIt)
	{
		if (voxelEdgeIt.Get() == Air)
		{
			airMaskIt.Set(1);
		} else {
			airMaskIt.Set(0);
		}
	}

	WriteITK(airMask, "airMask.hdr");

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

	WriteITK( airMask, "airMaskEdges.hdr");

	std::cout << "Finding voxels by gradient" << std::endl;

	int countMatches = 0;
	int countEdges = 0;

	for (	airMaskIt.GoToBegin(), gradientIt.GoToBegin(), inputIt.GoToBegin(); 
			!airMaskIt.IsAtEnd() && !gradientIt.IsAtEnd() && !inputIt.IsAtEnd(); 
			++airMaskIt, ++gradientIt, ++inputIt)
	{
		
		if (airMaskIt.Get() == 2) { // at edge
			countEdges++;

			ImageType::IndexType idx = inputIt.GetIndex();
			CovariantVectorType grad = gradientIt.Get();

			// Output gradient from filter to text file
			myfile << "Index: (" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t" << inputIt.Get() << "\t";
			myfile << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\t";

			//std::cout << "Index: (" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t" << inputIt.Get() << "\t";
			//std::cout << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\t\n";

			std::vector<ImageType::IndexType> indexVector; //store indices of voxels

			// Make sure there is a non-zero gradient
			if ( !(grad[0]==0 && grad[1]== 0 && grad[2]==0) )
			{
				FindVoxelsByGradient2(voxelEdge, idx, startIndex, endIndex, grad, numOfVoxels, indexVector);

				std::stringstream ss;
				
				for (int i=0; i < indexVector.size() ; i++)
				{
					ss << VoxelTypeToNum( voxelEdge->GetPixel( indexVector[i] ) );
				}

				std::string vs = ss.str();

				myfile << vs;

				// Match pattern of voxels
				bool isMatch = MatchVoxels( vs );
				
				if( isMatch )
				{
					myfile << "\tMATCH";
					countMatches++;

					// Convert TA/TS to SA transitions
					for (int i=0; i < vs.length(); i++)
					{
						if (vs.compare(i,1,"6") == 0 || vs.compare(i,1,"7") == 0)
						{
							voxelEdgeUpdate->SetPixel( indexVector[i], StoolAir );
						}
					}
				}
				myfile << "\n";
			}
		}
	}
	
	myfile << "\n";
	myfile << "Number of air edge voxels patterns: " << countEdges << "\n";
	myfile << "Number of matched patterns: " << countMatches << " (" << 100*countMatches/countEdges << "%)\n";

	myfile.close();

	// Copy updated voxel info into old voxel image
	for(voxelEdgeIt.GoToBegin(), voxelEdgeUpdateIt.GoToBegin();
		!voxelEdgeIt.IsAtEnd() && !voxelEdgeUpdateIt.IsAtEnd(); 
		++voxelEdgeIt, ++voxelEdgeUpdateIt) 
	{
		voxelEdgeIt.Set(voxelEdgeUpdateIt.Get());
	}
}

void FindVoxelsByGradient(VoxelTypeImage::Pointer voxelEdge, ImageType::IndexType &index, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex, CovariantVectorType &grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector) {
	//Calculates indices of voxels in direction of gradient
	ImageType::IndexType idx = index;
	ContinuousIndexType offsetIndex;

	// Space increment
	double ds = 0.1;
	VoxelTypeImage::PixelType voxel = voxelEdge->GetPixel( index );
	
	int i=0;
	while ( (i < numOfVoxels) && (voxel != Stool) && (voxel != Tissue) ) // stop finding once you hit stool or tissue
	{
		// Reset new offset index
		offsetIndex[0] = idx[0];
		offsetIndex[1] = idx[1];
		offsetIndex[2] = idx[2];
		
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

		} while (	(idx[0] == round(offsetIndex[0])) &&
					(idx[1] == round(offsetIndex[1])) &&
					(idx[2] == round(offsetIndex[2]))	);	// Until a new discrete voxel has been found
		
		// Set new index and store it
		for (int k=0; k < 3; k++) { idx[k] = round(offsetIndex[k]); }

		voxel = voxelEdge->GetPixel( idx );
		
		indexVector.push_back(idx);

		i++;
	}	
}

void FindVoxelsByGradient2(VoxelTypeImage::Pointer voxelEdge, ImageType::IndexType &index, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex, CovariantVectorType &grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector) {
	ImageType::IndexType idx = index;
	VoxelType voxel = voxelEdge->GetPixel( idx );

	int i=0;
	while ( (i < numOfVoxels) && (voxel != Stool) && (voxel != Tissue) ) // stop finding once you hit stool or tissue
	{
		idx[0] += round( grad[0] );
		idx[1] += round( grad[1] );
		//idx[2] += round( grad[2] );

		// Check bounds
		//for (int j=0; j < 3; j++)
		for (int j=0; j < 2; j++)
		{
			if (idx[j] < startIndex[j]) {	idx[j] = startIndex[j]; }
			if (idx[j] > endIndex[j])   {	idx[j] = endIndex[j];	}
		}

		voxel = voxelEdge->GetPixel(idx);
		indexVector.push_back(idx);

		i++;
	}	
}

double round(float d)
{
  return floor(d + 0.5);
}

void WriteITK(ImageType::Pointer image, std::string name) {
	// Use ITK Writer
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if (Modified)
		ss << "Modified_";
	if ( !note.empty() )
		ss << note << "_";
	ss << writeCount++ << "_" << name;
	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
	writer.~SmartPointer();
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
	ss << writeCount++ << "_" << name;
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
	ss << writeCount++ << "_" << name;
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
	ss << writeCount++ << "_" << name;
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

ImageType::Pointer ReadDicom( std::string path )
{	
	// Create reader
	itk::ImageSeriesReader<ImageType>::Pointer reader = itk::ImageSeriesReader<ImageType>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	// Truncate data
	//names.erase( names.begin()+10, names.end() );

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return 0;
	}

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI);
    orienter->SetInput(reader->GetOutput());
    orienter->Update();

    return orienter->GetOutput();

}

void SmoothPartialVector(ImageVectorType::Pointer pv, ByteImageType::Pointer chamfer_colon, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex)
{
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,pv->GetLargestPossibleRegion());

	// Find edge of chamfer colon
	//ByteImageType::Pointer chamfer_colon_edge = FindBinaryEdge(chamfer_colon,startIndex,endIndex);
	//IteratorTypeByteWithIndex chamfer_colon_edge_iter(chamfer_colon_edge,pv->GetLargestPossibleRegion());
	//WriteITK(chamfer_colon_edge,"chamfer_colon_edge.hdr");

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

		//WriteITK(part,"part.hdr");

		// Apply smoothing filter to one dimension
		DiscreteGaussianFilterType::Pointer dgFilter = DiscreteGaussianFilterType::New();
		dgFilter->SetInput(part);

		std::cout << "Smoothing partial by " << pv->GetSpacing()[0] << "mm" << std::endl;

		dgFilter->SetVariance(pv->GetSpacing()[0]*pv->GetSpacing()[0]);
		dgFilter->Update();
		part=dgFilter->GetOutput();
		part_iter=IteratorTypeFloat4WithIndex(part,part->GetLargestPossibleRegion());

		//WriteITK(part,"part_smooth.hdr");

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

ImageType::Pointer ComputeHessianResponse(ImageType::Pointer input)
{

	// Resample input to make it anisotropic
	//input = ResampleImage(input);
	//WriteITK(input, "input_resample.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	// Threshold to compute only in tagged regions
	IteratorTypeFloat4WithIndex input_iter(input,region);

	ByteImageType::Pointer tagged = AllocateNewByteImage(region);
	IteratorTypeByteWithIndex tagged_iter(tagged,region);

	for (input_iter.GoToBegin(), tagged_iter.GoToBegin();
		 !input_iter.IsAtEnd() && !tagged_iter.IsAtEnd();
		 ++input_iter, ++tagged_iter)
	{
		if (input_iter.Get() > 200) 
		{
			tagged_iter.Set(1);
		} else {
			tagged_iter.Set(0);
		}
	}

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[7] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0], 6*spacing[0], 7*spacing[0]};

	for (int k=0; k<5; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);
		std::vector< ImageType::Pointer > FImageVector(7); // A, B, C1, C2, Cup, Rut, Ralpha

		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);
		std::vector< itk::ImageRegionIterator<ImageType> > FIterVector(7);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		for (int i=0; i<7; i++)
		{
			FImageVector[i] = AllocateNewImage(region);
			FIterVector[i] = itk::ImageRegionIterator<ImageType>( FImageVector[i], region );
			FIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		tagged_iter.GoToBegin();
		Himage_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			//if (tagged_iter.Get() == 1) { // compute only in tagged region
				
				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				//std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				if (count % 50000 == 0)
				{
					std::cout << "what" << std::endl;
				}	

				for (int i=0; i<3; i++)
				{
					LambdaIterVector[i].Set( lambda[i] );
				}

				if ( lambda[2] < 0 )
				{
					// A, B, C1, C2, Cup, Rut
					FIterVector[0].Set( fA(lambda, ALPHA) );
					FIterVector[1].Set( fB(lambda, BETA, GAMMA) );
					FIterVector[2].Set( fC(lambda[0], lambda[1], ETA) );
					FIterVector[3].Set( fC(lambda[1], lambda[2], ETA) );
					FIterVector[4].Set( fCup(lambda, ETA) );
					FIterVector[5].Set( fRut(lambda, ALPHA, BETA, GAMMA) );


					FIterVector[6].Set( abs(lambda[0])/sqrt(abs(lambda[1]*lambda[2]) ));

					Himage_iter.Set( vnl_math_max( fRut(lambda, ALPHA, BETA, GAMMA), fCup(lambda, ETA) ) );	
				} else {
					Himage_iter.Set(0);
//
					for (int i=0; i<6; i++) { FIterVector[i].Set(0); }
				}
				

							
			//} else {
			//	Himage_iter.Set(0);
			//	for (int i=0; i<3; i++) { LambdaIterVector[i].Set(0); }
			//	for (int i=0; i<7; i++) { FIterVector[i].Set(0); }
			//}

			++hessian_iter;
			++tagged_iter;
			++Himage_iter;
			count++;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
			for (int i=0; i<7; i++) { ++FIterVector[i]; }
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "H.hdr";
		WriteITK(Himage, ss.str());

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			WriteITK(LambdaImageVector[i], ss.str());
		}
	
		/*
		for (int i=0; i<3; i++)
		{
			RescaleIntensityFilterType::Pointer rescaleFilter = RescaleIntensityFilterType::New();
			rescaleFilter->SetInput( LambdaImageVector[i] );
			rescaleFilter->SetOutputMinimum(0);
			rescaleFilter->SetOutputMaximum(1);
			rescaleFilter->Update();
			LambdaImageVector[i]=rescaleFilter->GetOutput();
			LambdaIterVector[i]=itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);

			for (LambdaIterVector[i].GoToBegin(), tagged_iter.GoToBegin();
				 !LambdaIterVector[i].IsAtEnd() && !tagged_iter.IsAtEnd();
				 ++LambdaIterVector[i], ++tagged_iter)
			{
				if (tagged_iter.Get() != 1)
				{
					LambdaIterVector[i].Set(0);
				}
			}

			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_rescale_" << i << ".hdr";
			WriteITK(LambdaImageVector[i], ss.str());
		}
		*/

		for (int i=0; i<7; i++)
		{
			//if (i==0 || i==6) {
				std::stringstream ss;
				ss << "sigma_" << sigma[k] << "_F_";
				
				// A, B, C1, C2, Cup, Rut
				switch (i)
				{
					case 0: ss<<"a_alpha_"<<ALPHA; break;
					case 1: ss<<"b"; break;
					case 2: ss<<"c1"; break;
					case 3: ss<<"c2"; break;
					case 4: ss<<"cup"; break;
					case 5: ss<<"rut"; break;
					case 6: ss<<"Ralpha"; break;
				}

				ss << ".hdr";
				
				WriteITK(FImageVector[i], ss.str());
			//}
		}
	}

	

	return Himage;
}

void MeasureObjectness(ImageType::Pointer input)
{
	// test over sigma
	// test scale by largest eigenvalue

	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	float sigma=0;

	// Compute smoothed Hessian
	HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
	hessianFilter->SetInput(input);
	hessianFilter->SetNormalizeAcrossScale(true);
	hessianFilter->SetSigma(sigma);
	hessianFilter->Update();

	// Get object measures at each sigma
	for (int i=0; i<2; i++)
	{
		sigma = spacing[0]*(i+1);

		for (int m=0; m<3; m++) // object type
		{
			ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
			objectnessFilter->SetInput(hessianFilter->GetOutput());
			objectnessFilter->SetBrightObject(true);
			objectnessFilter->SetObjectDimension(m);
			objectnessFilter->SetScaleObjectnessMeasure(true);

			objectnessFilter->Update();

			std::stringstream ss;
			ss << "sigma_" << sigma << "_";

			switch (m)
			{
				case 0:
					ss << "blob";
					break;
				case 1:
					ss << "vessel";
					break;
				case 2:
					ss << "plate";
					break;
			}

			ss << "_scaled";

			ss << ".hdr";

			WriteITK(objectnessFilter->GetOutput(),ss.str());

		}	
	}
}

void ComputeVesselness(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	for (int k=0; k<2; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);
		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		Himage_iter.GoToBegin();
		
		float max=0;
		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			// Get eigenvalues
			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);
			//std::sort(lambda.Begin(),lambda.End(),OrderByValue);

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			Himage_iter.Set( vesselness(lambda[0], lambda[1], lambda[2], 0.25, 1, 1) );	

			count++;

			++hessian_iter;
			++Himage_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_Vesselness.hdr";
		WriteITK(Himage,ss.str());
	}
}

ImageType::Pointer ComputeHessianResponse2(ImageType::Pointer input)
{
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType::IndexType origin = fullRegion.GetIndex();
	ImageType::SizeType size = fullRegion.GetSize();
	ImageType::SpacingType spacing = input->GetSpacing();

	ImageType::Pointer output = AllocateNewImage(fullRegion);
	
	IteratorTypeFloat4WithIndex input_iter(input,fullRegion);
	IteratorTypeFloat4WithIndex output_iter(output,fullRegion);	
	
	NeighborhoodIteratorType::RadiusType radius;
	radius[0] = radius[1] = radius[2] = 1;

	DerivativeImageFilterType::Pointer derivative_filter = DerivativeImageFilterType::New();
	derivative_filter->SetInput(input);
	derivative_filter->SetOrder(1);
	derivative_filter->SetDirection(0);
	derivative_filter->SetUseImageSpacing(true);
	derivative_filter->Update();

	ImageType::Pointer x_partial = derivative_filter->GetOutput();
	IteratorTypeFloat4WithIndex x_iter(x_partial,fullRegion);

	derivative_filter.~SmartPointer();
	derivative_filter = DerivativeImageFilterType::New();
	derivative_filter->SetInput(input);
	derivative_filter->SetOrder(1);
	derivative_filter->SetDirection(1);
	derivative_filter->SetUseImageSpacing(true);
	derivative_filter->Update();

	ImageType::Pointer y_partial = derivative_filter->GetOutput();
	IteratorTypeFloat4WithIndex y_iter(y_partial,fullRegion);

	derivative_filter.~SmartPointer();
	derivative_filter = DerivativeImageFilterType::New();
	derivative_filter->SetInput(input);
	derivative_filter->SetOrder(1);
	derivative_filter->SetDirection(2);
	derivative_filter->SetUseImageSpacing(true);
	derivative_filter->Update();

	ImageType::Pointer z_partial = derivative_filter->GetOutput();
	IteratorTypeFloat4WithIndex z_iter(z_partial,fullRegion);

	derivative_filter.~SmartPointer();

	WriteITK(x_partial,"x_partial.hdr");
	WriteITK(y_partial,"y_partial.hdr");
	WriteITK(z_partial,"z_partial.hdr");

	float weights[3];
	weights[0] = 1/spacing[0];
	weights[1] = 1/spacing[1];
	weights[2] = 1/spacing[2];

	ImageType::Pointer lamda1 = AllocateNewImage(fullRegion);
	ImageType::Pointer lamda2 = AllocateNewImage(fullRegion);
	ImageType::Pointer lamda3 = AllocateNewImage(fullRegion);

	IteratorTypeFloat4WithIndex lamda1_iter(lamda1,fullRegion);
	IteratorTypeFloat4WithIndex lamda2_iter(lamda2,fullRegion);
	IteratorTypeFloat4WithIndex lamda3_iter(lamda3,fullRegion);

	int count=0;
	float max=0;

	input_iter.GoToBegin();
	output_iter.GoToBegin();
	x_iter.GoToBegin();
	y_iter.GoToBegin();
	z_iter.GoToBegin();
	lamda1_iter.GoToBegin(); lamda2_iter.GoToBegin(); lamda3_iter.GoToBegin();

	while (!input_iter.IsAtEnd())
	{
		vnl_matrix_fixed<float ,3,3> Next;
		vnl_matrix_fixed<float ,3,3> Previous;
		for (int i=0; i<3; i++) {
			ImageType::IndexType next_index = input_iter.GetIndex();
			next_index[i]=next_index[i]+1;
			if (next_index[i]>size[i]-1+origin[i]) {
				next_index[i]=size[i]-1+origin[i];
			}

			ImageType::IndexType previous_index = input_iter.GetIndex();
			previous_index[i]=previous_index[i]-1;
			if (previous_index[i]<origin[i]){
				previous_index[i]=origin[i];
			}

			for (int j=0; j<3; j++) {
				float nextvalue=0;
				float previousvalue=0;
				switch (j) {
					case 0:
						nextvalue = x_partial->GetPixel(next_index);
						previousvalue = x_partial->GetPixel(previous_index);
						break;
					case 1:
						nextvalue = y_partial->GetPixel(next_index);
						previousvalue = y_partial->GetPixel(previous_index);
						break;
					case 2:
						nextvalue = z_partial->GetPixel(next_index);
						previousvalue = z_partial->GetPixel(previous_index);
						break;
				}
				Next[i][j]=nextvalue;
				Previous[i][j]=previousvalue;
			}
		}
		
		vnl_matrix<float> temp_return = EvaluateAtNeighborhood(Next, Previous, weights);
		Next.~vnl_matrix_fixed();
		Previous.~vnl_matrix_fixed();

		vnl_symmetric_eigensystem<float> eigensystem(temp_return);

		float lamda[3]={eigensystem.get_eigenvalue(2), eigensystem.get_eigenvalue(1), eigensystem.get_eigenvalue(0)};

		for (int i=0;i<3;i++) {
			for (int j=i;j<3;j++) {
				if (abs(lamda[j])<abs(lamda[i])) {
					float temp = lamda[j];
					lamda[j] = lamda[i];
					lamda[i] = temp;
				}
			}
		}

		lamda1_iter.Set(lamda[0]);
		lamda2_iter.Set(lamda[1]);
		lamda3_iter.Set(lamda[2]);

		//float mean =(lamda[0]+lamda[1]+lamda[2])/3;
		//float variance = (vnl_math_sqr(lamda[0]-mean)+vnl_math_sqr(lamda[2]-mean)+vnl_math_sqr(lamda[1]-mean))/3;
		

		//std::cerr<<std::endl;
		//if (abs(lamda[1]/lamda[2])<.01 && abs(lamda[2])>.1) {
		//	//output_iter.Set(1);
		//	output_iter.Set(abs(lamda[1]/lamda[2]));
		//} else 
		
		//if (abs(lamda[0]/lamda[1])<.5 && abs(lamda[0])>.05 && lamda[2]<.1) {
		//	output_iter.Set(abs(lamda[0]/lamda[1]));
		//}
		
		if (abs(lamda[0])>20 && lamda[2]>0) {
			output_iter.Set(abs(lamda[1]/lamda[2]));		//lamda[0]/lamda[1]?
		} else {
			output_iter.Set(0);
		}
		
		if (max<output_iter.Get()) {
			max = output_iter.Get();
		}
		
			//if ((0.5<lamda_1/lamda_2 && lamda_1/lamda_2<2) && (lamda_1>0 && lamda_3/lamda_1<3)) {
				
			//}
		/*	else if (lamda_1>-0.5 && lamda_1<0.5 && (lamda_3-lamda_2)>100) {
				curvature_image_iter.Set(0);
			}*/
		//}

		count++;

		++input_iter;
		++output_iter;
		++x_iter;
		++y_iter;
		++z_iter;
		++lamda1_iter; ++lamda2_iter; ++lamda3_iter;
	}

	WriteITK(lamda1,"lamda1.hdr");
	WriteITK(lamda2,"lamda2.hdr");
	WriteITK(lamda3,"lamda3.hdr");

	WriteITK(output, "eig_test1.hdr");

	x_partial.~SmartPointer();
	y_partial.~SmartPointer();
	z_partial.~SmartPointer();
	for (output_iter.GoToBegin(); !output_iter.IsAtEnd(); ++output_iter) {
		if (output_iter.Get() > 0)
		{
			output_iter.Set(1-output_iter.Get()/max);
		}
	}

	return output;
}

ImageType::Pointer ComputeHessianResponse3(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	for (int k=0; k<3; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);
		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		float alpha = 0.25;
		float gamma = 0.5;

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		Himage_iter.GoToBegin();
		
		float max=0;
		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			// Get eigenvalues
			hessian_iter.Get().ComputeEigenValues(lambda);
			//std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);
			std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			float l[3];
			l[0] = lambda[0];
			l[1] = lambda[1];
			l[2] = lambda[2];

			Himage_iter.Set(S_fold(l, alpha, gamma));		

			count++;

			++hessian_iter;
			++Himage_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			WriteITK(LambdaImageVector[i], ss.str());
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_S_fold.hdr";
		WriteITK(Himage,ss.str());
	}

	return Himage;
}

void ComputeThinness(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);
	Himage->FillBuffer(0.0);

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	std::ofstream file;
	file.open("h.txt");

	for (int k=0; k<1; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);
		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		Himage_iter.GoToBegin();
		
		float max=0;
		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			// Get eigenvalues
			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

			
			//std::sort(lambda.Begin(),lambda.End(),OrderByValue);

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			float l[3];
			l[0] = lambda[0];
			l[1] = lambda[1];
			l[2] = lambda[2];
			
			if (l[2] < -100.0)
			{
				float t = thinness(l);

				file <<l[0]<<"\t"<<l[1]<<"\t"<<l[2]<<"\t"<<t<<"\n";

				Himage_iter.Set(t);

			}

			if ( count % 50000 == 0) {
				std::cout << "wait!" << std::endl;
			}

			count++;

			++hessian_iter;
			++Himage_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			WriteITK(LambdaImageVector[i], ss.str());
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_Thinness.hdr";
		WriteITK(Himage,ss.str());

		
	}
	file.close();
}

ImageType::Pointer ComputeHessianTest(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	ImageType::SizeType size = region.GetSize();

	// Create cone
	ImageType::Pointer cone = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex cone_iter(cone,region);

	for (cone_iter.GoToBegin(); !cone_iter.IsAtEnd(); ++cone_iter)
	{
		ImageType::IndexType idx = cone_iter.GetIndex();
		float x,y;
		x = ((float) idx[0]-size[0]/2) / size[0];
		y = ((float) idx[1]-size[0]/2) / size[1];

		float val = exp( -(vnl_math_sqr(x)+vnl_math_sqr(y)) / (2*0.5*0.5) );

		cone_iter.Set( val );
	}

	WriteITK(cone,"cone.hdr");

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	for (int k=1; k<2; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(cone);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);

		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		Himage_iter.GoToBegin();
		
		float max=0;
		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			// Get eigenvalues
			hessian_iter.Get().ComputeEigenValues(lambda);
			//std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);
			std::sort(lambda.Begin(),lambda.End(),OrderByValue);

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			float val = abs(lambda[1]/lambda[2]);
			if ( val > 0 && abs(lambda[0])>20 && lambda[2]>0)
			{
				Himage_iter.Set(val);

				if ( count == 0)
				{
					max = val;
				}

				if ( val > max )
				{
					max = val;
				}
			} else {
				Himage_iter.Set(0);
			}			

			count++;

			++hessian_iter;
			++Himage_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			WriteITK(LambdaImageVector[i], ss.str());
		}

		for (Himage_iter.GoToBegin(); !Himage_iter.IsAtEnd(); ++Himage_iter)
		{
			if (Himage_iter.Get() > 0)
			{
				Himage_iter.Set( 1 - (Himage_iter.Get()/max) );
			}
		}
	}

	return Himage;
}

bool OrderByMagnitude (double a, double b) {
	return ( abs(a) < abs(b) );
}

bool OrderByValue (double a, double b) {
	return ( a < b );
}

bool OrderByValueDesc (double a, double b) {
	return ( a > b );
}

double fA(EigenValueArrayType lambda, double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs(lambda[1]*lambda[2]));
	return exp(-vnl_math_sqr(Ra)/(2*vnl_math_sqr(alpha)));
}

double fB(EigenValueArrayType lambda, double beta, double gamma) {
	double Rb;
	Rb = abs(lambda[1]/lambda[2]);
	return exp(-vnl_math_sqr(Rb - gamma)/(2*vnl_math_sqr(beta)));
}

double fC(double ev1, double ev2, double eta) {
	double Rc;
	Rc = abs(ev1)/abs(ev2);
	return 1.0 - exp(-vnl_math_sqr(Rc)/(2*vnl_math_sqr(eta)));
}

double fRut(EigenValueArrayType lambda, double alpha, double beta, double gamma) {
	return fA(lambda, alpha)*fB(lambda, beta, gamma);
}

double fCup(EigenValueArrayType lambda, double eta) {
	return fC(lambda[0], lambda[1], eta)*fC(lambda[1], lambda[2], eta);
}

ImageType::Pointer ResampleImage(ImageType::Pointer input)
{
	ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
	
	typedef itk::AffineTransform< double, 3> TransformType;
	TransformType::Pointer transform = TransformType::New();
	resampleFilter->SetTransform(transform);

	/*
	typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	resampleFilter->SetInterpolator( interpolator );
	*/
	
	resampleFilter->SetDefaultPixelValue( -1024 );

	// Use isotropic spacing
	ImageType::SpacingType spacing = input->GetSpacing();
	spacing[1] = spacing[0];
	spacing[2] = spacing[0];

	resampleFilter->SetOutputSpacing( spacing );
	resampleFilter->SetOutputOrigin( input->GetOrigin() );
	resampleFilter->SetOutputDirection( input->GetDirection() );
	resampleFilter->SetSize( input->GetLargestPossibleRegion().GetSize() );
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
				if (i != n/2 )
				{
					VoxelTypeImage::IndexType idx = nit.GetIndex(i);
					
					if (idx[0] >= sidx[0] && idx[0] <= eidx[0] && idx[1] >= sidx[1] && idx[1] <= eidx[1] && idx[2] >= sidx[2] && idx[2] <= eidx[2])
					{
						
						if (nit.GetPixel(i) == Stool /*|| val > 400 */)
						{
							float val = input->GetPixel(idx);

							if ( val > max )
							{
								max = val;
							}	
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

		ss << "partial_" << i << "_" << name << ".hdr";
		
		WriteITK(partial_image, ss.str() );
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

float vesselness(const float lambda1, const float lambda2, const float lambda3, const float alpha, const float gamma12, const float gamma23)
{
	float vesselnessMeasure;
	if((0.0 >= lambda1) && (lambda1 > lambda2) && (lambda2 > lambda3))
	{
		vesselnessMeasure = (abs(lambda3) * pow(lambda2 / lambda3, gamma23)) * pow(1 + (lambda1 / abs(lambda2)), gamma12);
	}

	else if(((abs(lambda2) / alpha) > lambda1) && (lambda1 > 0.0) && (0.0 > lambda2) && (lambda2 > lambda3))
	{
		vesselnessMeasure = (abs(lambda3) * pow(lambda2 / lambda3, gamma23)) * pow(1 - (alpha * (lambda1 / abs(lambda2))), gamma12);
	}

	else
	{
		vesselnessMeasure = 0.0;
	}

	return vesselnessMeasure;
}

float thinness( float l[3] )
{
	//if ( l[0] < 0 && l[1] < 0 && l[2] < 0 )
	//{
		return exp((abs(l[0]/l[1]) - 1)*(abs(l[1])+abs(l[2])))*exp(-l[0]/l[2]);
	//} else {
	//	return 0;
	//}
}

void HessianMeasure(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	ImageType::SizeType size = region.GetSize();

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	for (int k=0; k<1; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);

		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		Himage_iter.GoToBegin();
		
		float max=0;
		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			// Get eigenvalues
			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);
			//std::sort(lambda.Begin(),lambda.End(),OrderByValue);

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			float val = abs(lambda[1]/lambda[2]);
			Himage_iter.Set( exp(-val) );

			count++;

			++hessian_iter;
			++Himage_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			//WriteITK(LambdaImageVector[i], ss.str());
		}
		
		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "ratio.hdr";
		WriteITK(Himage,ss.str());
		
	}
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

float weight(float ls, float lt, float alpha, float gamma)
{
	if ( (lt <= ls) && (ls <= 0 ) )
	{
		return pow(1 + (ls/abs(lt)),gamma);
	
	} else if ( (abs(lt)/alpha > ls) && (ls > 0) )
	{
		return pow(1-(alpha*ls/abs(lt)),gamma);
	} 

	return 0;
}

float omega(float ls, float lt, float gamma)
{
	if ( (lt <= ls) && (ls < 0) )
	{
		return pow(ls/lt,gamma);
	} 

	return 0;
}

float S_line(float lambda[3], float alpha, float gamma)
{
	if ( (lambda[2] <= lambda[1]) && (lambda[1] < 0) )
	{
		return abs(lambda[2])*omega(lambda[1],lambda[2],gamma)*weight(lambda[0],lambda[1],alpha,gamma);
	}

	return 0;
}

float S_blob(float lambda[3], float gamma)
{
	if ( (lambda[2] <= lambda[1]) && (lambda[1] <= lambda[0]) && (lambda[0] < 0) )
	{
		return abs(lambda[2])*omega(lambda[1], lambda[2], gamma)*omega(lambda[0], lambda[1], gamma);
	} 

	return 0;
}

float S_sheet(float lambda[3], float alpha, float gamma)
{
	if (lambda[2] < 0)
	{
		return abs(lambda[2])*weight(lambda[1], lambda[2], alpha, gamma)*weight(lambda[0], lambda[2], alpha, gamma);
	}

	return 0;
}

float S_fold(float lambda[3], float alpha, float gamma)
{
	if (lambda[2] < 0)
	{
		return abs(lambda[2])*weight(lambda[1], lambda[2], alpha, gamma)*weight(lambda[0], lambda[1], alpha, gamma);
	}
	
	return 0;
}

void ComputeSatoHessian(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	ImageType::SizeType size = region.GetSize();

	// Create response image
	ImageType::Pointer S = AllocateNewImage(region);
	

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	float alpha = 0.25;
	float gamma = 0.5;

	for (int k=1; k<4; k++)
	{
		
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		for (int j=0; j<3; j++) // loop over type: line, blob, sheet
		{
			S->FillBuffer(0.0);
			IteratorTypeFloat4WithIndex S_iter(S,region);

			hessian_iter.GoToBegin();
			S_iter.GoToBegin();
			
			int count=0;

			while (!hessian_iter.IsAtEnd()) 
			{
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				float l[3];
				l[0] = lambda[0]; l[1] = lambda[1]; l[2] = lambda[2];

				float val = 0;

				switch (j)
				{
				case 0:
					val = S_line(l, alpha, gamma);
					break;
				case 1:
					val = S_blob(l, gamma);
					break;
				case 2:
					val = S_sheet(l, alpha, gamma);
					break;
				}

				S_iter.Set(val);

				count++;

				++hessian_iter;
				++S_iter;

			}

			std::string type;
			switch (j)
			{
			case 0: type="line";break;
			case 1: type="blob";break;
			case 2: type="sheet";break;
			}

			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_type_" << type << ".hdr";
			WriteITK(S,ss.str());
		}
	}
}

void ComputeFrangiHessian(ImageType::Pointer input)
{
	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	ImageType::SizeType size = region.GetSize();

	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	for (int k=0; k<3; k++)
	{
		
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();

		for (int j=0; j<3; j++) // loop over type: blobs, vessels, plates
		{
			ObjectnessFilterType::Pointer ofilter = ObjectnessFilterType::New();
			ofilter->SetInput(hessianFilter->GetOutput());
			ofilter->SetBrightObject(false);
			ofilter->SetObjectDimension(j);

			std::string type;
			switch (j)
			{
				case 0: type="blob";break;
				case 1: type="vessels";break;
				case 2: type="plates";break;
			}

			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_type_" << type << ".hdr";
			WriteITK(ofilter->GetOutput(),ss.str());
		}
	}
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
		ss << "fuzzy_scene_mean" << mean << "_variance_" << variance[i] << "_Tmax_" << Tmax[i] << ".hdr";

		WriteITK(rescaler->GetOutput(), ss.str());
	}

	/*

	typedef itk::ImageFileWriter< FuzzySceneType > FuzzyWriterType;
	FuzzyWriterType::Pointer fwriter = FuzzyWriterType::New();
	fwriter->SetFileName("fuzzy_scene.hdr");
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
	WriteITK(smoother->GetOutput(),"input_sub_smoothed.hdr");

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
	ss << "connected_threshold_lower" << lower << "_upper_" << upper << ".hdr";
	WriteITK(connecter->GetOutput(),ss.str());

}