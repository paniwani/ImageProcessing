#include "RemoveStool3.h"

float Modified=false;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
double beta=.7;
double weight_sum=2.5;
int writeCount=1;
boost::xpressive::sregex rex = boost::xpressive::sregex::compile( "^2*[5-7]+1+$" );

int main(int argc, char * argv[])
{
	//-------------------------------------------BEGIN SETUP-------------------------------------------------------

	std::cerr<<"Started"<<std::endl;

	// Use parameter changes
	Modified = true;

	// Read input
	ImageType::Pointer input = ImageType::New();
	ReadITK(input, "C:/ImageData/mr10_092_13p.i0344_100-105.hdr");

	// Set region
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();
	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);

	//declares a temporary storage of the image with the selected region
	ImageType::Pointer input_temp = AllocateNewImage(fullRegion);
	input_temp->SetSpacing(input->GetSpacing());

	//declares two iterators, one for the original image, one for the sub image
	IteratorTypeFloat4WithIndex input_temp_iter(input,fullRegion);
	IteratorTypeFloat4WithIndex input_iter(input_temp,fullRegion);

	//Copy Image the sub region into the temporary container
	for(input_iter.GoToBegin() , input_temp_iter.GoToBegin() ; !input_iter.IsAtEnd() && !input_temp_iter.IsAtEnd(); ++input_iter, ++input_temp_iter) {
		input_iter.Set(input_temp_iter.Get());
	}

	//Print out the sub region of the image, un modified
	WriteITK(input_temp,"input_temp.hdr");

	//Creates an image to store the type of voxel for each voxel in the sub image
	VoxelTypeImage::Pointer voxel_type = VoxelTypeImage::New();
	voxel_type->SetRegions(fullRegion);
	voxel_type->Allocate();

	IteratorTypeVoxelType voxel_type_iter(voxel_type,fullRegion);

	//computes the gradient magnitude for the image and set it to a temporary variable
    GradientMagnitudeFilterType::Pointer gradient_filter = GradientMagnitudeFilterType::New();
    gradient_filter->SetInput(input_temp);
    gradient_filter->Update();
    ImageType::Pointer temp = gradient_filter->GetOutput();					//this variable will be reused later for other purpose, meaning change will be commented
	temp->SetSpacing(input->GetSpacing());
	gradient_filter.~SmartPointer();										//Deleting Gradient Filter

	//defines an iterator for the gradient magnitude
	IteratorTypeFloat4WithIndex temp_iter(temp,fullRegion);
	
	//writes out the gradient magnitude for the input image
	WriteITK(temp, "gradient.hdr");

	//stores the largest gradient magnitude, this is purely for reference on the validity of the image
	float gradient_max=0;

	//Computes an initial hard threshold for the voxel type and sets initial partials
	for(input_iter.GoToBegin(), temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !input_iter.IsAtEnd() && !temp_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd();
        ++input_iter, ++temp_iter, ++voxel_type_iter) 
	{
		//updates for maximum gradient
		if (gradient_max<temp_iter.Get()) { gradient_max=temp_iter.Get(); }

		//computes the hard threshold classification for a voxel in the image: Air/Tissue/Stool/Unclassified
		//std::cerr << "input_iter.Get(): " << input_iter.Get() << std::endl;
        VoxelType type = SingleMaterialClassification(input_iter.Get(), temp_iter.Get());
		voxel_type_iter.Set(type);		//sets the type in the voxel_type structure
    }

	WriteITK(voxel_type, "voxel_type_start.hdr");

	// outputs the gradient magnitude purely for reference
	std::cerr<<"Gradient max: "<<gradient_max<<std::endl;

	//Filling in Stool Holes
	for (voxel_type_iter.GoToBegin(), temp_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++temp_iter)
    {
        //creates the morphological images to fill holes
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
					voxel_type_iter.Set(Stool);			//Modify unclassified into stool holes
                }
				temp_iter.Set(0);
                break;
            case Air:
                temp_iter.Set(-1);
                break;
        }
    }

	WriteITK(voxel_type, "voxel_type_stool_closed.hdr");
	WriteITK(temp, "temp_tissue.hdr");

	closingFilter = ClosingFilterType::New();
	closingFilter->SetKernel(structuringElement);
	closingFilter->SetInput(temp);
	closingFilter->Update();
	temp=closingFilter->GetOutput();
	temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	closingFilter.~SmartPointer();

	WriteITK(temp, "temp_tissue_closed.hdr");

	ByteImageType::Pointer chamfer_colon=AllocateNewByteImage( fullRegion );
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);

    // Updates Information after morphological operations
    for (chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(), temp_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++chamfer_colon_iter, ++temp_iter) 
    {
        switch(voxel_type_iter.Get()) {
            case Unclassified:
                if (temp_iter.Get()==1) {					// unclassified into tissue
					voxel_type_iter.Set(Tissue);		
					chamfer_colon_iter.Set(0);
				} else {
					chamfer_colon_iter.Set(1);
				}
                break;
			case Tissue:
				chamfer_colon_iter.Set(0);
				break;
			default:
				chamfer_colon_iter.Set(1);
				break;
        }
    }

	WriteITK(voxel_type, "voxel_type_tissue_closed.hdr");
	
	WriteITK(chamfer_colon, "chamfer_colon_air_pre.hdr");

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
	
	labelVector.erase( labelVector.begin() );
	sort(labelVector.begin(), labelVector.end(), compareSize);
	int colonLabel = labelVector[0]->GetLabel();

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

	WriteITK(chamfer_colon,"chamfer_colon_no_bkg.hdr");

	// Smooth colon mask using median filter
	BinaryMedianFilterType::Pointer medianFilter = BinaryMedianFilterType::New();
	medianFilter->SetInput(chamfer_colon);
	medianFilter->SetForegroundValue(1);
	medianFilter->Update();
	chamfer_colon=medianFilter->GetOutput();
	medianFilter.~SmartPointer();
	chamfer_colon_iter=IteratorTypeByteWithIndex(chamfer_colon,fullRegion);

	WriteITK(chamfer_colon,"chamfer_colon_median.hdr");

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
	WriteITK(chamfer_colon, "distance_map_to_air.hdr");

	chamfer_colon_iter.GoToBegin();
	chamfer_colon_iter.Set(-1);	//first pixel == -1

	for (chamfer_colon_iter.GoToBegin(), temp_iter.GoToBegin(), input_iter.GoToBegin();
        !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd() && ! input_iter.IsAtEnd();  
        ++chamfer_colon_iter, ++temp_iter, ++input_iter) 
    {
		if (chamfer_colon_iter.Get()<20 && chamfer_colon_iter.Get()>-1) {
			chamfer_colon_iter.Set(1);
		} else {
			chamfer_colon_iter.Set(0);
		}
    }

	WriteITK(chamfer_colon, "chamfer_colon_air_mask.hdr");

	// Read chamfer colon
	//ReadITK(chamfer_colon, "mr10_092_13p.i0344_100-105_chamfer_colon_air_mask.hdr");
	//chamfer_colon_iter=IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	//WriteITK(chamfer_colon,"chamfer_colon.hdr");
	

	//-------------------------------------------END SETUP-------------------------------------------------------

	//-------------------------------------------BEGIN QR--------------------------------------------------------

	//ReadITK(voxel_type, "C:/ImageData/mr10_092_13p.i0344_100-105_voxel_edge_optimized.hdr");
	//voxel_type_iter = IteratorTypeVoxelType(voxel_type,fullRegion);
	//WriteITK(voxel_type, "voxel_input.hdr");

	
	
    // get cubic bspline image interpolation
    InterpolationType::Pointer input_interpolator = InterpolationType::New();
    input_interpolator->SetSplineOrder(3);
    input_interpolator->SetInputImage(input);
	// Storage for Smax
	ImageType::Pointer input_smax = AllocateNewImage(fullRegion);

	// Load smax
	//ReadITK(input_smax, "mr10_092_13p.i0344_100-105_smax.hdr");
	//WriteITK(input_smax, "smax.hdr");

	IteratorTypeFloat4WithIndex input_smax_iter(input_smax,fullRegion);

	std::cerr << "Running voxel edge classification..." << std::endl;

	std::ofstream myfile;
	myfile.open("gradient.txt");

	//Runs edge classification with 3 different increments inside the colon only
    for (voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), temp_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
		!voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !temp_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++voxel_type_iter, ++input_smax_iter, ++temp_iter, ++chamfer_colon_iter) {
		
		if (chamfer_colon_iter.Get() == 1 && voxel_type_iter.Get()==Unclassified) {
			
			//Initialization
			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
			float temp_threshold=-1;
			VoxelType temp_type = Unclassified;
			VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,input_interpolator,input_smax_iter,voxel_index);
			VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,input_interpolator,input_smax_iter,voxel_index);
			VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,input_interpolator,input_smax_iter,voxel_index);

			voxel_type_iter.Set(temp_type);
		}
    }

	myfile.close();

	std::cerr<<"Done Edge"<<std::endl;

	WriteITK(voxel_type,"voxel_edge_class.hdr");
	WriteITK(input_smax,"smax.hdr");

	// Optimize SA transitions using gradient information
	std::cerr<<"Running voxel edge optimization"<<std::endl;
	OptimizeVoxelEdge(input, input_iter, voxel_type, voxel_type_iter);
	
	WriteITK(voxel_type, "voxel_edge_optimized.hdr");

	// Deletes the Interpolators
	input_interpolator.~SmartPointer();

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

	WriteITK(chamfer_from_stool_involvement,"non_stool_mask.hdr");

	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_from_stool_involvement);
	int weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_from_stool_involvement=chamfer_filter->GetOutput();
	chamfer_from_stool_involvement_iter = IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	chamfer_filter.~SmartPointer();

	WriteITK(chamfer_from_stool_involvement, "non_stool_distance.hdr");

	// Find thin stool voxels <=8 and no neighbors >=8
	int chamfercutoff = 8;
	int thinStoolCount = 0;

	for(voxel_type_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
		!voxel_type_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd(); 
		++voxel_type_iter, ++chamfer_from_stool_involvement_iter) 
    {
		
		if (chamfer_from_stool_involvement_iter.Get() <= chamfercutoff && chamfer_from_stool_involvement_iter.Get() > 0)
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
									
									if (chamfer_from_stool_involvement->GetPixel(neighbor_index) >= chamfercutoff) {
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

	std::cout << "Number of thin stool voxels: " << thinStoolCount << std::endl;

	WriteITK(voxel_type,"voxel_type_thinstool.hdr");










	//-------------------------------------------END QR-----------------------------------------------------

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


	//-------------------------------------------BEGIN EM--------------------------------------------------------

 //   int counter=0;							//used to total certain values (description given when set)
 //   double mean[3]={0, 0, 0};				//stores the mean of the air/tissue/stool classes
 //   double sum[3]={0, 0, 0};				//stores the total # of partials of the three classes
 //   double variance[3]={0, 0, 0};			//stores the variance of the air/tissue/stool classes
	//float weight[3]={0,0,0};				//stores the weights of the air/tissue/stool classes

	////Declares a structure to store the air/tissue/stool partial volumes for a voxel in the image
	//ImageVectorType::Pointer partialVector2 = ImageVectorType::New();
 //   partialVector2->SetRegions(fullRegion);
 //   partialVector2->Allocate();
	//partialVector2->SetSpacing(input->GetSpacing());

	////defines an iterator over the partial volume structure
 //   IteratorImageVectorType partialVector_iter(partialVector2,fullRegion);

	////computes the gradient magnitude for the image and set it to a temporary variable
 //   gradient_filter = GradientMagnitudeFilterType::New();
 //   gradient_filter->SetInput(input_temp);
 //   gradient_filter->Update();
 //   temp = gradient_filter->GetOutput();					//this variable will be reused later for other purpose, meaning change will be commented
	//temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	//gradient_filter.~SmartPointer();	

	//std::cerr<<"Allocated All Spaces"<<std::endl;

	////computes the total number of voxels in the sub image
	//int total_vertices=0;
	//for(partialVector_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin();
 //       !partialVector_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd(); 
 //       ++partialVector_iter, ++input_iter, ++chamfer_colon_iter, ++voxel_type_iter) 
	//{
	//	total_vertices++;
	//	VoxelType type=voxel_type_iter.Get();
	//	float value[3]={.1,.45,.45};		//declares default partials for unclassified class
	//	
	//	// air/tissue/stool 0/1/2
	//	switch (type) {
	//		case Tissue:
	//			value[0]=0;
	//			value[1]=1;					//sets the partial to all tissue
	//			value[2]=0;
	//			mean[1]+=input_iter.Get();	//updates the tissue mean
	//			sum[1]++;					//updates the total partial count
	//			break;
	//		case Stool:
	//			value[0]=0;
	//			value[1]=0;
	//			value[2]=1;					//sets the partial to all stool
	//			mean[2]+=input_iter.Get();	//updates the stool mean
	//			sum[2]++;					//updates the total stool partial
	//			break;
	//		case Air:
	//			value[0]=1;					//sets the partial to all Air
	//			value[1]=0;
	//			value[2]=0;
	//			mean[0]+=input_iter.Get();	//updates the air mean
	//			sum[0]++;					//updates the total air partial
	//			break;
	//		case Unclassified:
	//			break;
	//	}   

	//	CovariantVectorType data(value);	//sets the partial to be stored in the actual structure
	//	partialVector_iter.Set(data);		//stores the partial in the structure
 //   }

	//std::cerr<<"All partials Initialized"<<std::endl;

	////computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
	//for (int i=0;i<3;i++) { mean[i]=mean[i]/sum[i]; }    
 //
 //   for(input_iter.GoToBegin(), voxel_type_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
	//	!input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
	//	++input_iter, ++voxel_type_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	//{
	//	if (chamfer_colon_iter.Get()==1) {		
	//		CovariantVectorType data = partialVector_iter.Get();		//gets the partials
	//		float Z[3]={data[0],data[1],data[2]};

	//		for (int i=0;i<3;i++) {
	//			variance[i]+=Z[i]*vnl_math_sqr(input_iter.Get()-mean[i]);	//computes the variance sum(partial[i]*(value[i]-mean[i])^2)/sum(partial[i])
	//		}
	//	}
 //   }

	////computes the weights and varaince
	//for (int i=0;i<3;i++) {
	//	variance[i]=variance[i]/sum[i];
	//	weight[i]=sum[i]/total_vertices;
	//}

	////outputs the initial mean and variance of each class
	//std::cerr<<std::endl;
	//std::cerr<<"EM0"<<std::endl;
	//std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	//std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	//std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;


	////declares an image used to store the output image at the end of the i-th iteration
	//ImageType::Pointer em_output_image = ImageType::New();
	//em_output_image->SetRegions(fullRegion);
	//em_output_image->Allocate();
	//em_output_image->SetSpacing(input->GetSpacing());
	//IteratorTypeFloat4WithIndex em_output_image_iter(em_output_image,fullRegion);

	////used 20 in original test case
	////computes 5 iterations of the 
	//// Maximization Algorithm
 //   for(int i=0;i<5;i++) {
	//	//gives the temporary storage of the variables corresponding to sum, variance, and mean for each class on the i+1 iteration
	//	double sum_temp[3]={0,0,0};
 //       double variance_temp[3]={0,0,0};
	//	double mean_temp[3]={0,0,0};

	//	for(input_iter.GoToBegin(), voxel_type_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), temp_iter.GoToBegin(); 
	//		!input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd(); 
	//		++input_iter, ++partialVector_iter, ++voxel_type_iter, ++chamfer_colon_iter, ++temp_iter) 
	//	{
	//		CovariantVectorType data = partialVector_iter.Get();		//retrieves the partial informations
	//		float Z[3]={data[0],data[1],data[2]};

	//		vnl_vector<float> Z_update=expectation(input_iter.Get(),mean, variance, weight, GetNeighbor(partialVector2,input_iter.GetIndex()), Z, temp_iter.Get());
	//		data[0]=Z_update[0];										//updates the partial values
	//		data[1]=Z_update[1];
	//		data[2]=Z_update[2];
	//		partialVector_iter.Set(data);

	//		//updates the new mean, variance, and total partial sum for each class accordingly
	//		for (int k=0;k<3;k++) {
	//			mean_temp[k]+=data[k]*input_iter.Get();	
	//			variance_temp[k]+=data[k]*vnl_math_sqr(input_iter.Get()-mean[k]);
	//			sum_temp[k]+=data[k];
	//		}
 //       }

	//	//updates the mean, variance, and weights for each class
	//	for (int k=0;k<3;k++) {
	//		mean[k]=mean_temp[k]/sum_temp[k];
	//		variance[k]=variance_temp[k]/sum_temp[k];
	//		weight[k]=sum_temp[k]/total_vertices;
	//	}
	//	std::cerr<<std::endl;
	//	std::cerr<<"EM"<<i+1<<std::endl;
	//	std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	//	std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	//	std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;
 //   }

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
			return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance))*exp(-vnl_math_sqr(current_partial-local_mean)/(2*local_variance))/(sqrt(2*PI*local_variance));
		} else {
			return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance));
		}
	} else {
		return 0;
	}
}

vnl_vector<float> expectation(double Y, double mean[], double variance[], float weight[],  vnl_matrix<float> neighbor, float current_partial[], float gradient) {
	int nKernal =3;
	float pFK[3]={0,0,0};	//probabilities of each class: air/tissue/stool per voxel
	float pKF[4]={0,0,0};   //""	"", [3] = sum of probabilities of air/tissue/stool
	float sumPFK=0.0;
	for (int i=0;i<3;i++) {
		double total_neighbor=0;
		
		if (current_partial[0]<.6 && current_partial[1]<.6 && current_partial[2]<.6) {
			beta=.7;
		} else {
			beta=0;
		}
		
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

                                      InterpolationType::Pointer input_interpolator, 

                                      IteratorTypeFloat4WithIndex input_smax,

                                      ImageType::IndexType index) 

{

	//if (index[0] == 90 && index[1] == 512-274 && index[2] == 2) {
	//	std::cout<<"stop"<<std::endl;
	//}

	ImageType::RegionType fullRegion = input_interpolator->GetInputImage()->GetLargestPossibleRegion();
    InterpolationType::ContinuousIndexType offset_index;
    offset_index[0]=index[0];
    offset_index[1]=index[1];
    offset_index[2]=index[2];
	ImageType::IndexType endIndex = input_interpolator->GetEndIndex();
	ImageType::IndexType startIndex = input_interpolator->GetStartIndex();
    InterpolationType::CovariantVectorType gradient = input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index);
    
	double norm = gradient.GetNorm();
	if (norm>0)
	{
		for(int i=0; i<3;i++)
		{
			gradient[i] /= norm;

			if ( fabs( gradient[i] ) < 1E-5 )
				gradient[i] = 0;
		}
	} 

	//myfile << "index: (" << index[0] << " " << 512-index[1] << " " << index[2] << ") grad: (" << gradient[0] << " " << -gradient[1] << " " << gradient[2] << ")\n";

    //stores the intensity
    float temp_intensity[5];
    float temp_gradient_magnitude[5];
	//stores the target voxel
    temp_intensity[2]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[2]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //+d1
    offset_index[0]=index[0]+gradient[0]*d1;
    offset_index[1]=index[1]+gradient[1]*d1;
    offset_index[2]=index[2]+gradient[2]*d1;
    temp_intensity[3]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[3]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //-d1
    offset_index[0]=index[0]-gradient[0]*d1;
    offset_index[1]=index[1]-gradient[1]*d1;
    offset_index[2]=index[2]-gradient[2]*d1;
    temp_intensity[1]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[1]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //+d2
    offset_index[0]=index[0]+gradient[0]*d2;
    offset_index[1]=index[1]+gradient[1]*d2;
    offset_index[2]=index[2]+gradient[2]*d2;
    temp_intensity[4]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[4]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //-d2
    offset_index[0]=index[0]-gradient[0]*d2;
    offset_index[1]=index[1]-gradient[1]*d2;
    offset_index[2]=index[2]-gradient[2]*d2;
    temp_intensity[0]=input_interpolator->EvaluateAtContinuousIndex(offset_index);
    temp_gradient_magnitude[0]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
	//std::cerr<<"vector created"<<std::endl;
	float min_distance=*threshold;
	float stool_tissue_Smax;
	float stool_air_Smax;
	//Tissue Stool
	stool_tissue_Smax=ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
	stool_air_Smax=Stool_Air_ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
    float distanceTS=AverageTissueStoolDist(stool_tissue_Smax, temp_intensity,temp_gradient_magnitude);
	float distanceSA=AverageStoolAirDist(stool_air_Smax, temp_intensity,temp_gradient_magnitude);
	float max_gradient =0;
	
	for (int i=0;i<5;i++) {
		if (max_gradient<temp_gradient_magnitude[i]) {
			max_gradient=temp_gradient_magnitude[i];
		}
	}

	//for (int i=0;i<5;i++) {
	//	temp_gradient_magnitude[i]=temp_gradient_magnitude[i]/max_gradient*1000;
	//	
	//}


	float distanceTA=AverageTissueAirDist(temp_intensity,temp_gradient_magnitude);

	

	/*bool containPositive=false;
	bool containNegative=false;
	for (int i=0;i<5;i++) {
		if (50<temp_intensity[i]) {
			containPositive=true;
		} else if (temp_intensity[i]<-50) {
			containNegative=true;
		}
	}*/

	float smax=0;
	float distance=0;
	VoxelType voxel_type=Unclassified;
	
	if (distanceSA<=distanceTS && distanceSA<=distanceTA){
		//if (containNegative && !containPositive) {
		//	distance=distanceTA;
		//	voxel_type=TissueAir;
		//	smax=stool_air_Smax;
		//} else if (containPositive && !containNegative) {
		//	distance=distanceTS;
		//	voxel_type=TissueStool;
		//	smax=stool_tissue_Smax;
		//} else {
		//if (stool_air_Smax>600 || !Modified) {
			distance=distanceSA;
			voxel_type=StoolAir;
			smax=stool_air_Smax;
		//} else {
		//	distance=distanceSA;
		//	voxel_type=TissueStool;
		//	smax=stool_tissue_Smax;
		//}
		//}
	} else if (distanceTS<=distanceTA && distanceTS<=distanceSA ) {
		distance=distanceTS;
		voxel_type=TissueStool;
		smax=stool_tissue_Smax;
	} else {
		smax=stool_air_Smax;
		distance=distanceTA;
		voxel_type=TissueAir;
	}

	/*
	if (smax<150 && Modified) {
		smax=stool_air_Smax;
		distance=distanceTA;
		voxel_type=TissueAir;
	}
	*/


    if (min_distance>=distance || min_distance==-1 || *previous==Unclassified) {
		*threshold=distance;
        *previous=voxel_type;
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
	double coefficients[3]={-4.0/1000,-4.0,0};;
	
	/*if (Modified) {
		coefficients[0]= -6.0/1000;
		coefficients[1]= -6.0;
		coefficients[2]=	0;
	}*/

	if (Modified) {
		coefficients[0]= -2.4/1000;
		coefficients[1]= -2.4;
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
    double coefficients[3] ={-4/(Smax+1000),-4/(Smax+1000)*(1000-Smax),4/(Smax+1000)*1000*Smax};
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
	/*if (Modified) {
		if ((input_pixel >=400 && input_gradient_pixel<=0.8*input_pixel) || input_pixel >= 1000) {
			return Stool;
		} else if (input_pixel<=-800 && input_gradient_pixel<=250) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
			return Tissue;
		} else {
			return Unclassified;
		}
	} else {*/
		if ((input_pixel >=180 && input_gradient_pixel<=0.8*input_pixel)) {
			return Stool;
		} else if (input_pixel<=-800 && input_gradient_pixel<=250) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
			return Tissue;
		} else {
			return Unclassified;
		}
	//}
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

void OptimizeVoxelEdge(ImageType::Pointer input, IteratorTypeFloat4WithIndex inputIt, VoxelTypeImage::Pointer &voxelEdge, IteratorTypeVoxelType &voxelEdgeIt ) {
	// PARAMETERS
	int numOfVoxels = 10; // to move in direction of gradient

	// Pattern matching
	char *regex = "^2*[5-7]+1+$";

	// Setup text file
	std::ofstream myfile;
	myfile.open ("OptimizeVoxelEdgeOutput.txt");

	myfile << "Number of voxels to move in direction of gradient: " << numOfVoxels << "\n";
	myfile << "Pattern: " << regex << "\n\n";

	// Get subregions
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

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
	
	// Calculate gradient
	GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->Update();
	ImageVectorType::Pointer gradient = gradientFilter->GetOutput();
	IteratorImageVectorType gradientIt( gradient, region );

	/*
	// Test gradient
	for (gradientIt.GoToBegin(); !gradientIt.IsAtEnd(); ++gradientIt)
	{
		ImageVectorType::IndexType idx = gradientIt.GetIndex();
		CovariantVectorType grad = gradientIt.Get();

		// Normalize gradient vector
		float norm = grad.GetNorm();
		if (norm>0)
		{
			for(int i=0; i<3;i++)
			{
				grad[i] /= norm;
			}
		}

		//grad.Normalize();

		// Output gradient from filter to text file
		//myfile << "Index: (" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t" << inputIt.Get() << "\t";
		myfile << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\n";

		//std::cout << "Index: (" << idx[0] << ", " << 511-idx[1] << ", " << idx[2] << ")\t" << inputIt.Get() << "\t";
		//std::cout << "Gradient: (" << grad[0] << ", " << grad[1] << ", " << grad[2] << ")\t\n";
	}
	myfile.close();
	*/
	
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
			
			float norm = grad.GetNorm();
			if (norm>0)
			{
				for(int i=0; i<3;i++)
				{
					grad[i] /= norm;
				}
			}

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
				boost::xpressive::smatch match;

				if( boost::xpressive::regex_match( vs, match, rex ) )
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
		idx[2] += round( grad[2] );

		// Check bounds
		for (int j=0; j < 3; j++)
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


