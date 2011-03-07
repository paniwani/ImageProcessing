#include "LaxativeFreePrep.h"


#include <itkImageFileReader.h>

float Modified=false;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
double beta=.7;
double weight_sum=2.5;

double Probability(double Y, double mean, double variance,  double current_partial, double local_variance, double local_mean) {  
	if (variance>0.01) {
		if (local_variance>0.01) {
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
	float pFK[3]={0,0,0};
	float pKF[4]={0,0,0};
	float sumPFK=0.0;
	for (int i=0;i<3;i++) {
		double total_neighbor=0;
		if (current_partial[0]<.6 && current_partial[1]<.6 && current_partial[2]<.6) {
			beta=.7;
		} else {
			beta=0;
		}
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
		pFK[i]=Probability(Y,mean[i],variance[i], current_partial[i], local_variance, local_mean)*weight[i];
		sumPFK+=pFK[i];
	}
	for(int i=0;i<3;i++) {
		pKF[i]=pFK[i]/sumPFK;
	}
	pKF[3]=sumPFK;
	vnl_vector<float> return_value(4, 4, pKF);
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
					//	voxel_type_iter.Set(Air);        //sets the voxel type
					} else if (data[1]>.8 && temp_iter.Get()<300*data[1] && (input_iter.Get()-1024)<600) {
						voxel_type_iter.Set(Tissue);        //sets the voxel type
					} else if (data[2]>.9 && temp_iter.Get()<0.8*(input_iter.Get()-1024)) {
					//	voxel_type_iter.Set(Stool);        //sets the voxel type
					} else {
						//voxel_type_iter.Set(Unclassified);        //sets the voxel type
					}
					break;
				case Tissue:
					if (data[1]<.1 /*&& temp_iter.Get()>300*/) {
						voxel_type_iter.Set(Unclassified);
					}
					break;
			}
		}
		switch (voxel_type_iter.Get()) {
			case Unclassified:
				temp_iter.Set(0);
				break;
			case Stool:
				temp_iter.Set(-2);
				break;
			case Tissue:
				temp_iter.Set(1);
				break;
			case Air:
				temp_iter.Set(-1);
				break;
		}
	}

}

ImageType::Pointer longruntests(ImageType::Pointer temp, int size, int direction, float modifier) {
	ImageType::RegionType fullRegion=temp->GetLargestPossibleRegion();
	temp->SetRequestedRegion(fullRegion);
	
	ImageType::Pointer output = AllocateNewImage(fullRegion);

	IteratorTypeFloat4WithIndex temp_iter(output, fullRegion);
	IteratorTypeFloat4WithIndex input_iter(temp, fullRegion);

	for(input_iter.GoToBegin(), temp_iter.GoToBegin(); !input_iter.IsAtEnd() && !temp_iter.IsAtEnd();++input_iter, ++temp_iter) {
		int count=0;
		float mean=0;
		ImageType::IndexType index = input_iter.GetIndex();
		
		float max=0;
		floatDynArray data((2*size+1)*(2*size+1));
		int x_start=-size-1;
		int y_start=-size-1;

		int x_end=x_start;
		int y_end=y_start;

		for(int i=-size;i<=size;i++) {
			ImageType::IndexType temp_index;
			for (int j=-size;j<=size;j++) {
				switch(direction) {
					case 0:
						//yz
						temp_index[0]=index[0];
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2]+j;
						break;
					case 1:
						//xz
						temp_index[0]=index[0]+i;
						temp_index[1]=index[1];
						temp_index[2]=index[2]+j;
						break;
					case 2:
						//xy
						temp_index[0]=index[0]+j;
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2];
						break;
					case 3:
						//z 45%
						temp_index[0]=index[0]+i;
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2]+j;
						break;
					case 4:
						//z -45%
						temp_index[0]=index[0]-i;
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2]+j;
						break;
					case 5:
						//x 45%
						temp_index[0]=index[0]+j;
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2]+i;
						break;
					case 6:
						//x -45%
						temp_index[0]=index[0]+j;
						temp_index[1]=index[1]+i;
						temp_index[2]=index[2]-i;
						break;
					case 7:
						//y 45%
						temp_index[0]=index[0]+i;
						temp_index[1]=index[1]+j;
						temp_index[2]=index[2]+i;
						break;
					case 8:
						//y -45%
						temp_index[0]=index[0]-i;
						temp_index[1]=index[1]+j;
						temp_index[2]=index[2]+i;
						break;
				}

				if (fullRegion.IsInside(temp_index)) {
					if (y_start<-size) {
						y_start=j;
					}
					if (x_start<-size) {
						x_start=i;
					}
					data.SetAt(count,temp->GetPixel(temp_index));
					if (max<data.GetAt(count)){
						max = data.GetAt(count);
					}
					if (data.GetAt(count)==-1) {
						data.SetAt(count,0);
					}
					count++;
					if (i>x_end) {
						x_end=i;
					}
					if (j>y_end) {
						y_end=j;
					}
					//std::cerr<<data.GetAt(i)<<std::endl;
				}
			}
		}

		if (count>0) {
			//std::cerr<<"short: "<<count<<" "<<x_end-x_start<<" "<<y_end-y_start<<std::endl;
			data.SetSize(count);

			shortDynArray short_data(data.GetSize());
			for (int i=0;i<count;i++) {
				short_data.SetAt(i, data.GetAt(i)*255.0/max);
			//	std::cerr<<data.GetAt(i)<<std::endl;
			}
			TextureFeature_RunLength texture_run_length;

			//
			//int returned_value = NIH_Algo_Texture_RunLength_2D(short_data,x_end-x_start+1,y_end-y_start+1,0,255,15,texture_run_length,-1);
			TextureFeature_StatisticalMoment texture_stat_momment;

			int returned_value =NIH_Algo_Texture_StatisticalMoments_2D(short_data,x_end-x_start+1, y_end-y_start+1,texture_stat_momment,-1);
			if (returned_value==0) {
				temp_iter.Set(abs(texture_stat_momment.variance)*modifier);
			} else {
				std::cerr<<"fucked up somewher"<<std::endl;
				break;
			}
		}
	}

	return output;
}

ImageType::Pointer levelset_comp(ImageType::Pointer inputImage) {

	//WriteITK(inputImage, "level set input.mhd");
	//WriteITK(chamfer, "level set chamfer.mhd");


	//typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType,
 //   ImageType > GradientImageType;

	//GradientImageType::Pointer gradMagnitude = GradientImageType::New();
	//gradMagnitude->SetInput( inputImage );
	//gradMagnitude->SetSigma( 1.0 );

	//typedef itk::SigmoidImageFilter< ImageType, ImageType >
 //   SigmoidFilterType;
	//SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	//sigmoid->SetOutputMinimum( 0.0 );
	//sigmoid->SetOutputMaximum( 1.0 );
	//sigmoid->SetAlpha( -0.4 );
	//sigmoid->SetBeta( 200 );
	//sigmoid->SetInput( inputImage );
	//sigmoid->Update();

	//WriteITK(sigmoid->GetOutput(), "sigmoid.mhd");

	LocalRoughnessImageFilterType::Pointer test = LocalRoughnessImageFilterType::New();
	test->SetInput(inputImage);
	test->Update();
	ImageType::Pointer edge_potential=test->GetOutput();


	ImageType::Pointer marching = AllocateNewImage(inputImage->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex marching_iter(marching,marching->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_iter(inputImage,inputImage->GetLargestPossibleRegion());
	
	for (marching_iter.GoToBegin(), input_iter.GoToBegin(); !marching_iter.IsAtEnd() && !input_iter.IsAtEnd(); ++marching_iter, ++input_iter) {
		if (input_iter.Get()<300) {
			marching_iter.Set(15);
		} else {
			marching_iter.Set(-15);
		}
	}

	typedef itk::CurvesLevelSetImageFilter< ImageType, ImageType > CurvesFilterType;
  
	CurvesFilterType::Pointer curvesFilter = CurvesFilterType::New();
  
	// set the initial level set
	curvesFilter->SetInput( marching );

//   set the edge potential image
	curvesFilter->SetFeatureImage( edge_potential );

  // set the weights between the propagation, curvature and advection terms
	curvesFilter->SetPropagationScaling( 0.1 );
	curvesFilter->SetCurvatureScaling( 0.5 );
	curvesFilter->SetAdvectionScaling( 0.4 );

  // set the convergence criteria
	curvesFilter->SetNumberOfIterations(50);
  
	curvesFilter->SetNumberOfThreads(2);
	curvesFilter->Update();

	
	


	return curvesFilter->GetOutput();
}

ImageType::Pointer RemoveStool2(ImageType::Pointer input, int x_offset, int y_offset, int z_offset, int x_size, int y_size, int z_size, SSlice*& DAB) {

	std::cerr<<"Started"<<std::endl;

	Modified=true;

	ImageType::SizeType size = {x_size, y_size, z_size};
	ImageType::IndexType index = {x_offset,y_offset,z_offset};
	ImageType::RegionType fullRegion(index,size);
	input->SetRequestedRegion(fullRegion);
	input->Update();

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
	WriteITK(input_temp,"input_unmodified.mhd");

	//Creates an image to store the type of voxel for each voxel in the sub image
	VoxelTypeImage::Pointer voxel_type = VoxelTypeImage::New();
    voxel_type->SetRegions(fullRegion);
    voxel_type->Allocate();

	//Create an iterator for the voxel types
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
	WriteITK(temp, "input gradient.mhd");

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
        VoxelType type=SingleMaterialClassification(input_iter.Get()-1024,temp_iter.Get());
		voxel_type_iter.Set(type);		//sets the type in the voxel_type structure
    }

	//outputs the gradient magnitude purely for reference
	std::cerr<<"Gradient max: "<<gradient_max<<std::endl;

	for (temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++temp_iter)
    {
        switch(voxel_type_iter.Get()) {
            case Tissue:
                temp_iter.Set(1);
				break;
			case Stool:
                temp_iter.Set(2);
                break;
            case Unclassified:
                temp_iter.Set(-1);
                break;
            case Air:
                temp_iter.Set(0);
                break;
        }
    }

	WriteITK(temp,"hard threshold classification.mhd");

	//Filling in Tissue Holes
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

    StructuringElementType  structuringElement;
    structuringElement.SetRadius( .5 );  // 1
    structuringElement.CreateStructuringElement();
	ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
	closingFilter->SetKernel(structuringElement);	 
	closingFilter->SetInput(temp);
	closingFilter->Update();
	temp=closingFilter->GetOutput();
	temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	closingFilter.~SmartPointer();

    //Updates Information after morphological operations of tissue in stool
    //Begins update the unclassified in tissue morphological operation
	for (temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
        ++voxel_type_iter, ++temp_iter)
    {
        switch(voxel_type_iter.Get()) {
            case Tissue:
                //tissue into stool
                if (temp_iter.Get()==1) {
                //    voxel_type_iter.Set(Stool);			//Modify filled tissue holes
                }
                temp_iter.Set(1);
				break;
			case Stool:
                temp_iter.Set(-1);
                break;
            case Unclassified:
                temp_iter.Set(0);
                break;
            case Air:
                temp_iter.Set(-1);
                break;
        }
    }

	closingFilter = ClosingFilterType::New();
	closingFilter->SetKernel(structuringElement);
	closingFilter->SetInput(temp);
	closingFilter->Update();
	temp= closingFilter->GetOutput();
	temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	closingFilter.~SmartPointer();

    ByteImageType::Pointer chamfer_colon=AllocateNewByteImage( fullRegion );
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,fullRegion);

    //Updates Information after morphological operations
    for (chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd();  
        ++voxel_type_iter, ++chamfer_colon_iter) 
    {
        switch(voxel_type_iter.Get()) {
            case Unclassified:
			case Stool:
                if (chamfer_colon_iter.Get()==1) {
                //    voxel_type_iter.Set(Tissue);			//unclassified into tissue
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

	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_colon);
	int chamfer_weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();
	chamfer_colon=chamfer_filter->GetOutput();
	chamfer_colon_iter = IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	chamfer_filter.~SmartPointer();


	ImageType::IndexType endIndex = fullRegion.GetIndex();
	ImageType::IndexType startIndex = fullRegion.GetIndex();	
	endIndex[0]+=(fullRegion.GetSize()[0]-1);
	endIndex[1]+=(fullRegion.GetSize()[1]-1);
	endIndex[2]+=(fullRegion.GetSize()[2]-1);


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

	WriteITK(chamfer_colon,"masked colon.mhd");

    int counter=0;							//used to total certain values (description given when set)
    double mean[3]={0, 0, 0};				//stores the mean of the air/tissue/stool classes
    double sum[3]={0, 0, 0};					//stores the total partials of the three classes
    double variance[3]={0, 0, 0};			//stores the variance of the air/tissue/stool classes
	float weight[3]={0,0,0};				//stores the weights of the air/tissue/stool classes

	//Declares a structure to store the air/tissue/stool partial volumes for a voxel in the image
	ImageVectorType::Pointer partialVector2 = ImageVectorType::New();
    partialVector2->SetRegions(fullRegion);
    partialVector2->Allocate();
	partialVector2->SetSpacing(input->GetSpacing());

	//defines an iterator over the partial volume structure
    IteratorImageVectorType partialVector_iter(partialVector2,fullRegion);

	//computes the gradient magnitude for the image and set it to a temporary variable
    gradient_filter = GradientMagnitudeFilterType::New();
    gradient_filter->SetInput(input_temp);
    gradient_filter->Update();
    temp = gradient_filter->GetOutput();					//this variable will be reused later for other purpose, meaning change will be commented
	temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
	gradient_filter.~SmartPointer();	

	std::cerr<<"Allocated All Spaces"<<std::endl;

	//computes the total number of voxels in the sub image
	int total_vertices=0;
	for(partialVector_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !partialVector_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd(); 
        ++partialVector_iter, ++input_iter, ++chamfer_colon_iter, ++voxel_type_iter) 
	{
		total_vertices++;
		VoxelType type=voxel_type_iter.Get();
		float value[3]={.1,.45,.45};		//declares default partials for unclassified class
		switch (type) {
			case Tissue:
				value[0]=0;
				value[1]=1;					//sets the partial to all tissue
				value[2]=0;
				mean[1]+=input_iter.Get();	//updates the tissue mean
				sum[1]++;					//updates the total partial count
				break;
			case Stool:
				value[0]=0;
				value[1]=0;
				value[2]=1;					//sets the partial to all stool
				mean[2]+=input_iter.Get();	//updates the stool mean
				sum[2]++;					//updates the total stool partial
				break;
			case Air:
				value[0]=1;					//sets the partial to all Air
				value[1]=0;
				value[2]=0;
				mean[0]+=input_iter.Get();	//updates the air mean
				sum[0]++;					//updates the total air partial
				break;
			case Unclassified:
				break;
		}   

		CovariantVectorType data(value);	//sets the partial to be stored in the actual structure
		partialVector_iter.Set(data);		//stores the partial in the structure
    }

	std::cerr<<"All partial Initlaized"<<std::endl;

	//computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
	for (int i=0;i<3;i++) { mean[i]=mean[i]/sum[i]; }    
 
    for(input_iter.GoToBegin(), voxel_type_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
		!input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++input_iter, ++voxel_type_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	{
		if (chamfer_colon_iter.Get()==1) {
			CovariantVectorType data = partialVector_iter.Get();		//gets the partials
			float Z[3]={data[0],data[1],data[2]};

			for (int i=0;i<3;i++) {
				variance[i]+=Z[i]*vnl_math_sqr(input_iter.Get()-mean[i]);	//computes the variance sum(partial[i]*(value[i]-mean[i])^2)/sum(partial[i])
			}
		}
    }

	//computes the weights and varaince
	for (int i=0;i<3;i++) {
		variance[i]=variance[i]/sum[i];
		weight[i]=sum[i]/total_vertices;
	}

	//outputs the initial mean and variance of each class
	std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;


	//declares an image used to store the output image at the end of the i-th iteration
	ImageType::Pointer em_output_image = ImageType::New();
	em_output_image->SetRegions(fullRegion);
	em_output_image->Allocate();
	em_output_image->SetSpacing(input->GetSpacing());
	IteratorTypeFloat4WithIndex em_output_image_iter(em_output_image,fullRegion);

	//used 20 in original test case
	//computes 5 iterations of the Expectation Maximization Algorithm
    for(int i=0;i<5;i++) {
		//gives the temporary storage of the variables corresponding to sum, variance, and mean for each class on the i+1 iteration
		double sum_temp[3]={0,0,0};
        double variance_temp[3]={0,0,0};
		double mean_temp[3]={0,0,0};

		for(input_iter.GoToBegin(), voxel_type_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), temp_iter.GoToBegin(); 
			!input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !temp_iter.IsAtEnd(); 
			++input_iter, ++partialVector_iter, ++voxel_type_iter, ++chamfer_colon_iter, ++temp_iter) 
		{
			CovariantVectorType data = partialVector_iter.Get();		//retrieves the partial informations
			float Z[3]={data[0],data[1],data[2]};

			vnl_vector<float> Z_update=expectation(input_iter.Get(),mean, variance, weight, GetNeighbor(partialVector2,input_iter.GetIndex()), Z, temp_iter.Get());
			data[0]=Z_update[0];										//updates the partial values
			data[1]=Z_update[1];
			data[2]=Z_update[2];
			partialVector_iter.Set(data);

			//updates the new mean, variance, and total partial sum for each class accordingly
			for (int k=0;k<3;k++) {
				mean_temp[k]+=data[k]*input_iter.Get();	
				variance_temp[k]+=data[k]*vnl_math_sqr(input_iter.Get()-mean[k]);
				sum_temp[k]+=data[k];
			}
        }

		//updates the mean, variance, and weights for each class
		for (int k=0;k<3;k++) {
			mean[k]=mean_temp[k]/sum_temp[k];
			variance[k]=variance_temp[k]/sum_temp[k];
			weight[k]=sum_temp[k]/total_vertices;
		}
		std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
		std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
		std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;
    }

	//set modification images type.
	//Updates the existing voxel type classifications with knowledge found in the EM
	EMClassification(input_iter, voxel_type_iter, partialVector_iter, chamfer_colon_iter, temp_iter);
	WriteITK(temp,"Global EM classification.mhd");

	////std::cerr<<"Starting Localized EM"<<std::endl;

	////Declares a neighborhood to grad local neighborhood from with a 1 voxel by 1 voxel by 1 voxel cube
	//NeighborhoodVectorIteratorType::SizeType neighborhood_size={1,1,1};
	//NeighborhoodVectorIteratorType partial_neighborhood(neighborhood_size, partialVector2, fullRegion);
	//NeighborhoodIteratorType input_neighbor(neighborhood_size, input_temp, fullRegion);

	//
	////for(em_output_image_iter.GoToBegin();!em_output_image_iter.IsAtEnd(); ++em_output_image_iter) {
	////	em_output_image_iter.Set(0);
	////}
	//
	////WriteITK(input_temp, "Localized input input.mhd");
	//
	////iterates through the partial image, along with the voxel types
	//for(partialVector_iter.GoToBegin(), voxel_type_iter.GoToBegin()/*, em_output_image_iter.GoToBegin()*/; 
	//	!partialVector_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd()/* && !em_output_image_iter.IsAtEnd()*/; 
	//	++partialVector_iter, ++voxel_type_iter/*, ++em_output_image_iter*/) 
	//{
	//	if (voxel_type_iter.Get()==Unclassified) {					//if the voxel type is unclassified do a local partial
	//		partial_neighborhood.SetLocation(partialVector_iter.GetIndex());	//declares an image region iterator for the neighborhood
	//		input_neighbor.SetLocation(partialVector_iter.GetIndex());
	//
	//		if (partial_neighborhood.InBounds()) {
	//			float local_vertices = partial_neighborhood.GetNeighborhood().Size();
	//			//ImageType::RegionType local_region = partial_neighborhood.GetBoundingBoxAsImageRegion();

	//			//ImageVectorType::Pointer local_partialVector = ImageVectorType::New();
	//			//local_partialVector->SetRegions(local_region);
	//			//local_partialVector->Allocate();
	//			//local_partialVector->SetSpacing(input->GetSpacing());

	//			////defines an iterator over the partial volume structure
	//			//IteratorImageVectorType local_partialVector_iter(local_partialVector,local_region);
	//			//for (local_partialVector_iter.GoToBegin(); !local_partialVector_iter.IsAtEnd(); ++local_partialVector_iter) {
	//			//	float value[2]={.5,.5};				//declares default partials for unclassified class
	//			//	CovariantVectorType data(value);	//sets the partial to be stored in the actual structure
	//			//	local_partialVector_iter.Set(data);
	//			//}

	//			//float local_variance[2]={0,0};
	//			//float local_mean[2]={0,0};
	//			//float local_weight[2]={.5,.5};
	//			//
	//			//for (int j=0; j<local_vertices; j++) {
	//			//	local_mean[0]+=input_neighbor.GetPixel(j);
	//			//}
	//			//float cut_off=local_mean[0]/local_vertices;
	//			//for (int j=0; j<local_vertices; j++) {
	//			//	if(input_neighbor.GetPixel(j)>cut_off) {
	//			//		local_mean[1]+=input_neighbor.GetPixel(j);
	//			//		++local_weight[1];
	//			//	} else {
	//			//		local_mean[0]+=input_neighbor.GetPixel(j);
	//			//		++local_weight[0];
	//			//	}
	//			//}
	//			//local_mean[0]/=local_weight[0];
	//			//local_mean[1]/=local_weight[1];
	//			//for (int j=0; j<local_vertices; j++) {
	//			//	if(input_neighbor.GetPixel(j)>cut_off) {
	//			//		local_variance[1]+=vnl_math_sqr(local_mean[1]-input_neighbor.GetPixel(j));
	//			//	} else {
	//			//		local_variance[0]+=vnl_math_sqr(local_mean[0]-input_neighbor.GetPixel(j));
	//			//	}
	//			//}
	//			//local_variance[0]/=local_weight[0];
	//			//local_variance[1]/=local_weight[1];
	//			//local_weight[0]/=local_vertices;
	//			//local_weight[1]/=local_vertices;

	//			//if (!(local_mean[0]>0)) {
	//			//	std::cerr<<local_mean[0]<<" "<< local_mean[1]<<std::endl;
	//			//	std::cerr<<local_variance[0]<<" "<< local_variance[1]<<std::endl;
	//			//	std::cerr<<local_weight[0]<<" "<< local_weight[1]<<std::endl;
	//			//}

	//			////std::cerr<<cut_off<<std::endl;
	//			//CovariantVectorType local_center_partial;

	//			//for (int i=1;i<=20;i++) {
	//			////	//declares the temp sum/weight/variance for the iteration
	//			//	float sum_temp[2]={0,0};
	//			//	float variance_temp[2]={0,0};
	//			//	float mean_temp[2]={0,0};

	//			//	int j=0;
	//			//	for (j=0, local_partialVector_iter.GoToBegin(); j<local_vertices && !local_partialVector_iter.IsAtEnd(); j++, ++local_partialVector_iter) {
	//			//		CovariantVectorType data = local_partialVector_iter.Get();	//get the data
	//			//		float Z[2]={data[0],data[1]};
	//			//		if (Z[0]==1 || Z[1]==1) {
	//			//			break;
	//			//		}

	//			////		//this is exactly like the above expectation except no local weighting as it is already a local neighborhood
	//			////		//computes the updated partials
	//			//		vnl_vector<float> Z_update=expectation(input_neighbor.GetPixel(j), local_mean, local_variance, local_weight, Z);	
	//			////		//std::cerr<<Z_update[0]<<" "<<Z_update[1]<<" "<<Z_update[2]<<std::endl;
	//			////		//updates for the i+1 th iteration's mean, variance, and sum
	//			//		for (int k=0;k<2;k++) {
	//			//			mean_temp[k]+=Z_update[k]*input_neighbor.GetPixel(j);
	//			//			variance_temp[k]+=Z_update[k]*vnl_math_sqr(input_neighbor.GetPixel(j)-local_mean[k]);
	//			//			sum_temp[k]+=Z_update[k];
	//			//		}
	//			//		data[0]=Z_update[0];										//updates the partial values
	//			//		data[1]=Z_update[1];
	//			//		local_partialVector_iter.Set(data);

	//			//		//if (!(data[0]>=0 && data[0]<=1) || !(data[1]>=0 && data[1]<=1)) {
	//			//		//	std::cerr<<input_neighbor.GetPixel(j)<<" "<<Z[0]<<" "<<Z[1]<<std::endl;
	//			//		//	std::cerr<<local_mean[0]<<" "<< local_mean[1]<<std::endl;
	//			//		//	std::cerr<<local_variance[0]<<" "<< local_variance[1]<<std::endl;
	//			//		//	std::cerr<<local_weight[0]<<" "<< local_weight[1]<<std::endl;
	//			//		//}

	//			//		if (j==ceil(local_vertices/2)) {
	//			//			local_center_partial = data;
	//			//		}
	//			//	}
	//			//	//updates the new local partials
	//			//	for (int k=0;k<2;k++) {
	//			//		local_mean[k]=mean_temp[k]/sum_temp[k];
	//			//		local_variance[k]=variance_temp[k]/sum_temp[k];
	//			//		local_weight[k]=sum_temp[k]/local_vertices;
	//			//	}
	//			//	//std::cerr<<"Iteration: "<<i<<std::endl;
	//			//	//std::cerr<<"Mean: "<<local_mean[0]<<" "<<local_mean[1]<<std::endl;
	//			//	//std::cerr<<"Variance: "<<local_variance[0]<<" "<<local_variance[1]<<std::endl;
	//			//	//std::cerr<<"Weight: "<<local_weight[0]<<" "<<local_weight[1]<<std::endl;
	//			//}

	//			//CovariantVectorType data = partialVector_iter.Get();	//get the data
	//			//if (local_mean[0]<500) {
	//			//	data[0]=local_center_partial[0];
	//			//	if (local_variance[1]>100000 || (local_mean[1]-1024)>800) {
	//			//		data[1]=0;
	//			//		data[2]=local_center_partial[1];
	//			//	} else {
	//			//		data[1]=local_center_partial[1];
	//			//		data[2]=0;			
	//			//	}
	//			//} else {
	//			//	data[0]=0;
	//			//	data[1]=local_center_partial[0];
	//			//	data[2]=local_center_partial[1];
	//			//}
	//			//partialVector_iter.Set(data);
	//			//local_partialVector.~SmartPointer();

	//			//declares local neighborhood variables for the sum/variance/mean of each partial class
	//			float local_variance[3]={variance[0],variance[1],variance[2]};
	//			float local_mean[3]={mean[0],mean[1],mean[2]};
	//			float local_weight[3]={.3,.3,.3};

	//			//outputs the mean and variance of each class on the (i+1)-th iteration
	//			//std::cerr<<"Start (total verticies) "<<local_vertices<<std::endl;
	//			//std::cerr<<"Mean: "<<local_mean[0]<<" "<<local_mean[1]<<" "<<local_mean[2]<<std::endl;
	//			//std::cerr<<"Variance: "<<local_variance[0]<<" "<<local_variance[1]<<" "<<local_variance[2]<<std::endl;
	//			//std::cerr<<"Weight: "<<local_weight[0]<<" "<<local_weight[1]<<" "<<local_weight[2]<<std::endl;

	//			for (int i=1;i<=20;i++) {
	//				//declares the temp sum/weight/variance for the iteration
	//				float sum_temp[3]={0,0,0};
	//				float variance_temp[3]={0,0,0};
	//				float mean_temp[3]={0,0,0};
	//		
	//				for (int j=0;j<local_vertices;j++) {
	//					CovariantVectorType data = partial_neighborhood.GetPixel(j);	//get the data
	//					float Z[3]={data[0],data[1],data[2]};
	//					if (Z[0]==1 || Z[1]==1 || Z[2]==1) {
	//						break;
	//					}
	//					//this is exactly like the above expectation except no local weighting as it is already a local neighborhood
	//					//computes the updated partials
	//					vnl_vector<float> Z_update=expectation(input_neighbor.GetPixel(j), local_mean, local_variance, local_weight, Z);	
	//					//std::cerr<<Z_update[0]<<" "<<Z_update[1]<<" "<<Z_update[2]<<std::endl;
	//					//updates for the i+1 th iteration's mean, variance, and sum
	//					for (int k=0;k<3;k++) {
	//						if (local_weight[k]!=0) {
	//							mean_temp[k]+=Z_update[k]*input_neighbor.GetPixel(j);
	//							variance_temp[k]+=Z_update[k]*vnl_math_sqr(input_neighbor.GetPixel(j)-local_mean[k]);
	//							sum_temp[k]+=Z_update[k];
	//						}
	//					}

	//					if (j==ceil(local_vertices/2)) {
	//						data[0]=Z_update[0];			//updates the result
	//						data[1]=Z_update[1];
	//						data[2]=Z_update[2];
	//						partial_neighborhood.SetPixel(j, data);
	//					}
	//				}

	//				//updates the new local partials
	//				for (int k=0;k<3;k++) {
	//					if (sum_temp[k]>0 && local_weight[k]!=0) {
	//						local_mean[k]=mean_temp[k]/sum_temp[k];
	//						local_variance[k]=variance_temp[k]/sum_temp[k];
	//						local_weight[k]=sum_temp[k]/local_vertices;
	//					} else {
	//						local_weight[k]=0;
	//					}
	//				}
	//			}
	//			//std::cerr<<"Iteration: "<<i<<std::endl;
	//			//std::cerr<<"Mean: "<<local_mean[0]<<" "<<local_mean[1]<<" "<<local_mean[2]<<std::endl;
	//			//std::cerr<<"Variance: "<<local_variance[0]<<" "<<local_variance[1]<<" "<<local_variance[2]<<std::endl;
	//			//std::cerr<<"Weight: "<<local_weight[0]<<" "<<local_weight[1]<<" "<<local_weight[2]<<std::endl;
	//		}
	//	}
	//}
	////WriteITK(em_output_image, "Localized input.mhd");

	for(input_iter.GoToBegin(), em_output_image_iter.GoToBegin(), partialVector_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
		!input_iter.IsAtEnd() && !em_output_image_iter.IsAtEnd() && !partialVector_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++input_iter, ++em_output_image_iter, ++partialVector_iter, ++chamfer_colon_iter) 
	{
		CovariantVectorType data = partialVector_iter.Get();		//retrieves the partial informations
		em_output_image_iter.Set(data[1]);		//sets the voxel intensity for the i-th iteration output image
	}

	WriteITK(em_output_image, "tissue partial.mhd");

	//em_output_image.~SmartPointer();
	std::cerr<<"Finished EM"<<std::endl;
	em_output_image.~SmartPointer();
	partialVector2.~SmartPointer();
	std::cerr<<"partialVector2 has been destroyed: "<<partialVector2.IsNull()<<std::endl;


    //get cubic bspline image interpolation
    InterpolationType::Pointer input_interpolator = InterpolationType::New();
    input_interpolator->SetSplineOrder(3);
    input_interpolator->SetInputImage(input);
	//Storage for Smax
    ImageType::Pointer input_smax = AllocateNewImage(fullRegion);
    IteratorTypeFloat4WithIndex input_smax_iter(input_smax,fullRegion);
	//Debugging Flags
	bool StoolAirbool=false;
	bool TissueStoolbool=false;
	bool TissueAirbool = false;


	//Distribution Code Here

	//ofstream distribution_record;
	//distribution_record.open("C:/ds/distribution_record.csv");
	//for (voxel_type_iter.GoToBegin(), input_iter.GoToBegin(); 
	//	!voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd(); 
	//	++voxel_type_iter, ++input_iter) {
	//		switch (voxel_type_iter.Get()) {
	//			case Tissue:
	//				distribution_record<<"1";
	//				break;
	//			case Stool:
	//				distribution_record<<"2";
	//				break;
	//			case Air:
	//				distribution_record<<"0";
	//				break;
	//			default:
	//				distribution_record<<"3";
	//				break;
	//		}
	//		distribution_record<<","<<input_iter.Get()<<endl;
	//}
	//distribution_record.close();

	//Runs edge classification with 3 different increments
    for (voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), temp_iter.GoToBegin(); 
		!voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !temp_iter.IsAtEnd(); 
		++voxel_type_iter, ++input_smax_iter, ++temp_iter) {
		if (voxel_type_iter.Get()==Unclassified) {
			//Initialization
			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
			float temp_threshold=-1;
			VoxelType temp_type = Unclassified;
			//+/- 1.5, +/-1.0
			VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,
				input_interpolator,input_smax_iter,voxel_index);
			//+/- 1.0, +/-0.5
			VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,
				input_interpolator,input_smax_iter,voxel_index);
			//+/- .6, +/-.3
			VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,
				input_interpolator,input_smax_iter,voxel_index);
			//std::cerr<<"result: "<<temp_type<<" "<<temp_threshold<<std::endl;
			if (StoolAirbool==false) {
				StoolAirbool=true;
				//std::cerr<<"Stool Air Region Detected"<<std::endl;
			}
			if (TissueStoolbool==false) {
				TissueStoolbool=true;
				//std::cerr<<"Tissue Stool Region Detected"<<std::endl;
			}

			//if (TissueAirbool==false) {
			//	TissueAirbool=true;
			//	std::cerr<<"Tissue Air Region Detected"<<std::endl;
			//}

			voxel_type_iter.Set(temp_type);
		}

		

		//Set up data for heterogenous stool

        switch(voxel_type_iter.Get()) {
			case TissueStool:
                temp_iter.Set(0);
                break;
            case Stool:
                temp_iter.Set(1);
                break;
			case Tissue:
				temp_iter.Set(3);
				break;
			case Air:
				temp_iter.Set(-2);
				break;
			case TissueAir:
				temp_iter.Set(2);
				break;
			case StoolAir:
				temp_iter.Set(-1);
				break;
            default:
                temp_iter.Set(-3);
                break;
        }
    }

	std::cerr<<"Done Edge"<<std::endl;
	//Write out unmodified image
	WriteITK(temp,"tissue stool.mhd");
	WriteITK(input_smax,"smax.mhd");


	//Deletes the Interpolators
	input_interpolator.~SmartPointer();
	
	//std::cerr<<"set partial volume "<<std::endl;
    //0: ps
    //1: pt
    //2: pa
    ImageType::Pointer partialVolume = AllocateNewImage(fullRegion);
    IteratorTypeFloat4WithIndex partialVolume_iter(partialVolume,fullRegion);
	//bool random=true;

	ByteImageType::Pointer chamfer_from_stool_involvement=AllocateNewByteImage(fullRegion);
	IteratorTypeByteWithIndex chamfer_from_stool_involvement_iter(chamfer_from_stool_involvement,fullRegion);

	//Sets the partial volumes
    for(voxel_type_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin(), input_smax_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
        ++voxel_type_iter, ++input_iter, ++partialVolume_iter, ++input_smax_iter, ++chamfer_from_stool_involvement_iter) 
    {
		//prepares information for chamfer stool calculation
		if (voxel_type_iter.Get()==Stool || voxel_type_iter.Get()==TissueStool || voxel_type_iter.Get()==StoolAir) {
			chamfer_from_stool_involvement_iter.Set(1);     //stool
		} else {
			chamfer_from_stool_involvement_iter.Set(0);     //non stool is the object
		}
		partialVolume_iter.Set(0);
		float value=0;
        switch(voxel_type_iter.Get()) {
			case Tissue:
				partialVolume_iter.Set(1.0);
				break;
			case Stool:
			case Air:
				partialVolume_iter.Set(0.0);
				break;
            case TissueAir:
				//std::cerr<<"TA reached"<<std::endl;
				
				partialVolume_iter.Set(1.0);	
				/*value = 1.0+(input_iter.Get()-1024)/1000;
				if (value<=1 && value>=0) {
					partialVolume_iter.Set(value);	
				} else if (value<0) {
					partialVolume_iter.Set(0);	
				} else {
					partialVolume_iter.Set(1);	
				}*/
                break;
            case TissueStool:
					value=1.0-(input_iter.Get()-1024)/(input_smax_iter.Get());
					if (value<=1 && value>=0) {
						partialVolume_iter.Set(value);
					} else if (value<0) {
						partialVolume_iter.Set(0);	
					} else {
						partialVolume_iter.Set(1);	
					}

					////special case when there is nothing really to fit
					//if (input_smax_iter.Get()-input_iter.Get()<50 && input_iter.Get()-1024<500) {
					//	partialVolume_iter.Set(1);	
					//}
                break;
            case StoolAir:
                partialVolume_iter.Set(0);
                break;
        }
		if (-partialVolume_iter.Get()==std::numeric_limits<float>::infinity()) {
			//std::cerr<<partialVolume_iter.GetIndex()[0]<<" "<<partialVolume_iter.GetIndex()[1]<<" "<<partialVolume_iter.GetIndex()[2]<<std::endl;
		}
    }

	input_smax.~SmartPointer();

	//Write out partial volume
	//WriteITK(partialVolume,"tissue partial one.mhd");

	//Write out unmodified image
	//WriteITK(temp,"chamfer from stool input.mhd");

	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_from_stool_involvement);
	int weights[3]={3,4,5};	//3d distance weight recommended by julian
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(false);
	chamfer_filter->Update();
	chamfer_from_stool_involvement=chamfer_filter->GetOutput();
	chamfer_from_stool_involvement_iter = IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	chamfer_filter.~SmartPointer();


	//Write out chamfer thickness of stool
	//WriteITK(chamfer_from_stool_involvement,"chamfer from stool.mhd");
	//Updates chamfer by setting all stool attached to an 8+ chamfer stool to also 8 chamfer units
	int chamfercutoff=8;
	for (int i=0;i<2;i++) {
		//Forward
		for(chamfer_from_stool_involvement_iter.GoToBegin(); !chamfer_from_stool_involvement_iter.IsAtEnd(); ++chamfer_from_stool_involvement_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//Write the modified chamfer map
		//WriteITK(chamfer_from_stool_involvement,"chamfer from stool forward.mhd");
		
		//Backward
		for(chamfer_from_stool_involvement_iter.GoToReverseBegin(); !chamfer_from_stool_involvement_iter.IsAtReverseEnd(); --chamfer_from_stool_involvement_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//write the modified chamfer map
		//WriteITK(chamfer_from_stool_involvement,"chamfer from stool back.mhd");
	}


    ByteImageType::Pointer chamfer_air=AllocateNewByteImage( fullRegion );
    IteratorTypeByteWithIndex chamfer_air_iter(chamfer_air,fullRegion);

    for(voxel_type_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd();
        ++voxel_type_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter, ++input_iter, ++partialVolume_iter)
    {
		//prepare for computation of chamfer to air
		if (voxel_type_iter.Get()==Air) {
            chamfer_air_iter.Set(0);					//distance to AIR
        } else {
            chamfer_air_iter.Set(1);					//non important locations
        }
		float temp_value=chamfer_from_stool_involvement_iter.Get();
		//the temp_value>3 is uncessary
		if ((temp_value>=chamfercutoff && temp_value>3 && voxel_type_iter.Get()!=TissueStool) || voxel_type_iter.Get()==Stool) {
			partialVolume_iter.Set(0);
		} else if (temp_value<chamfercutoff && temp_value>3) {
			voxel_type_iter.Set(ThinStool);		//if the chamfer < 8 then it is set to thin stool
		}
    }

	//compute chamfer for air
	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_air);
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(false);
	chamfer_filter->Update();
	chamfer_air=chamfer_filter->GetOutput();
	chamfer_air_iter= IteratorTypeByteWithIndex(chamfer_air,fullRegion);
	chamfer_filter.~SmartPointer();
	//Writes out the chamfer to air
	//WriteITK(chamfer_air,"chamfer from air.mhd");
	//Normalize the data

    for(chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
        !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
        ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter)
    {
        chamfer_from_stool_involvement_iter.Set(chamfer_from_stool_involvement_iter.Get()/chamfercutoff);
        if (chamfer_air_iter.Get()>=chamfercutoff) {
            chamfer_air_iter.Set(1);
        } else {
            chamfer_air_iter.Set(chamfer_air_iter.Get()/8);
        }
    }

    //Calculating partial air
    for(voxel_type_iter.GoToBegin(), partialVolume_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
        !voxel_type_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
		++voxel_type_iter, ++partialVolume_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter) 
	{
        if(voxel_type_iter.Get()==ThinStool) {
			float value=0;
			float ps=1/2*(1+vnl_erf((chamfer_from_stool_involvement_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
			float pa=1/2*(1+vnl_erf((chamfer_air_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
			value=1-ps-pa;
			partialVolume_iter.Set(value);
        }
	}

	for(partialVolume_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin(); 
		!partialVolume_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd(); 
		++partialVolume_iter, ++chamfer_colon_iter, ++voxel_type_iter) 
	{
		if (chamfer_colon_iter.Get()==1) { 
			if (partialVolume_iter.Get()<.3) {
				partialVolume_iter.Set(0);
			}
			if (voxel_type_iter.Get() == TissueAir) {
				partialVolume_iter.Set(1);
			}
		}
	}

	//Write out the partial tissue with thinstool modification done
	WriteITK(partialVolume,"tissue partial two.mhd");

	//Prepares the output image
	ImageType::Pointer output = AllocateNewImage(fullRegion);
	output->SetSpacing(input->GetSpacing());
	IteratorTypeFloat4WithIndex output_iter(output,fullRegion);
    //Modifies to intensity for output
    for(output_iter.GoToBegin(), input_iter.GoToBegin(), temp_iter.GoToBegin(), voxel_type_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), partialVolume_iter.GoToBegin();
		!output_iter.IsAtEnd() && !input_iter.IsAtEnd() && !temp_iter.IsAtEnd() &&  !voxel_type_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd();
		++output_iter, ++input_iter, ++temp_iter, ++voxel_type_iter, ++chamfer_colon_iter, ++partialVolume_iter) 
    {	
		if (voxel_type_iter.Get()==Air|| voxel_type_iter.Get()==StoolAir) {
			output_iter.Set(mean[0]);
		} else if (chamfer_colon_iter.Get()==1 || voxel_type_iter.Get()==TissueStool  || voxel_type_iter.Get()==Stool) { 
	 //       ImageType::IndexType index =input_iter.GetIndex();
			output_iter.Set(partialVolume_iter.Get()*input_iter.Get());
		} else {
			output_iter.Set(input_iter.Get());
		}
	}
	//voxel_type.~SmartPointer();
	chamfer_air.~SmartPointer();
	//Writes out the output preclosing
	WriteITK(output,"output one.mhd");
	for(output_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !output_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd();
		++output_iter, ++chamfer_from_stool_involvement_iter, ++voxel_type_iter) 
	{
		if (output_iter.Get()>300 || voxel_type_iter.Get()==Tissue || voxel_type_iter.Get()==TissueAir) {
			chamfer_from_stool_involvement_iter.Set(1);
		} else {
			chamfer_from_stool_involvement_iter.Set(0);
		}
	}
	//Write out unmodified image
	WriteITK(temp,"chamfer from tissue input.mhd");
	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_from_stool_involvement);
	chamfer_filter->SetWeights(weights, weights+3);
	chamfer_filter->SetDistanceFromObject(false);
	chamfer_filter->Update();

	chamfer_from_stool_involvement = chamfer_filter->GetOutput();
	chamfer_from_stool_involvement_iter=IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
	chamfer_filter.~SmartPointer();

	//Write out chamfer thickness of stool
	WriteITK(chamfer_from_stool_involvement,"chamfer tissue stool.mhd");

	
	ImageType::Pointer min = longruntests(output, 2, 0, input->GetSpacing().GetElement(2));
	IteratorTypeFloat4WithIndex min_iter(min,min->GetLargestPossibleRegion());

	float modifiers[9] = { input->GetSpacing().GetElement(2), //requires z axis 
		input->GetSpacing().GetElement(2),	//requires z axis
		input->GetSpacing().GetElement(0),

		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1)),
		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1)),

		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1)),
		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1)),

		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1)),
		sqrt(input->GetSpacing().GetElement(2)*input->GetSpacing().GetElement(2)+input->GetSpacing().GetElement(1)*input->GetSpacing().GetElement(1))};

	std::cerr<<modifiers[0]<<" "<<modifiers[1]<<" "<<modifiers[2]<<" "<<modifiers[3]<<std::endl;
	for (int i=1;i<9;i++) {
		ImageType::Pointer temp_image = longruntests(output, 2, i, 1/modifiers[i]);
		
		IteratorTypeFloat4WithIndex axial_iter=IteratorTypeFloat4WithIndex(temp_image,temp_image->GetLargestPossibleRegion());
		for (axial_iter.GoToBegin(), min_iter.GoToBegin();
			!axial_iter.IsAtEnd()&& !min_iter.IsAtEnd();
			++axial_iter, ++min_iter) 
		{
			if (axial_iter.Get()<min_iter.Get()) {
				min_iter.Set(axial_iter.Get());
			}
		}
		temp_image.~SmartPointer();
	}


	WriteITK(min, "min variance.mhd");

	chamfercutoff=20;
	for(chamfer_from_stool_involvement_iter.GoToBegin(); !chamfer_from_stool_involvement_iter.IsAtEnd(); ++chamfer_from_stool_involvement_iter) {
		if (chamfer_from_stool_involvement_iter.Get()>chamfercutoff) {
			chamfer_from_stool_involvement_iter.Set(chamfercutoff);
		}
	}

	//Updates chamfer by setting all stool attached to an 20+ chamfer stool to also 8 chamfer units
	
	for (int i=0;i<2;i++) {
		//Forward
		for(chamfer_from_stool_involvement_iter.GoToBegin(), min_iter.GoToBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtEnd() && !min_iter.IsAtEnd(); 
			++chamfer_from_stool_involvement_iter, ++min_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff && min_iter.Get()<500 /*|| temp_value>chamfercutoff*/) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if ((chamfer_from_stool_involvement->GetPixel(temp_index)>0)) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//Write the modified chamfer map
		//WriteITK(chamfer_from_stool_involvement,"chamfer from tissue forward.mhd");
		
		//Backward
		for(chamfer_from_stool_involvement_iter.GoToReverseBegin(), min_iter.GoToReverseBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtReverseEnd() && !min_iter.IsAtReverseEnd(); 
			--min_iter, --chamfer_from_stool_involvement_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff && min_iter.Get()<500 /*|| temp_value>chamfercutoff*/) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if ((chamfer_from_stool_involvement->GetPixel(temp_index)>0)) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//write the modified chamfer map
		WriteITK(chamfer_from_stool_involvement,"chamfer from tissue back.mhd");
	}

	temp.~SmartPointer();


	ImageType::Pointer curvesFilter_image = levelset_comp(output);
	WriteITK(curvesFilter_image, "levelset_output.mhd");
	std::cerr<<curvesFilter_image->GetLargestPossibleRegion().GetIndex()[0]<<" "<<curvesFilter_image->GetLargestPossibleRegion().GetIndex()[1]<<" "<<curvesFilter_image->GetLargestPossibleRegion().GetIndex()[2]<<std::endl;
	std::cerr<<output->GetLargestPossibleRegion().GetIndex()[0]<<" "<<output->GetLargestPossibleRegion().GetIndex()[1]<<" "<<output->GetLargestPossibleRegion().GetIndex()[2]<<std::endl;
	std::cerr<<chamfer_from_stool_involvement->GetLargestPossibleRegion().GetIndex()[0]<<" "<<chamfer_from_stool_involvement->GetLargestPossibleRegion().GetIndex()[1]<<" "<<chamfer_from_stool_involvement->GetLargestPossibleRegion().GetIndex()[2]<<std::endl;
	IteratorTypeFloat4WithIndex curvesFilter_iter(curvesFilter_image,curvesFilter_image->GetLargestPossibleRegion());


	for (int i=0;i<2;i++) {
		//Forward
		for(chamfer_from_stool_involvement_iter.GoToBegin(), min_iter.GoToBegin(), curvesFilter_iter.GoToBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtEnd() && !min_iter.IsAtEnd() && !curvesFilter_iter.IsAtEnd(); 
			++chamfer_from_stool_involvement_iter, ++min_iter, ++curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff-2/* && curvesFilter_iter.Get()>-.5 && curvesFilter_iter.Get()<.5)*/) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if ((chamfer_from_stool_involvement->GetPixel(temp_index)>0 && chamfer_from_stool_involvement->GetPixel(temp_index)<temp_value) && min->GetPixel(temp_index)<10000) {
											float pixel_value = curvesFilter_image->GetPixel(temp_index);
											//if ((pixel_value-curvesFilter_iter.Get()>-.5 && curvesFilter_iter.Get()> -1 && pixel_value < 1) ||
											//	(pixel_value>-2 && curvesFilter_iter.Get()<=-1 && pixel_value<-1)) {
											//	chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff-1);
											//}
											if (pixel_value-curvesFilter_iter.Get()>-.5 && curvesFilter_image->GetPixel(temp_index)> -1 && pixel_value < 1 ) {
												chamfer_from_stool_involvement->SetPixel(temp_index, chamfercutoff - 2);
											}

											if (pixel_value>-2 && curvesFilter_iter.Get()<=-1 && chamfer_from_stool_involvement_iter.Get()>=chamfercutoff-1) {
												chamfer_from_stool_involvement->SetPixel(temp_index, chamfercutoff-1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//Write the modified chamfer map
		//WriteITK(chamfer_from_stool_involvement,"chamfer from tissue forward 2.mhd");
		
		//Backward
		for(chamfer_from_stool_involvement_iter.GoToReverseBegin(), min_iter.GoToReverseBegin(), curvesFilter_iter.GoToReverseBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtReverseEnd() && !min_iter.IsAtReverseEnd() && !curvesFilter_iter.IsAtReverseEnd(); 
			--min_iter, --chamfer_from_stool_involvement_iter, --curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value>=chamfercutoff-2) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if ((chamfer_from_stool_involvement->GetPixel(temp_index)>0) && chamfer_from_stool_involvement->GetPixel(temp_index)<temp_value && min->GetPixel(temp_index)<10000) {
											float pixel_value = curvesFilter_image->GetPixel(temp_index);
											//if ((pixel_value-curvesFilter_iter.Get()>-0.5  && curvesFilter_iter.Get()> -1 && pixel_value < 1) ||
											//	(pixel_value>-2 && curvesFilter_iter.Get()<=-1 && pixel_value<-1)) {
											//	chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff-1);
											//}
											if (pixel_value-curvesFilter_iter.Get()>-.5 && curvesFilter_image->GetPixel(temp_index)> -1 && pixel_value < 1 ) {
												chamfer_from_stool_involvement->SetPixel(temp_index, chamfercutoff - 2);
											}

											if (pixel_value>-2 && curvesFilter_iter.Get()<=-1 && chamfer_from_stool_involvement_iter.Get()>=chamfercutoff-1) {
												chamfer_from_stool_involvement->SetPixel(temp_index, chamfercutoff - 1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//write the modified chamfer map
		WriteITK(chamfer_from_stool_involvement,"chamfer from tissue back 2.mhd");
	}
	
	WriteITK(chamfer_colon,"chamfer_colon.mhd");


	ByteImageType::Pointer chamfer_heterogenous = AllocateNewByteImage(fullRegion);
	IteratorTypeByteWithIndex chamfer_heterogenous_iter(chamfer_heterogenous,fullRegion);
	for(chamfer_heterogenous_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), curvesFilter_iter.GoToBegin(), min_iter.GoToBegin(); 
		!chamfer_heterogenous_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !curvesFilter_iter.IsAtEnd() && !min_iter.IsAtEnd(); 
		++chamfer_heterogenous_iter, ++chamfer_from_stool_involvement_iter, ++chamfer_colon_iter, ++curvesFilter_iter, ++min_iter){
		if ((chamfer_from_stool_involvement_iter.Get()<chamfercutoff-2 && 
			(chamfer_from_stool_involvement_iter.Get()>=10 || ( chamfer_from_stool_involvement_iter.Get()>=7 && (curvesFilter_iter.Get()<=-4))))  
			/*|| (chamfer_from_stool_involvement_iter.Get()==19 && chamfer_from_stool_involvement_iter.Get()>5 && chamfer_colon_iter.Get()==0)*/) {
			chamfer_colon_iter.Set(1);
			chamfer_heterogenous_iter.Set(1);
		}
	}
		
	WriteITK(chamfer_heterogenous,"chamfer partial.mhd");

	chamfer_filter = ChamferDistanceFilterType::New();
	chamfer_filter->SetInput(chamfer_heterogenous);
	chamfer_filter->SetWeights(chamfer_weights, chamfer_weights+3);
	chamfer_filter->SetDistanceFromObject(true);
	chamfer_filter->Update();

    chamfer_heterogenous=chamfer_filter->GetOutput();
	chamfer_heterogenous_iter=IteratorTypeByteWithIndex(chamfer_heterogenous,fullRegion);
	chamfer_filter.~SmartPointer();

	for(chamfer_heterogenous_iter.GoToBegin(); 
		!chamfer_heterogenous_iter.IsAtEnd(); 
		++chamfer_heterogenous_iter){
		if (chamfer_heterogenous_iter.Get()>=10) {
			chamfer_heterogenous_iter.Set(10);
		}
	}
	WriteITK(chamfer_heterogenous,"chamfer hetero.mhd");
	WriteITK(chamfer_colon,"chamfer_colon_1.mhd");

	int heterogeneous_cuttoff = 5;
	for (int i=0;i<2;i++) {
		//Forward
		for(chamfer_from_stool_involvement_iter.GoToBegin(), chamfer_heterogenous_iter.GoToBegin(), curvesFilter_iter.GoToBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtEnd() && !chamfer_heterogenous_iter.IsAtEnd() && !curvesFilter_iter.IsAtEnd(); 
			++chamfer_from_stool_involvement_iter, ++chamfer_heterogenous_iter, ++curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (chamfer_heterogenous_iter.Get()<=heterogeneous_cuttoff) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};

										float pixel_value = curvesFilter_image->GetPixel(temp_index);
										if (chamfer_from_stool_involvement->GetPixel(temp_index)==19 || chamfer_from_stool_involvement->GetPixel(temp_index)==18) {
											if (pixel_value-curvesFilter_iter.Get()>-.1 && curvesFilter_image->GetPixel(temp_index)> -1) {
												chamfer_heterogenous->SetPixel(temp_index, heterogeneous_cuttoff);
												chamfer_from_stool_involvement->SetPixel(temp_index, 10);
												chamfer_colon->SetPixel(temp_index, 1);
											}
											if (curvesFilter_iter.Get()<0.1 && chamfer_heterogenous_iter.Get()<heterogeneous_cuttoff && min->GetPixel(index)<5000) {
												chamfer_heterogenous->SetPixel(temp_index, heterogeneous_cuttoff-1);
												chamfer_from_stool_involvement->SetPixel(temp_index, 10);
												chamfer_colon->SetPixel(temp_index, 1);
											}

										}
									}
								}
							}
						}
					}
				}
			}
		}

		//Write the modified chamfer map
		WriteITK(chamfer_heterogenous,"chamfer from heterogenous forward.mhd");
		
		//Backward
		for(chamfer_from_stool_involvement_iter.GoToReverseBegin(), chamfer_heterogenous_iter.GoToReverseBegin(), curvesFilter_iter.GoToReverseBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtReverseEnd() && !chamfer_heterogenous_iter.IsAtReverseEnd() && !curvesFilter_iter.IsAtReverseEnd(); 
			--chamfer_heterogenous_iter, --chamfer_from_stool_involvement_iter, --curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (chamfer_heterogenous_iter.Get()<=heterogeneous_cuttoff) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										float pixel_value = curvesFilter_image->GetPixel(temp_index);
										if (chamfer_from_stool_involvement->GetPixel(temp_index)==19 || chamfer_from_stool_involvement->GetPixel(temp_index)==18) {
											if (pixel_value-curvesFilter_iter.Get()>-.1 && curvesFilter_image->GetPixel(temp_index)> -1) {
												chamfer_heterogenous->SetPixel(temp_index, heterogeneous_cuttoff);
												chamfer_from_stool_involvement->SetPixel(temp_index, 10);
												chamfer_colon->SetPixel(temp_index, 1);
											}
											if (curvesFilter_iter.Get()<0.1 && chamfer_heterogenous_iter.Get()<heterogeneous_cuttoff && min->GetPixel(index)<5000) {
												chamfer_heterogenous->SetPixel(temp_index, heterogeneous_cuttoff-1);
												chamfer_from_stool_involvement->SetPixel(temp_index, 10);
												chamfer_colon->SetPixel(temp_index, 1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//write the modified chamfer map
		WriteITK(chamfer_heterogenous,"chamfer from heterogenous back.mhd");
	}
	
	//NOTE: the following two execution seems collapsible into the previous loop.

	for(chamfer_heterogenous_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), partialVolume_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); 
		!chamfer_heterogenous_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd(); 
		++chamfer_heterogenous_iter, ++chamfer_from_stool_involvement_iter, ++partialVolume_iter, ++chamfer_colon_iter)
	{
		if (chamfer_heterogenous_iter.Get()<=heterogeneous_cuttoff) {
			chamfer_from_stool_involvement_iter.Set(0);
			chamfer_colon_iter.Set(1);
			partialVolume_iter.Set(0);
		} else if (chamfer_from_stool_involvement_iter.Get()<chamfercutoff-2) {
			chamfer_from_stool_involvement_iter.Set(0);
			chamfer_colon_iter.Set(1);
		}
	}
	std::cerr<<chamfercutoff<<std::endl;
	WriteITK(chamfer_from_stool_involvement,"chamfer updated.mhd");

	for (int i=0;i<2;i++) {
		//Forward
		for(chamfer_from_stool_involvement_iter.GoToBegin(), min_iter.GoToBegin(), curvesFilter_iter.GoToBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtEnd() && !min_iter.IsAtEnd() && !curvesFilter_iter.IsAtEnd(); 
			++chamfer_from_stool_involvement_iter, ++min_iter, ++curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value<=chamfercutoff/* && curvesFilter_iter.Get()>-.5 && curvesFilter_iter.Get()<.5)*/) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff+chamfer_from_stool_involvement->GetPixel(temp_index));
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//Write the modified chamfer map
		WriteITK(chamfer_from_stool_involvement,"chamfer from tissue forward 3.mhd");
		
		//Backward
		for(chamfer_from_stool_involvement_iter.GoToReverseBegin(), min_iter.GoToReverseBegin(), curvesFilter_iter.GoToReverseBegin(); 
			!chamfer_from_stool_involvement_iter.IsAtReverseEnd() && !min_iter.IsAtReverseEnd() && !curvesFilter_iter.IsAtReverseEnd(); 
			--min_iter, --chamfer_from_stool_involvement_iter, --curvesFilter_iter){
			ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
			float temp_value=chamfer_from_stool_involvement_iter.Get();
			if (temp_value<=chamfercutoff/* && curvesFilter_iter.Get()>-.5 && curvesFilter_iter.Get()<.5)*/) {
				for(int i=-1;i<=1;i++) {
					if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
						for (int j=-1;j<=1;j++) {
							if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
								for (int k=-1;k<=1;k++) {
									if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
										ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
										if (chamfer_from_stool_involvement->GetPixel(temp_index)>=chamfercutoff-1) {
											chamfer_from_stool_involvement->SetPixel(temp_index,chamfercutoff);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//write the modified chamfer map
		WriteITK(chamfer_from_stool_involvement,"chamfer from tissue back 3.mhd");
	}

	for(partialVolume_iter.GoToBegin(),  chamfer_from_stool_involvement_iter.GoToBegin() ;
        !partialVolume_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
		++partialVolume_iter, ++chamfer_from_stool_involvement_iter) 
	{
		partialVolume_iter.Set(partialVolume_iter.Get()*(chamfer_from_stool_involvement_iter.Get())/chamfercutoff);	
	}

	curvesFilter_image.~SmartPointer();
	min.~SmartPointer();
	chamfer_from_stool_involvement.~SmartPointer();
	chamfer_heterogenous.~SmartPointer();

	WriteITK(partialVolume,"partialVolume three.mhd");


	chamfer_from_stool_involvement.~SmartPointer();
	min.~SmartPointer();
	chamfer_heterogenous.~SmartPointer();
	curvesFilter_image.~SmartPointer();
	chamfer_air.~SmartPointer();

	//Does the Gaussian Blurring
    GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
	//GaussianFilterType2::Pointer gaussianFilter = GaussianFilterType2::New();
    double gaussian_sigma[3];
    ImageType::SpacingType input_spacing=input->GetSpacing();   //Gets the spacing

	//adjusted based on spacing of image
    gaussian_sigma[0]=(0.7/input_spacing.GetElement(0));
    gaussian_sigma[1]=(0.7/input_spacing.GetElement(1));
    gaussian_sigma[2]=(0.7/input_spacing.GetElement(2));

	gaussianFilter->SetSigma(gaussian_sigma);      //0.7mm not sure how many voxel this is
	gaussianFilter->SetInputImage(partialVolume);
	std::cerr<<"Preparing output"<<std::endl;

	output.~SmartPointer();
	/*ImageType::Pointer */output = AllocateNewImage( input->GetLargestPossibleRegion());
	output->SetSpacing(input->GetSpacing());
	output_iter = IteratorTypeFloat4WithIndex(output, fullRegion);

	IteratorTypeFloat4WithIndex output_iter_2(output, input->GetLargestPossibleRegion());
	IteratorTypeFloat4WithIndex input_iter_2(input, input->GetLargestPossibleRegion());

	for(output_iter_2.GoToBegin(), input_iter_2.GoToBegin(); !output_iter_2.IsAtEnd() && !input_iter_2.IsAtEnd(); ++output_iter_2, ++input_iter_2) 
	{	
		output_iter_2.Set(input_iter_2.Get());
	}

	
	WriteITK(chamfer_colon,"chamfer_colon_2.mhd");

	//
	//std::cerr<<output_iter.GetRegion().GetIndex()[0]<<" "<<output_iter.GetRegion().GetIndex()[1]<<" "<<output_iter.GetRegion().GetIndex()[2]<<" "<<std::endl;
	//std::cerr<<output_iter.GetRegion().GetSize()[0]<<" "<<output_iter.GetRegion().GetSize()[1]<<" "<<output_iter.GetRegion().GetSize()[2]<<" "<<std::endl;
	//std::cerr<<input_iter.GetRegion().GetIndex()[0]<<" "<<input_iter.GetRegion().GetIndex()[1]<<" "<<input_iter.GetRegion().GetIndex()[2]<<" "<<std::endl;
	//std::cerr<<input_iter.GetRegion().GetSize()[0]<<" "<<input_iter.GetRegion().GetSize()[1]<<" "<<input_iter.GetRegion().GetSize()[2]<<" "<<std::endl;
	//
	//int count_1=0;

	int ending_z = input->GetLargestPossibleRegion().GetSize()[2];

	for(output_iter.GoToBegin(), input_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(), voxel_type_iter.GoToBegin();
        !output_iter.IsAtEnd() && !input_iter.IsAtEnd() && !chamfer_colon_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd();
		++output_iter, ++input_iter, ++chamfer_colon_iter, ++voxel_type_iter) 
	{
	//	count_1++;
		ImageType::IndexType temp_index = output_iter.GetIndex();

		//this is for the padding (5 on each side)
		if ((temp_index[2]<5 || temp_index[2] >= z_offset+5) && (temp_index[2]+1>=ending_z-5 || temp_index[2] < z_offset+z_size-5)) {
			if (voxel_type_iter.Get() == Air || voxel_type_iter.Get() == StoolAir) {
				output_iter.Set(mean[0]);
			}
			if (chamfer_colon_iter.Get()==1 || voxel_type_iter.Get()==TissueStool  || voxel_type_iter.Get()==Stool) {

				float value1 = (input_iter.Get())*gaussianFilter->EvaluateAtIndex(output_iter.GetIndex());
				if ((input_iter.Get()<400 && value1>input_iter.Get()) || input_iter.Get()>=400) {
			////colon mask should be used here
					output_iter.Set(value1);
				} else {
					output_iter.Set(mean[0]);
				}
			}
		}
		//output_iter.Set(2);
	}
	input.~SmartPointer();
	chamfer_colon.~SmartPointer();
	voxel_type.~SmartPointer();
	partialVolume.~SmartPointer();
	input_temp.~SmartPointer();

	//std::cerr<<count_1<<std::endl;

	gaussianFilter.~SmartPointer();
	//Writes out the final output	
	WriteITK(output,"output three.mhd");
	std::cerr<<"Ended"<<std::endl;
	//Sets the input to be the output
	return output;
}

//
//ImageType::Pointer RemoveStool(ImageType::Pointer input) {
//	Modified=false;
//	//Region Setup
//	ImageType::SizeType size = {512,512,40};
//	ImageType::IndexType index = {0,0,100};
//	ImageType::RegionType fullRegion(index,size);
//	//ImageType::RegionType fullRegion=input->GetLargestPossibleRegion();
//	ImageType::Pointer input_temp = AllocateNewImage(fullRegion);
//	IteratorTypeFloat4WithIndex input_temp_iter(input,fullRegion);
//	IteratorTypeFloat4WithIndex input_iter(input_temp,fullRegion);
//	//Copy Image
//	for(input_iter.GoToBegin() , input_temp_iter.GoToBegin() ; !input_iter.IsAtEnd() && !input_temp_iter.IsAtEnd(); 
//		++input_iter, ++input_temp_iter) {
//		input_iter.Set(input_temp_iter.Get());
//	}
//	//print out unmodified
//	WriteITK(input_temp,"input_unmodified.mhd");
//	//Creates an image to store the descriptions
//	VoxelTypeImage::Pointer voxel_type = VoxelTypeImage::New();
//    voxel_type->SetRegions(fullRegion);
//    voxel_type->Allocate();
//	IteratorTypeVoxelType voxel_type_iter(voxel_type,fullRegion);
//
//    //computes the gradient magnitude
//    GradientMagnitudeFilterType::Pointer gradient_filter = GradientMagnitudeFilterType::New();
//    gradient_filter->SetInput(input_temp);
//    gradient_filter->Update();
//    ImageType::Pointer temp = gradient_filter->GetOutput();
//	//Shan's Heuristic Gradient
//	//ImageType::Pointer temp = ComputePseudoGradient(input_temp);
//	IteratorTypeFloat4WithIndex temp_iter(temp,fullRegion);
//	gradient_filter.~SmartPointer();										//Deleting Gradient Filter
//    //Does the Single Material Types Classification (Step1)
//    for (input_iter.GoToBegin(), temp_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !input_iter.IsAtEnd() && !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  ++voxel_type_iter, ++input_iter, ++temp_iter) {
//       voxel_type_iter.Set(SingleMaterialClassification(input_iter.Get()-1024,temp_iter.Get()));        //sets the voxel type
//	}
//	//Filling in Tissue Holes
//	for (voxel_type_iter.GoToBegin(), temp_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
//        ++voxel_type_iter, ++temp_iter)
//    {
//        //creates the morphological images to fill holes
//        switch (voxel_type_iter.Get()) {
//            case Tissue:
//                temp_iter.Set(0);
//                break;
//            case Stool:
//                temp_iter.Set(1);
//                break;
//            case Unclassified:
//                temp_iter.Set(0);
//                break;
//            case Air:
//                temp_iter.Set(0);
//                break;
//        }
//    }
//	//Writing out un filled image
//	WriteITK(temp,"tissue holes.mhd");
//	//Shan's One Voxel Closing Filter
//
//    //StructuringElementType  structuringElement;
//    //structuringElement.SetRadius( 1 );  // 1
//    //structuringElement.CreateStructuringElement();
//	//ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
//	//closingFilter->SetKernel(structuringElement);	 
//	//closingFilter->SetInput(temp);
//	//closingFilter->Update();
//	//temp=closingFilter->GetOutput();
//	//temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
//	//closingFilter.~SmartPointer();
//	//Writing out filled tissue hole image
//	WriteITK(temp,"tissue holes filled.mhd");
//    //Updates Information after morphological operations of tissue in stool
//    //Begins update the unclassified in tissue morphological operation
//	for (temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
//        ++voxel_type_iter, ++temp_iter)
//    {
//        switch(voxel_type_iter.Get()) {
//            case Tissue:
//                //tissue into stool
//                if (temp_iter.Get()==1) {
//                    voxel_type_iter.Set(Stool);			//Modify filled tissue holes
//                }
//                temp_iter.Set(1);
//				break;
//			case Stool:
//                temp_iter.Set(0);
//                break;
//            case Unclassified:
//                temp_iter.Set(0);
//                break;
//            case Air:
//                temp_iter.Set(0);
//                break;
//        }
//    }
//	//Write out unmodified image
//	WriteITK(temp,"unclassified holes.mhd");
//
//	//Shan's Heuristic hole filter
//	//closingFilter = ClosingFilterType::New();
//	//closingFilter->SetKernel(structuringElement);
//	//closingFilter->SetInput(temp);
//	//closingFilter->Update();
//	//temp= closingFilter->GetOutput();
//	//temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
//	//closingFilter.~SmartPointer();
//	//Write out unclassified holes filled
//	WriteITK(temp,"unclassified holes filled.mhd");
//
//    //Updates Information after morphological operations
//    for (temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !temp_iter.IsAtEnd();  
//        ++voxel_type_iter, ++temp_iter) 
//    {
//        switch(voxel_type_iter.Get()) {
//            case Unclassified:
//                if (temp_iter.Get()==1) {
//                    voxel_type_iter.Set(Tissue);			//unclassified into tissue
//                }
//                break;
//        }
//    }
//    //get cubic bspline image interpolation
//    InterpolationType::Pointer input_interpolator = InterpolationType::New();
//    input_interpolator->SetSplineOrder(3);
//    input_interpolator->SetInputImage(input);
//	//Storage for Smax
//    ImageType::Pointer input_smax = AllocateNewImage(fullRegion);
//    IteratorTypeFloat4WithIndex input_smax_iter(input_smax,fullRegion);
//	//Debugging Flags
//	bool StoolAirbool=false;
//	bool TissueStoolbool=false;
//	bool TissueAirbool = false;
//
//	//Runs edge classification with 3 different increments
//    for (voxel_type_iter.GoToBegin(), input_smax_iter.GoToBegin(), temp_iter.GoToBegin(); 
//		!voxel_type_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !temp_iter.IsAtEnd(); 
//		++voxel_type_iter, ++input_smax_iter, ++temp_iter) {
//		if (voxel_type_iter.Get()==Unclassified) {
//			//Initialization
//			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
//			float temp_threshold=-1;
//			VoxelType temp_type = Unclassified;
//			//+/- 1.5, +/-1.0
//			VoxelEdgeClassification(&temp_threshold,&temp_type,1.5,1.0,
//				input_interpolator,input_smax_iter,voxel_index);
//			//+/- 1.0, +/-0.5
//			VoxelEdgeClassification(&temp_threshold,&temp_type,1.0,0.5,
//				input_interpolator,input_smax_iter,voxel_index);
//			//+/- .6, +/-.3
//			VoxelEdgeClassification(&temp_threshold,&temp_type,0.6,0.3,
//				input_interpolator,input_smax_iter,voxel_index);
//			//std::cerr<<"result: "<<temp_type<<" "<<temp_threshold<<std::endl;
//			if (StoolAirbool==false) {
//				StoolAirbool=true;
//				//std::cerr<<"Stool Air Region Detected"<<std::endl;
//			}
//			if (TissueStoolbool==false) {
//				TissueStoolbool=true;
//				//std::cerr<<"Tissue Stool Region Detected"<<std::endl;
//			}
//			if (TissueAirbool==false) {
//				TissueAirbool=true;
//				//std::cerr<<"Tissue Air Region Detected"<<std::endl;
//			}
//			voxel_type_iter.Set(temp_type);
//		}
//		//Set up data for heterogenous stool
//        switch(voxel_type_iter.Get()) {
//			case TissueStool:
//                temp_iter.Set(0);
//                break;
//            case Stool:
//                temp_iter.Set(1);
//                break;
//			case Tissue:
//				temp_iter.Set(3);
//				break;
//			case Air:
//				temp_iter.Set(-2);
//				break;
//			case TissueAir:
//				temp_iter.Set(2);
//				break;
//			case StoolAir:
//				temp_iter.Set(-1);
//				break;
//            default:
//                temp_iter.Set(-3);
//                break;
//        }
//    }
//	//Write out unmodified image
//	WriteITK(temp,"tissue stool.mhd");
//	//get bspline voxeltype interpolation
//    InterpolationType::Pointer voxel_type_interpolator = InterpolationType::New();
//    voxel_type_interpolator->SetSplineOrder(0);
//    voxel_type_interpolator->SetInputImage(temp);
//	for (int i=0;i<3;i++) {
//		bool heteroBool=false;
//		for(voxel_type_iter.GoToBegin(), input_iter.GoToBegin(), temp_iter.GoToBegin();
//			!voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd() && !temp_iter.IsAtEnd();
//			++voxel_type_iter, ++input_iter, ++temp_iter) 
//		{
//			if (!heteroBool) {
//				heteroBool=true;
//				//std::cerr<<"heteogenuous Reclassification ran"<<std::endl;
//			}
//			ImageType::IndexType voxel_index = voxel_type_iter.GetIndex();
//			switch(voxel_type_iter.Get()) {
//				case TissueStool:
//				/*	if (HeterogeneousStoolHeuristic(input_interpolator,voxel_type_interpolator,
//						voxel_index,input_iter.Get())) 
//					{
//						voxel_type_iter.Set(Stool);
//						temp_iter.Set(1);
//					}*/
//					break;
//			}
//		}
//	}
//	//Deletes the Interpolators
//	input_interpolator.~SmartPointer();
//	voxel_type_interpolator.~SmartPointer();
//	//std::cerr<<"set partial volume "<<std::endl;
//    //0: ps
//    //1: pt
//    //2: pa
//    ImageType::Pointer partialVolume = AllocateNewImage(fullRegion);
//    IteratorTypeFloat4WithIndex partialVolume_iter(partialVolume,fullRegion);
//
//
//    ByteImageType::Pointer chamfer_from_stool_involvement=AllocateNewByteImage(fullRegion);
//	IteratorTypeByteWithIndex chamfer_from_stool_involvement_iter(chamfer_from_stool_involvement,fullRegion);
//	//bool random=true;
//	//Sets the partial volumes
//    for(voxel_type_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin(), input_smax_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !input_smax_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
//        ++voxel_type_iter, ++input_iter, ++partialVolume_iter, ++input_smax_iter, ++chamfer_from_stool_involvement_iter) 
//    {
//		//prepares information for chamfer stool calculation
//		if (voxel_type_iter.Get()==Stool || voxel_type_iter.Get()==TissueStool || voxel_type_iter.Get()==StoolAir) {
//			chamfer_from_stool_involvement_iter.Set(1);     //stool
//		} else {
//			chamfer_from_stool_involvement_iter.Set(0);     //non stool is the object
//		}
//		partialVolume_iter.Set(0);
//		float value=0;
//        switch(voxel_type_iter.Get()) {
//			case Tissue:
//				partialVolume_iter.Set(1.0);
//				break;
//			case Stool:
//			case Air:
//				partialVolume_iter.Set(0.0);
//				break;
//            case TissueAir:
//				//std::cerr<<"TA reached"<<std::endl;
//				value = 1.0+(input_iter.Get()-1024)/1000;
//				if (value<=1 && value>=0) {
//					partialVolume_iter.Set(value);	
//				} else {
//					partialVolume_iter.Set(0);	
//				}
//                break;
//            case TissueStool:
//				//std::cerr<<"TS reached"<<std::endl;
//				value=1.0-(input_iter.Get()-1024)/1000;
//				if (value<=1 && value>=0) {
//					partialVolume_iter.Set(value);
//				} else {
//					partialVolume_iter.Set(0);	
//				}
//                break;
//            case StoolAir:
//                partialVolume_iter.Set(0);
//                break;
//        }
//		if (-partialVolume_iter.Get()==std::numeric_limits<float>::infinity()) {
//			//std::cerr<<partialVolume_iter.GetIndex()[0]<<" "<<partialVolume_iter.GetIndex()[1]<<" "<<partialVolume_iter.GetIndex()[2]<<std::endl;
//		}
//    }
//	//Write out partial volume
//	WriteITK(partialVolume,"tissue partial one.mhd");
//	//Write out unmodified image
//	WriteITK(temp,"chamfer from stool input.mhd");
//	ChamferDistanceFilterType::Pointer chamfer_filter = ChamferDistanceFilterType::New();
//	chamfer_filter->SetInput(chamfer_from_stool_involvement);
//	int weights[3]={3,4,5};	//3d distance weight recommended by julian
//	chamfer_filter->SetWeights(weights, weights+3);
//	chamfer_filter->SetDistanceFromObject(false);
//	chamfer_filter->Update();
//
//    chamfer_from_stool_involvement=chamfer_filter->GetOutput();
//	chamfer_from_stool_involvement_iter= IteratorTypeByteWithIndex(chamfer_from_stool_involvement,fullRegion);
//
//	chamfer_filter.~SmartPointer();
//	//Write out chamfer thickness of stool
//	WriteITK(chamfer_from_stool_involvement,"chamfer from stool.mhd");
//
//	//Updates chamfer by setting all stool attached to an 8+ chamfer stool to also 8 chamfer units
//	ImageType::IndexType endIndex = fullRegion.GetIndex();
//	ImageType::IndexType startIndex = fullRegion.GetIndex();	
//	endIndex[0]+=(fullRegion.GetSize()[0]-1);
//	endIndex[1]+=(fullRegion.GetSize()[1]-1);
//	endIndex[2]+=(fullRegion.GetSize()[2]-1);
//	//Forward
//	for(chamfer_from_stool_involvement_iter.GoToBegin(); !chamfer_from_stool_involvement_iter.IsAtEnd(); ++chamfer_from_stool_involvement_iter){
//		ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
//		float temp_value=chamfer_from_stool_involvement_iter.Get();
//		if (temp_value>=8) {
//			for(int i=-1;i<=1;i++) {
//				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
//					for (int j=-1;j<=1;j++) {
//						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
//							for (int k=-1;k<=1;k++) {
//								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
//									ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
//									//std::cerr<<temp_index[0]<<" "<<temp_index[1]<<" "<<temp_index[2]<<" "<<std::endl;
//									if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
//										chamfer_from_stool_involvement->SetPixel(temp_index,8);
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
//	WriteITK(chamfer_from_stool_involvement,"chamfer from stool forward.mhd");
//	//Backward
//	for(chamfer_from_stool_involvement_iter.GoToReverseBegin(); !chamfer_from_stool_involvement_iter.IsAtReverseEnd(); --chamfer_from_stool_involvement_iter){
//		ImageType::IndexType index = chamfer_from_stool_involvement_iter.GetIndex();
//		float temp_value=chamfer_from_stool_involvement_iter.Get();
//		if (temp_value>=8) {
//			for(int i=-1;i<=1;i++) {
//				if (index[0]+i<=endIndex[0] && index[0]+i>=startIndex[0]) {
//					for (int j=-1;j<=1;j++) {
//						if (index[1]+j<=endIndex[1] && index[1]+j>=startIndex[1]) {
//							for (int k=-1;k<=1;k++) {
//								if (index[2]+k<=endIndex[2] && index[2]+k>=startIndex[2]) {
//									ImageType::IndexType temp_index={index[0]+i,index[1]+j,index[2]+k};
//									//std::cerr<<temp_index[0]<<" "<<temp_index[1]<<" "<<temp_index[2]<<" "<<std::endl;
//									if (chamfer_from_stool_involvement->GetPixel(temp_index)>0) {
//										chamfer_from_stool_involvement->SetPixel(temp_index,8);
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
//	WriteITK(chamfer_from_stool_involvement,"chamfer from stool back.mhd");
//
//	ByteImageType::Pointer chamfer_air=chamfer_filter->GetOutput();
//    IteratorTypeByteWithIndex chamfer_air_iter(chamfer_air,fullRegion);
//
//	for(voxel_type_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin(), input_iter.GoToBegin(), partialVolume_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd() && !input_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd();
//        ++voxel_type_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter, ++input_iter, ++partialVolume_iter)
//    {
//		//prepare for computation of chamfer to air
//		if (voxel_type_iter.Get()==Air) {
//            chamfer_air_iter.Set(0);					//distance to AIR
//        } else {
//            chamfer_air_iter.Set(1);					//non important locations
//        }
//
//		float temp_value=chamfer_from_stool_involvement_iter.Get();
//		if (temp_value>=8) {
//			partialVolume_iter.Set(0);
//		} else if (temp_value<8 && temp_value>0) {
//			voxel_type_iter.Set(ThinStool);		//if the chamfer < 8 then it is set to thin stool
//		}
//    }
//
//	//compute chamfer for air
//	chamfer_filter = ChamferDistanceFilterType::New();
//	chamfer_filter->SetInput(chamfer_air);
//	chamfer_filter->SetWeights(weights, weights+3);
//	chamfer_filter->SetDistanceFromObject(false);
//	chamfer_filter->Update();
//
//    chamfer_air=chamfer_filter->GetOutput();
//    chamfer_air_iter=IteratorTypeByteWithIndex(chamfer_air,fullRegion);
//	
//	chamfer_filter.~SmartPointer();
//	//Writes out the chamfer to air
//	WriteITK(chamfer_air,"chamfer from air.mhd");
//	//Normalize the data
//    for(chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
//        !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
//        ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter)
//    {
//        chamfer_from_stool_involvement_iter.Set(chamfer_from_stool_involvement_iter.Get()/8);
//        if (chamfer_air_iter.Get()>=8) {
//            chamfer_air_iter.Set(1);
//        } else {
//            chamfer_air_iter.Set(chamfer_air_iter.Get()/8);
//        }
//    }
//
//    //Calculating partial air
//    for(voxel_type_iter.GoToBegin(), partialVolume_iter.GoToBegin(), chamfer_air_iter.GoToBegin(), chamfer_from_stool_involvement_iter.GoToBegin();
//        !voxel_type_iter.IsAtEnd() && !partialVolume_iter.IsAtEnd() && !chamfer_air_iter.IsAtEnd() && !chamfer_from_stool_involvement_iter.IsAtEnd();
//		++voxel_type_iter, ++partialVolume_iter, ++chamfer_air_iter, ++chamfer_from_stool_involvement_iter) 
//	{
//        if(voxel_type_iter.Get()==ThinStool) {
//			float value=0;
//			float ps=1/2*(1+vnl_erf((chamfer_from_stool_involvement_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
//			float pa=1/2*(1+vnl_erf((chamfer_air_iter.Get()-0.5)/(CDF_SIGMA*sqrtf(2))));
//			value=1-ps-pa;
//			partialVolume_iter.Set(value);
//        }
//	}
//
//	//Write out the partial tissue with thinstool modification done
//	WriteITK(partialVolume,"tissue partial two.mhd");
//
//	//Does the Gaussian Blurring
//    GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();
//    double gaussian_sigma[3];
//    ImageType::SpacingType input_spacing=input->GetSpacing();   //Gets the spacing
//	
//	//adjusted based on spacing of image
//    gaussian_sigma[0]=(0.7/input_spacing.GetElement(0));
//    gaussian_sigma[1]=(0.7/input_spacing.GetElement(1));
//    gaussian_sigma[2]=(0.7/input_spacing.GetElement(2));
//	gaussianFilter->SetSigma(gaussian_sigma);      //0.7mm not sure how many voxel this is
//    gaussianFilter->SetInputImage(partialVolume);
//	
//	//Prepares the output image
//	ImageType::Pointer output = AllocateNewImage(fullRegion);
//	IteratorTypeFloat4WithIndex output_iter(output,fullRegion);
//
//    //Modifies to intensity for output
//    for(output_iter.GoToBegin(), input_iter.GoToBegin(), temp_iter.GoToBegin(), voxel_type_iter.GoToBegin();
//		!output_iter.IsAtEnd() && !input_iter.IsAtEnd() && !temp_iter.IsAtEnd() &&  !voxel_type_iter.IsAtEnd();
//		++output_iter, ++input_iter, ++temp_iter, ++voxel_type_iter) 
//    {	
//        ImageType::IndexType index =input_iter.GetIndex();
//		output_iter.Set(gaussianFilter->EvaluateAtIndex(index)*(input_iter.Get()));
//		if (-output_iter.Get()== std::numeric_limits<float>::infinity()) {
//			output_iter.Set(0);
//		}
//
//		//Prepare for last closing of air enclosed voxels
//		if (voxel_type_iter.Get()==Air) {
//			temp_iter.Set(1);
//		} else {
//			 temp_iter.Set(0);
//		}
//    }
//
//	//Delete the Gaussian filter
//	gaussianFilter.~SmartPointer();
//
//	//Writes out the output preclosing
//	WriteITK(output,"output one.mhd");
//
//	//Shan's Heuristic closing
//
//	//Last Morphological Processing to remove anything contained completely in air
//	//closingFilter = ClosingFilterType::New();
//	//closingFilter->SetKernel(structuringElement);
//	//closingFilter->SetInput(temp);
//	//closingFilter->Update();
//	//temp=closingFilter->GetOutput();
//	//temp_iter=IteratorTypeFloat4WithIndex(temp,fullRegion);
//	//closingFilter.~SmartPointer();
//
//	//Updates the information for final output
//    for(temp_iter.GoToBegin(), output_iter.GoToBegin();
//        !temp_iter.IsAtEnd() && !output_iter.IsAtEnd();
//        ++temp_iter, ++output_iter) 
//    {
//        if (temp_iter.Get()==1) {
//            output_iter.Set(0); //sets for air
//        }
//		if (-output_iter.Get()== std::numeric_limits<float>::infinity()) {
//			output_iter.Set(0);
//		}
//    }
//
//	//Writes out the final output	
//	WriteITK(output,"output two.mhd");
//	//Sets the input to be the output
//	return output;
//}
//
//
//






void VoxelEdgeClassification(float * threshold, VoxelType * previous, double d2, double d1,

                                      InterpolationType::Pointer input_interpolator, 

                                      IteratorTypeFloat4WithIndex input_smax,

                                      ImageType::IndexType index) 

{
	ImageType::RegionType fullRegion = input_interpolator->GetInputImage()->GetLargestPossibleRegion();
    InterpolationType::ContinuousIndexType offset_index;
    offset_index[0]=index[0];
    offset_index[1]=index[1];
    offset_index[2]=index[2];
	ImageType::IndexType endIndex = input_interpolator->GetEndIndex();
	ImageType::IndexType startIndex = input_interpolator->GetStartIndex();
    InterpolationType::CovariantVectorType gradient = input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index);
    gradient.Normalize();
    //stores the intensity
    float temp_intensity[5];
    float temp_gradient_magnitude[5];
	//stores the target voxel
    temp_intensity[2]=input_interpolator->EvaluateAtContinuousIndex(offset_index)-1024;
    temp_gradient_magnitude[2]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //+d1
    offset_index[0]=index[0]+gradient[0]*d1;
    offset_index[1]=index[1]+gradient[1]*d1;
    offset_index[2]=index[2]+gradient[2]*d1;
    temp_intensity[3]=input_interpolator->EvaluateAtContinuousIndex(offset_index)-1024;
    temp_gradient_magnitude[3]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //-d1
    offset_index[0]=index[0]-gradient[0]*d1;
    offset_index[1]=index[1]-gradient[1]*d1;
    offset_index[2]=index[2]-gradient[2]*d1;
    temp_intensity[1]=input_interpolator->EvaluateAtContinuousIndex(offset_index)-1024;
    temp_gradient_magnitude[1]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //+d2
    offset_index[0]=index[0]+gradient[0]*d2;
    offset_index[1]=index[1]+gradient[1]*d2;
    offset_index[2]=index[2]+gradient[2]*d2;
    temp_intensity[4]=input_interpolator->EvaluateAtContinuousIndex(offset_index)-1024;
    temp_gradient_magnitude[4]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
    //-d2
    offset_index[0]=index[0]-gradient[0]*d2;
    offset_index[1]=index[1]-gradient[1]*d2;
    offset_index[2]=index[2]-gradient[2]*d2;
    temp_intensity[0]=input_interpolator->EvaluateAtContinuousIndex(offset_index)-1024;
    temp_gradient_magnitude[0]=input_interpolator->EvaluateDerivativeAtContinuousIndex(offset_index).GetNorm();
	//std::cerr<<"vector created"<<std::endl;
	float min_distance=*threshold;
	float stool_tissue_Smax;
	float stool_air_Smax;
	//Tissue Stool
	stool_tissue_Smax=ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
	stool_air_Smax=Stool_Air_ComputeSmax(temp_intensity,temp_gradient_magnitude, 5);
    float distanceTS=AverageTissueStoolDist(stool_tissue_Smax, temp_intensity,temp_gradient_magnitude);
	float distanceSA=1.2*AverageStoolAirDist(stool_air_Smax, temp_intensity,temp_gradient_magnitude);
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

	if (smax<150 && Modified) {
		smax=stool_air_Smax;
		distance=distanceTA;
		voxel_type=TissueAir;
	}


    if (min_distance>=distance || min_distance==-1 || *previous==Unclassified) {
		*threshold=distance;
        *previous=voxel_type;
        input_smax.Set(smax);
    }
}
float Stool_Air_ComputeSmax(float intensity[], float gradient_magnitude[], int size) {
	//std::cerr<<"Value ";
	float intensity2[7];
	for (int i=0;i<size;i++) {
		intensity2[i]=intensity[i]+1024;
	//	std::cerr<<intensity2[i]<<" ";
	}
	//std::cerr<<std::endl;
	return ComputeSmax(intensity2,gradient_magnitude, size)-1024;
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
	if (Modified) {
		coefficients[0]= -6.0/1000;
		coefficients[1]= -6.0;
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
	if (Modified) {
		if ((input_pixel >=400 && input_gradient_pixel<=0.8*input_pixel) || input_pixel >= 1000) {
			return Stool;
		} else if (input_pixel<=-800 && input_gradient_pixel<=250) {
			return Air;
		} else if (input_pixel<=150  && input_pixel>=-250 && input_gradient_pixel<=300) {
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

void WriteITK(ImageType::Pointer image, char * name) {
 //   typedef itk::ImageFileWriter< ImageType >  WriterType;
 //   WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName(name);
	//writer->SetInput(image);
	std::cerr<<"writing: "<<name<<std::endl;
	//writer->Update();
	//writer.~SmartPointer();
}

void WriteITK(ByteImageType::Pointer image, char * name) {
 //   typedef itk::ImageFileWriter< ByteImageType >  WriterType;
 //   WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName(name);
	//writer->SetInput(image);
	std::cerr<<"writing: "<<name<<std::endl;
	//writer->Update();
	//writer.~SmartPointer();
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