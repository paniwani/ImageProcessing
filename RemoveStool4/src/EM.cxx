vnl_matrix<float> GetNeighbor(ArrayImageType::Pointer partialVector, ImageType::IndexType index) {
    ImageType::RegionType fullregion = partialVector->GetLargestPossibleRegion();
    vnl_matrix<float> temp_return(3,6);
    float filler[3]={0,0,0};
    for(int i=0;i<3;i++) {
        ImageType::IndexType temp_index_1(index);
        temp_index_1[i]+=1;
        ArrayType data_1(filler);
        if (fullregion.IsInside(temp_index_1)) {
            data_1 = partialVector->GetPixel(temp_index_1);    //account for image issues
        }

        ImageType::IndexType temp_index_2(index);
        temp_index_2[i]-=1;
        ArrayType data_2(filler);
        if (fullregion.IsInside(temp_index_2)) {
            data_2 = partialVector->GetPixel(temp_index_2);
        } 

		//std::cout << "pause" << std::endl;

        for (int j=0;j<3;j++) {
            temp_return.put(j,i,data_1[j]);
			temp_return.put(j,i+3,data_2[j]);
        }
    }

    return temp_return;
}

double Probability(double Y, double mean, double variance,  double current_partial, double local_variance, double local_mean) {  

	if (mean == Y && variance == 0)
	{
		return 1;
	} else {
		if (local_mean == 0 && local_variance == 0)
		{
			return 0;
		} else {
			return exp(-vnl_math_sqr(Y-mean)/(2*variance))/(sqrt(2*PI*variance))*exp(-vnl_math_sqr(current_partial-local_mean)/(2*local_variance))/(sqrt(2*PI*local_variance));
		}
	}
}

vnl_vector<float> expectation(double Y, double mean[], double variance[], float weight[],  vnl_matrix<float> neighbor, float current_partial[]) {	
	int nKernal=3;
	float pFK[3]={0,0,0};	//probabilities of each class: air/tissue/stool per voxel
	float pKF[4]={0,0,0};   //""	"", [3] = sum of probabilities of air/tissue/stool
	float sumPFK=0.0;
	for (int i=0;i<3;i++) {
		double total_neighbor=0;
		
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

void EM(ArrayImageType::Pointer &partial, ByteImageType::Pointer &colon, ImageType::Pointer &input)
{
	IteratorType inputIt(input,REGION);
	ByteIteratorType colonIt(colon,REGION);
	ArrayIteratorType partialIt(partial,REGION);

	int counter=0;							//used to total certain values (description given when set)
    double mean[3]={0, 0, 0};				//stores the mean of the air/tissue/stool classes
    double sum[3]={0, 0, 0};				//stores the total # of partials of the three classes
    double variance[3]={0, 0, 0};			//stores the variance of the air/tissue/stool classes
	float weight[3]={0,0,0};				//stores the weights of the air/tissue/stool classes

	/*std::ofstream em;
	em.open("debugEM.csv");*/

	//int count = 0;

	// Computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
	for(inputIt.GoToBegin(), partialIt.GoToBegin(), colonIt.GoToBegin();
		!inputIt.IsAtEnd() && !partialIt.IsAtEnd() && !colonIt.IsAtEnd();
		++inputIt, ++partialIt, ++colonIt) 
	{
		if (colonIt.Get()==255) {		

			/*if ( count > 34585 )
			{
				std::cout << "DEBUG TIME" << std::endl;
			}*/

			ArrayType p = partialIt.Get();
			for (int i=0; i<3; i++)
			{

				sum[i] += p[i];
				mean[i] += p[i]*inputIt.Get();
			}

			/*em << sum[0] << "," << sum[1] << "," << sum[2] << "\n";

			count++;*/
		}
	}

	/*std::cout << sum[0] << " " << sum[1] << " " << sum[2] << std::endl;

	em.close();*/

	for (int i=0;i<3;i++) { mean[i]=mean[i]/sum[i]; } 

	// Compute variance and weights
	for(inputIt.GoToBegin(), partialIt.GoToBegin(), colonIt.GoToBegin();
		!inputIt.IsAtEnd() && !partialIt.IsAtEnd() && !colonIt.IsAtEnd();
		++inputIt, ++partialIt, ++colonIt) 
	{
		if (colonIt.Get()==255) {		
			ArrayType p = partialIt.Get();
			for (int i=0; i<3; i++)
			{
				variance[i] += p[i]*vnl_math_sqr(inputIt.Get()-mean[i]);
			}
		}
	}

	double sum_all = sum[0]+sum[1]+sum[2];

	for (int i=0;i<3;i++) 
	{ 
		variance[i]=variance[i]/sum[i];
		weight[i]=sum[i]/sum_all;
	} 

	std::cerr<<std::endl;
	std::cerr<<"EM0"<<std::endl;
	std::cerr<<"Mean: "<<mean[0]<<" "<<mean[1]<<" "<<mean[2]<<std::endl;
	std::cerr<<"Variance: "<<variance[0]<<" "<<variance[1]<<" "<<variance[2]<<std::endl;
	std::cerr<<"Std Dev: "<<sqrt(variance[0])<<" " <<sqrt(variance[1])<<" "<<sqrt(variance[2])<<std::endl;
	std::cerr<<"Weight: "<<weight[0]<<" "<<weight[1]<<" "<<weight[2]<<std::endl;

	ImageType::SpacingType spacing = input->GetSpacing();
	
	// Compute iterations of the Maximization Algorithm
    for (int emNum=0;emNum<5;emNum++) // used 20 in original test case
	{
		double sum_temp[3]={0,0,0};
        double variance_temp[3]={0,0,0};
		double mean_temp[3]={0,0,0};

		// Duplicate partials
		typedef itk::ImageDuplicator<ArrayImageType> ImageDuplicatorType;
		ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
		duplicator->SetInputImage(partial);
		duplicator->Update();
		ArrayImageType::Pointer partial2 = duplicator->GetOutput();
		ArrayIteratorType partial2It(partial2,REGION);

		for (inputIt.GoToBegin(), partialIt.GoToBegin(), partial2It.GoToBegin(), colonIt.GoToBegin();
			!inputIt.IsAtEnd() && !partialIt.IsAtEnd() && !partial2It.IsAtEnd() && !colonIt.IsAtEnd();
			++inputIt, ++partialIt, ++partial2It, ++colonIt) 
		{
			if (colonIt.Get()==255) {		
				ArrayType p = partialIt.Get();		//retrieves the partial informations
				float Z[3]={p[0],p[1],p[2]};

				if (p[1] > 0 && p[1] < 1)
				{
					ImageType::IndexType idx = inputIt.GetIndex();

					vnl_vector<float> Z_update=expectation(inputIt.Get(),mean, variance, weight, GetNeighbor(partial,idx), Z);	//updates the partial values
					
					ArrayType p2;

					for (int i=0; i<3; i++)
						p2[i] = Z_update[i];

					// set tolerance to allow uncertain probabilites to update to certain
					for (int i=0; i<3; i++)
					{
						p2[i] = p2[i] > 0.001 ? p2[i] : 0;
						p2[i] = p2[i] < 0.999 ? p2[i] : 1;
					}

					partial2It.Set(p2);
				}
			}
        }

		partial = partial2;

		partialIt = ArrayIteratorType(partial,REGION);

		// Computes the mean (expectation) for each class by using sum(partial[i]*value[i])/sum(partial[i])
		for(inputIt.GoToBegin(), partialIt.GoToBegin(), colonIt.GoToBegin();
			!inputIt.IsAtEnd() && !partialIt.IsAtEnd() && !colonIt.IsAtEnd();
			++inputIt, ++partialIt, ++colonIt) 
		{
			if (colonIt.Get()==255) 
			{		
				ArrayType p = partialIt.Get();
				for (int i=0; i<3; i++)
				{
					sum_temp[i] += p[i];
					mean_temp[i] += p[i]*inputIt.Get();
				}
			}
		}

		for (int i=0;i<3;i++) { mean_temp[i]=mean_temp[i]/sum_temp[i]; } 

		// Compute variance and weights
		for(inputIt.GoToBegin(), partialIt.GoToBegin(), colonIt.GoToBegin();
			!inputIt.IsAtEnd() && !partialIt.IsAtEnd() && !colonIt.IsAtEnd();
			++inputIt, ++partialIt, ++colonIt) 
		{
			if (colonIt.Get()==255) {		
				ArrayType p = partialIt.Get();
				for (int i=0; i<3; i++)
				{
					variance_temp[i] += p[i]*vnl_math_sqr(inputIt.Get()-mean_temp[i]);
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

		/*std::stringstream ss;
		ss<<emNum+1;
		ss<<"\t"<<mean[0]<<"\t"<<mean[1]<<"\t"<<mean[2]<<"\t"<<variance[0]<<"\t"<<variance[1]<<"\t"<<variance[2]<<"\t"<<weight[0]<<"\t"<<weight[1]<<"\t"<<weight[2]<<"\n";
		em<<ss.str();*/

		std::stringstream ss2;
		ss2<<"EM"<<emNum+1<<".nii";

		//Write(partial,ss2.str());
    }

	Write(partial,"EMPartial.nii");

	//debug.close();
}