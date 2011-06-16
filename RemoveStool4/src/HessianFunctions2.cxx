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

float weight_dark(float ls, float lt, float alpha, float gamma)
{
	if ( (ls >= lt) && (lt >= 0) )
	{
		return pow( 1 - ( lt / abs(ls) ) , gamma );
	} else if ( ( (-abs(ls)/alpha) < lt) & (lt < 0) )
	{
		return pow( 1 + alpha*( lt / abs(ls) ), gamma);
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

float omega_dark(float ls, float lt, float gamma)
{
	if ( (ls >= lt) && (lt > 0) )
	{
		return pow( lt / ls , gamma );
	}

	return 0;
}

float S_line(EigenValueArrayType lambda, float alpha, float gamma)
{
	if ( (lambda[2] <= lambda[1]) && (lambda[1] < 0) )
	{
		return abs(lambda[2])*omega(lambda[1],lambda[2],gamma)*weight(lambda[0],lambda[1],alpha,gamma);
	}

	return 0;
}

float S_line_dark(EigenValueArrayType lambda, float alpha, float gamma)
{
	if ( (lambda[0] >= lambda[1]) && (lambda[1] > 0) )
	{
		return abs(lambda[0])*omega_dark(lambda[0],lambda[1],gamma)*weight_dark(lambda[1],lambda[2],alpha,gamma);
	}

	return 0;
}

float S_blob(EigenValueArrayType lambda, float gamma)
{
	if ( (lambda[2] <= lambda[1]) && (lambda[1] <= lambda[0]) && (lambda[0] < 0) )
	{
		return abs(lambda[2])*omega(lambda[1], lambda[2], gamma)*omega(lambda[0], lambda[1], gamma);
	} 

	return 0;
}

float S_blob_dark(EigenValueArrayType lambda, float gamma)
{
	if ( lambda[0] >= lambda[1] && lambda[1] >= lambda[2] && lambda[2] > 0 )
	{
		return abs(lambda[0])*omega_dark(lambda[0], lambda[1], gamma)*omega_dark(lambda[1], lambda[2], gamma);
	} 

	return 0;
}

float S_sheet(EigenValueArrayType lambda, float alpha, float gamma)
{
	if (lambda[2] < 0)
	{
		return abs(lambda[2])*weight(lambda[1], lambda[2], alpha, gamma)*weight(lambda[0], lambda[2], alpha, gamma);
	}

	return 0;
}

float S_sheet_dark(EigenValueArrayType lambda, float alpha, float gamma)
{
	if (lambda[0] > 0)
	{
		return abs(lambda[0])*weight_dark(lambda[0], lambda[1], alpha, gamma)*weight_dark(lambda[0], lambda[2], alpha, gamma);
	}

	return 0;
}

float S_fold(EigenValueArrayType lambda, float alpha, float gamma)
{
	if (lambda[2] < 0)
	{
		return abs(lambda[2])*weight(lambda[1], lambda[2], alpha, gamma)*weight(lambda[0], lambda[1], alpha, gamma);
	}
	
	return 0;
}

float thinness( EigenValueArrayType l )
{
	return exp((abs(l[0]/l[1]) - 1)*(abs(l[1])+abs(l[2])))*exp(-l[0]/l[2]);
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

//FloatImageType::Pointer SatoResponse(ImageType::Pointer &inputS, ByteImageType::Pointer &colon, double alpha, double gamma, float T0)
//{	
//	// Cast input to float
//	typedef itk::CastImageFilter<ImageType,FloatImageType> CastType;
//	CastType::Pointer caster = CastType::New();
//	caster->SetInput(inputS);
//	caster->Update();
//	FloatImageType::Pointer input = caster->GetOutput();
//
//	// Set intensity and lambda thresholds
//	float intThreshold[2] = {T0, 1500};
//	float labmdaThreshold = 0;
//
//	// Get region
//	ImageType::RegionType region = input->GetLargestPossibleRegion();
//
//	// Create response image
//	FloatImageType::Pointer output = FloatImageType::New();
//	output->SetSpacing( input->GetSpacing() );
//	output->SetRegions( region );
//	output->CopyInformation( input );
//	output->Allocate();
//	output->FillBuffer(0);
//
//	// Get eigenvalues
//	FloatIteratorType inputIt(input,region);
//	FloatIteratorType outputIt(output,region);
//	ByteIteratorType colonIt(colon,region);
//
//	// Compute hessian across 2 sigma scales
//	double sigma[2] = {0.56,1.4};
//	//double sigma[2] = {3,5};
//
//	for (int k=0; k<2; k++)
//	{
//		typedef itk::HessianRecursiveGaussianImageFilter<FloatImageType> HessianType;
//		HessianType::Pointer hessianFilter = HessianType::New();
//		hessianFilter->SetInput(input);
//		hessianFilter->SetNormalizeAcrossScale(true);
//		hessianFilter->SetSigma(sigma[k]);
//		hessianFilter->Update();
//		itk::ImageRegionConstIterator<HessianType::OutputImageType> hessianIt(hessianFilter->GetOutput(),region);
//
//		hessianIt.GoToBegin();
//		inputIt.GoToBegin();
//		outputIt.GoToBegin();
//		colonIt.GoToBegin();
//
//		int count=0;
//
//		while (!hessianIt.IsAtEnd()) 
//		{
//			float I = (float) inputIt.Get();
//
//			if ( I >= intThreshold[0] && I <= intThreshold[1] && colonIt.Get() != 0)
//			{	
//				// Get eigenvalues
//				EigenValueArrayType lambda;
//				
//				hessianIt.Get().ComputeEigenValues(lambda);
//				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);
//				
//				if (lambda[0] >= labmdaThreshold)
//				{
//					/*for (int i=0; i<3; i++)
//					{
//						LambdaIterVector[i].Set( lambda[i] );
//					}*/
//
//					/*w1_iter.Set( weight_dark(Lambda1, Lambda2, alpha, gamma) );
//					w2_iter.Set( weight_dark(Lambda1, Lambda3, alpha, gamma) );*/
//
//					float val = S_sheet_dark(lambda, alpha, gamma);
//
//					if ( val > outputIt.Get() )
//					{
//						outputIt.Set( val );		
//					}
//				}
//
//
//			}
//
//			if (++count % 500000 == 0)
//			{
//				std::cout << count << std::endl;
//			}	
//
//			++hessianIt;
//			++inputIt;
//			++outputIt;
//			++colonIt;
//		}
//
//		std::stringstream ss;
//		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_sato_output.nii";
//		Write(output,ss.str());
//
//	}
//	
//	std::stringstream ss;
//	ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".nii";
//	//Write(output, ss.str());
//
//	return output;
//}

FloatImageType::Pointer ResampleImage(FloatImageType::Pointer input, FloatImageType::SpacingType output_spacing)
{
	typedef itk::ResampleImageFilter<FloatImageType,FloatImageType> ResampleImageFilterType;
	ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();

	resampleFilter->SetDefaultPixelValue( -1024 );

	FloatImageType::SpacingType input_spacing = input->GetSpacing();
	FloatImageType::SizeType input_size = input->GetLargestPossibleRegion().GetSize();

	FloatImageType::SizeType output_size;

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


FloatImageType::Pointer SatoResponse(ImageType::Pointer input_aniso_short, double alpha, double gamma, std::vector<float> sigmaVector)
{	
	// Cast input to float
	typedef itk::CastImageFilter<ImageType,FloatImageType> CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput(input_aniso_short);
	caster->Update();
	FloatImageType::Pointer input_aniso = caster->GetOutput();

	// Set intensity and lambda thresholds
	float intensity_threshold_lower = 200;
	float intensity_threshold_upper = 1200;
	float lambda1_threshold = 0;

	// Make image isotropic
	FloatImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();

	FloatImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	FloatImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	Write(input,"input_isotropic.nii");

	// Get input info
	FloatImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	FloatImageType::Pointer output = AllocateFloatImage(input);
	output->SetSpacing(iso_spacing);

	FloatIteratorType input_iter(input,region);
	FloatIteratorType output_iter(output,region);

	// Compute hessian across sigma scales
	//double sigma[2] = {0.56, 1.4};
	//double sigma[7] = {0.5,1,2,4,6,8,10};

	for (int k=0; k<sigmaVector.size(); k++)
	{
		// Reset buffer
		//output->FillBuffer(0);

		//// Keep track of weight functions
		//FloatImageType::Pointer w1 = AllocateFloatImage(input);
		//FloatImageType::Pointer w2 = AllocateFloatImage(input);

		//w1->SetSpacing( iso_spacing );
		//w2->SetSpacing( iso_spacing );

		//w1->FillBuffer(0);
		//w2->FillBuffer(0);

		//FloatIteratorType w1_iter(w1,region);
		//FloatIteratorType w2_iter(w2,region);

		//w1_iter.GoToBegin();
		//w2_iter.GoToBegin();

		// Compute smoothed Hessian
		typedef itk::HessianRecursiveGaussianImageFilter<FloatImageType> HessianGaussianFilterType;
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigmaVector[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		// Create helper images and iters to visualize algorithm
		std::vector< FloatImageType::Pointer > LambdaImageVector(3);
		std::vector< itk::ImageRegionIterator<FloatImageType> > LambdaIterVector(3);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateFloatImage(input);
			LambdaImageVector[i]->SetSpacing( iso_spacing );
			LambdaImageVector[i]->FillBuffer(0);
			LambdaIterVector[i] = itk::ImageRegionIterator<FloatImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			// Set submerged and non-submerged conditions

			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper )
			{				
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				double Lambda1 = lambda[0];
				double Lambda2 = lambda[1];
				double Lambda3 = lambda[2];

				for (int i=0; i<3; i++)
				{
					LambdaIterVector[i].Set( lambda[i] );
				}

				/*if (Lambda1 >= 0)
				{
					w1_iter.Set( weight_dark(Lambda1, Lambda2, alpha, gamma) );
					w2_iter.Set( weight_dark(Lambda1, Lambda3, alpha, gamma) );
				}*/

				//float eig[3] = {Lambda1, Lambda2, Lambda3};
				//double val = vnl_math_max( S_blob_dark(lambda, gamma), S_sheet_dark(lambda, alpha, gamma) );
				//double val = S_blob_dark(lambda, gamma);
				double val = S_sheet_dark(lambda, alpha, gamma);

				if ( val > output_iter.Get() )
				{
					output_iter.Set( val );		
				}
			}

			if (++count % 500000 == 0)
			{
				std::cout << count << std::endl;
			}	

			++hessian_iter;
			++input_iter;
			++output_iter;

			/*++w1_iter;
			++w2_iter;*/
			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
		}

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigmaVector[k] << "_lambda_" << i << ".nii";
			Write(  LambdaImageVector[i], ss.str());
		}

		//std::stringstream ss;
		//ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_w1.nii";
		////Write( ResampleImage( w1, aniso_spacing ) , ss.str() );

		//ss.str("");
		//ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_w2.nii";
		////Write( ResampleImage( w2, aniso_spacing ) , ss.str() );

		std::stringstream ss;
		ss << "sigma_" << sigmaVector[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_sato_output.nii";
		Write( output, ss.str() );

	}

	FloatImageType::Pointer output_aniso = ResampleImage(output, aniso_spacing);

	std::stringstream ss;
	ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".nii";
	Write(output_aniso, ss.str());

	return output_aniso;
}


FloatImageType::Pointer SatoResponse2(ImageType::Pointer input_aniso_short, double alpha, double gamma)
{	
	// Cast input to float
	typedef itk::CastImageFilter<ImageType,FloatImageType> CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput(input_aniso_short);
	caster->Update();
	FloatImageType::Pointer input_aniso = caster->GetOutput();

	// Set intensity and lambda thresholds
	float intensity_threshold_lower = 200;
	float intensity_threshold_upper = 1200;
	float lambda1_threshold = 0;

	// Make image isotropic
	FloatImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();

	FloatImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	FloatImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	Write(input,"input_isotropic.nii");

	// Get input info
	FloatImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response images
	FloatImageType::Pointer output1 = AllocateFloatImage(input);
	output1->SetSpacing(iso_spacing);

	FloatImageType::Pointer output2 = AllocateFloatImage(input);
	output2->SetSpacing(iso_spacing);



	FloatIteratorType input_iter(input,region);
	FloatIteratorType output1_iter(output1,region);
	FloatIteratorType output2_iter(output2,region);


	// Compute hessian across sigma scales
	double sigma[2] = {0.56, 1.4};

	for (int k=0; k<2; k++)
	{
		// Compute smoothed Hessian
		typedef itk::HessianRecursiveGaussianImageFilter<FloatImageType> HessianGaussianFilterType;
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output1_iter.GoToBegin();
		output2_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			// Set submerged and non-submerged conditions
			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper )
			{				
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				//double val = S_blob_dark(lambda, gamma);
				double val = S_sheet_dark(lambda, alpha, gamma);

				if (k == 0)
				{
					output1_iter.Set(val);
				} else {
					output2_iter.Set(val);
				}
			}

			if (++count % 500000 == 0)
			{
				std::cout << count << std::endl;
			}	

			++hessian_iter;
			++input_iter;
			++output1_iter;
			++output2_iter;
		}

		//std::stringstream ss;
		//ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_sato_output.nii";
		//Write( output, ss.str() );

	}

	// Get max of each response image scaled from [0,1]
	typedef itk::RescaleIntensityImageFilter<FloatImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(output1);
	rescaler->SetOutputMaximum(1);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	output1 = rescaler->GetOutput();

	rescaler = RescalerType::New();
	rescaler->SetInput(output2);
	rescaler->SetOutputMaximum(1);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	output2 = rescaler->GetOutput();

	typedef itk::MaximumImageFilter<FloatImageType> MaximumFilterType;
	MaximumFilterType::Pointer maxFilter = MaximumFilterType::New();
	maxFilter->SetInput1(output1);
	maxFilter->SetInput2(output2);
	maxFilter->Update();

	FloatImageType::Pointer output_aniso = ResampleImage(maxFilter->GetOutput(), aniso_spacing);

	std::stringstream ss;
	ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".nii";
	Write(output_aniso, ss.str());

	return output_aniso;
}



bool compare(ptype a, ptype b)
{
	if(a.size > b.size)
	{
		return true;
	}
	else return false;
}

bool compareintensity(ptype a, ptype b)
{
	if(a.intensity < b.intensity)
	{
		return true;
	}
	else return false;
}

void Relabel(LabelImageType::Pointer &image, unsigned int minSize=0)
{
	LabelIteratorType it(image, image->GetLargestPossibleRegion());
	unsigned int maxintensity = 0;
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxintensity)
		{
			maxintensity = it.Get();
		}
	}
	
	std::vector<ptype> count;
	for(unsigned int i = 1; i <= maxintensity; i++)
	{
		ptype a;
		a.intensity = i;
		a.size = 0;
		count.push_back(a);
	}

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
		{
			++count[it.Get() - 1].size;
		}
	}

	sort(count.begin(), count.end(), compare);

	for(unsigned int i = 0; i < count.size(); i++)
	{
		count[i].order = i + 1;
	}

	sort(count.begin(), count.end(), compareintensity);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0 && count[it.Get()-1].size > minSize)
		{
			it.Set(count[it.Get() - 1].order);
		} else {
			it.Set(0);
		}
	}
}

//void HessianAnalysis(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, ArrayImageType::Pointer &partial, PixelType tissueStoolThreshold)
//{
//	// get region
//	ImageType::RegionType region = input->GetLargestPossibleRegion();
//		
//	// get hessian response
//	double alpha = 0.25;
//	double gamma = 0.5;
//
//	FloatImageType::Pointer hessian = SatoResponse(input,alpha,gamma);
//	Write(hessian,"hessian.nii");
//
//	// mask with stool and tissue
//	hessian = Mask(hessian,BinaryOr(BinaryThreshold(vmap,Stool), BinaryThreshold(vmap,Tissue)));
//	Write(hessian,"hessianMaskedStoolAndTissue.nii");
//
//	//// mask with colon
//	//hessian = Mask(hessian,colon);
//	//Write(hessian,"hessianMaskedColon.nii");
//
//	// rescale to [0,1]
//	Rescale(hessian,0,1);
//	Write(hessian,"hessianRescaled.nii");
//
//	// remove lowest 20%
//	ByteImageType::Pointer hbin = BinaryThreshold(hessian,.2);
//	Write(hbin,"hbin20.nii");
//
//	hessian.~SmartPointer();
//
//	// ------------------------------------------------------------------------
//	// for each hessian component, count number of times it touches largest component of tissue
//	// ------------------------------------------------------------------------
//
//	// get connected components		
//	typedef itk::ConnectedComponentImageFilter<ByteImageType, LabelImageType> ConnectedType;
//	ConnectedType::Pointer connecter = ConnectedType::New();
//	connecter->SetInput(hbin);
//	connecter->SetBackgroundValue(0);
//	connecter->Update();
//	LabelImageType::Pointer cc = connecter->GetOutput();
//	
//	Relabel(cc,3);
//
//	LabelIteratorType ccIt(cc,region);
//
//	Write(cc,"ccRelabel.nii");
//	
//	/*typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType> RelabelType;
//	RelabelType::Pointer relabeler = RelabelType::New();
//	relabeler->SetInput(connecter->GetOutput());
//	relabeler->SetMinimumObjectSize(3);
//	relabeler->Update();
//
//	unsigned long originalNumOfObjects = relabeler->GetOriginalNumberOfObjects();
//	unsigned long numOfObjects = relabeler->GetNumberOfObjects();
//
//	std::cout << "Original # of objects: " << originalNumOfObjects << std::endl;
//	std::cout << "Number of objects: " << numOfObjects << std::endl;*/
//	
//	unsigned int numOfObjects = 0;
//
//	for (ccIt.GoToBegin(); !ccIt.IsAtEnd(); ++ccIt)
//	{
//		if (ccIt.Get() > numOfObjects)
//		{
//			numOfObjects = ccIt.Get();
//		}
//	}	
//
//	std::cout << "number of objects: " << numOfObjects << std::endl;
//
//	std::vector<unsigned int> countVector;
//	countVector.resize(numOfObjects+1);
//
//	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorByteType;
//
//	ByteImageType::SizeType radiusByte;
//	radiusByte.Fill(1);
//
//	ByteImageType::Pointer tissueLarge = BinaryKeeper(BinaryThreshold(vmap,Tissue), "Size", 1);
//	Write(tissueLarge,"tissueLarge.nii");
//
//	NeighborhoodIteratorByteType tIt(radiusByte,tissueLarge,region);
//
//	for (tIt.GoToBegin(), ccIt.GoToBegin(); !tIt.IsAtEnd(); ++tIt, ++ccIt)
//	{
//		if ( ccIt.Get() > 0 )
//		{
//			bool tissueNeighbor = false;
//
//			for (int i=0; i<tIt.Size(); i++)
//			{
//				if ( i != (tIt.Size()-1)/2 )
//				{
//					if ( tIt.GetPixel(i) != 0 )
//					{
//						tissueNeighbor = true;
//						break;
//					}
//				}
//			}
//
//			if (tissueNeighbor)
//				countVector[ ccIt.Get() ]++;
//		}
//	}
//
//	tissueLarge.~SmartPointer();
//
//	ByteIteratorType hbinIt(hbin,region);
//
//	for (ccIt.GoToBegin(), hbinIt.GoToBegin(); !ccIt.IsAtEnd(); ++ccIt, ++hbinIt)
//	{
//		if ( countVector[ ccIt.Get() ] < 6 )
//		{
//			hbinIt.Set( 0 );	
//		}
//	}
//
//	Write(hbin,"hbinConnected.nii");
//
//	VoxelIteratorType vmapIt(vmap,region);
//
//	for (vmapIt.GoToBegin(), hbinIt.GoToBegin(); !vmapIt.IsAtEnd(); ++vmapIt, ++hbinIt)
//	{
//		if (hbinIt.Get() != 0)
//			vmapIt.Set(Tissue);
//	}
//
//	Write(vmap,"vmapHessian.nii");
//
//	/* overlay with largest tissue component
//	hbin = BinaryOr(hbin,
//		BinaryKeeper( BinaryThreshold(vmap,Tissue) , "Size", 1) );
//	Write(hbin,"hbinWithTissue.nii");
//
//	 keep largest
//	hbin = BinaryKeeper(hbin,"Size",1);
//	Write(hbin,"hbinLargest.nii");*/
//
//	//IteratorType inputIt(input,region);
//	//FloatIteratorType smaxIt(smax,region);
//	//ArrayIteratorType partialIt(partial,region);
//
//	//for (inputIt.GoToBegin(), smaxIt.GoToBegin(), partialIt.GoToBegin(), hbinIt.GoToBegin(); !inputIt.IsAtEnd();
//	//	 ++inputIt, ++smaxIt, ++partialIt, ++hbinIt)
//	//{
//	//	if (hbinIt.Get() != 0)
//	//	{
//	//		float I = (float) inputIt.Get();
//	//		float S = (float) smaxIt.Get();
//
//	//		ArrayType p;
//
//	//		if ( S > 0 )
//	//		{
//	//			p[0] = 0;
//
//	//			p[1] = 1 - (I/S);
//
//	//			// bounds check
//	//			p[1] = (p[1] > 1) ? 1 : p[1];
//	//			p[1] = (p[1] < 0) ? 0 : p[1];
//
//	//			p[2] = 1 - p[1];
//
//	//		} else {
//	//			p[0] = 0;
//	//			p[1] = 1;
//	//			p[2] = 0;
//	//		}
//
//	//		partialIt.Set( p );
//	//	}
//	//}
//
//	//Write(partial,"hessianPartial.nii");	
//}