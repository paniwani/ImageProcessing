ImageType::Pointer SatoResponse(ImageType::Pointer &input, ByteImageType::Pointer &colon, double alpha, double gamma, float T0)
{	
	// Set intensity and lambda thresholds
	float intensity_threshold_lower = T0;
	float intensity_threshold_upper = 1500;
	float lambda1_threshold = 0;

	ImageType::Pointer input = input;

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	FloatImageType::Pointer output = FloatImageType::New();
	output->SetSpacing( input->GetSpacing() );
	output->SetRegions( region );
	output->CopyInformation( input );
	output->Allocate();
	output->FillBuffer(0);

	IteratorType input_iter(input,region);
	FloatIteratorType output_iter(output,region);
	ByteIteratorType colon_iter(colon,region);

	// Compute hessian across 2 sigma scales
	double sigma[2] = {0.56,1.4};

	for (int k=0; k<2; k++)
	{

		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output_iter.GoToBegin();
		colon_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			// Set submerged and non-submerged conditions

			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper && chamfer_colon_iter.Get() != 0)
			{				
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				double Lambda1 = lambda[0];
				double Lambda2 = lambda[1];
				double Lambda3 = lambda[2];
				
				
				if (Lambda1 >= lambda1_threshold)
				{
					/*for (int i=0; i<3; i++)
					{
						LambdaIterVector[i].Set( lambda[i] );
					}*/

					/*w1_iter.Set( weight_dark(Lambda1, Lambda2, alpha, gamma) );
					w2_iter.Set( weight_dark(Lambda1, Lambda3, alpha, gamma) );*/

					float eig[3] = {Lambda1, Lambda2, Lambda3};

					double val = S_sheet_dark(eig, alpha, gamma);

					if ( val > output_iter.Get() )
					{
						output_iter.Set( val );		
					}
				}


			}

			if (++count % 500000 == 0)
			{
				std::cout << count << std::endl;
			}	

			++hessian_iter;
			++input_iter;
			++output_iter;
			++chamfer_colon_iter;
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_sato_output.nii";
		Write(output,ss.str());

	}
	
	std::stringstream ss;
	ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".nii";
	Write(output, ss.str());

	return output;
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

float S_line(float lambda[3], float alpha, float gamma)
{
	if ( (lambda[2] <= lambda[1]) && (lambda[1] < 0) )
	{
		return abs(lambda[2])*omega(lambda[1],lambda[2],gamma)*weight(lambda[0],lambda[1],alpha,gamma);
	}

	return 0;
}

float S_line_dark(float lambda[3], float alpha, float gamma)
{
	if ( (lambda[0] >= lambda[1]) && (lambda[1] > 0) )
	{
		return abs(lambda[0])*omega_dark(lambda[0],lambda[1],gamma)*weight_dark(lambda[1],lambda[2],alpha,gamma);
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

float S_sheet_dark(float lambda[3], float alpha, float gamma)
{
	if (lambda[0] > 0)
	{
		return abs(lambda[0])*weight_dark(lambda[0], lambda[1], alpha, gamma)*weight_dark(lambda[0], lambda[2], alpha, gamma);
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

float thinness( float l[3] )
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