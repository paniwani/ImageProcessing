//void MeasureObjectness(ImageType::Pointer input)
//{
//	// test over sigma
//	// test scale by largest eigenvalue
//
//	ImageType::RegionType region = input->GetLargestPossibleRegion();
//	ImageType::SpacingType spacing = input->GetSpacing();
//
//	float sigma=0;
//
//	// Compute smoothed Hessian
//	HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
//	hessianFilter->SetInput(input);
//	hessianFilter->SetNormalizeAcrossScale(true);
//	hessianFilter->SetSigma(sigma);
//	hessianFilter->Update();
//
//	// Get object measures at each sigma
//	for (int i=0; i<2; i++)
//	{
//		sigma = spacing[0]*(i+1);
//
//		for (int m=0; m<3; m++) // object type
//		{
//			ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
//			objectnessFilter->SetInput(hessianFilter->GetOutput());
//			objectnessFilter->SetBrightObject(true);
//			objectnessFilter->SetObjectDimension(m);
//			objectnessFilter->SetScaleObjectnessMeasure(true);
//
//			objectnessFilter->Update();
//
//			std::stringstream ss;
//			ss << "sigma_" << sigma << "_";
//
//			switch (m)
//			{
//				case 0:
//					ss << "blob";
//					break;
//				case 1:
//					ss << "vessel";
//					break;
//				case 2:
//					ss << "plate";
//					break;
//			}
//
//			ss << "_scaled";
//
//			ss << ".hdr";
//
//			WriteITK(objectnessFilter->GetOutput(),ss.str());
//
//		}	
//	}
//}

/*
void ComputeSatoHessian(ImageType::Pointer input_aniso)
{
	// Make image isotropic
	ImageType::Pointer input = ResampleImage(input_aniso);
	WriteITK(input,"input_isotropic.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();
	ImageType::SizeType size = region.GetSize();

	// Only look in tagged regions
	ByteImageType::Pointer tagged = AllocateNewByteImage(region);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeByteWithIndex tagged_iter(tagged,region);

	for (input_iter.GoToBegin(), tagged_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++tagged_iter)
	{
		if (input_iter.Get() > 180)
		{
			tagged_iter.Set(1);
		} else {
			tagged_iter.Set(0);
		}
	}

	WriteITK(tagged,"tagged_200.hdr");

	// Create response image
	ImageType::Pointer S = AllocateNewImage(region);
	
	// Compute hessian across sigma scales
	double sigma[5] = {spacing[0], 2*spacing[0], 3*spacing[0], 4*spacing[0], 5*spacing[0]};

	float alpha = 0.25;
	float gamma = 1;

	// Save eigenvals
	std::vector< ImageType::Pointer > LambdaImageVector(3);
	std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

	for (int i=0; i<3; i++)
	{
		LambdaImageVector[i] = AllocateNewImage(region);
		LambdaImageVector[i]->FillBuffer(0);
		LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
		LambdaIterVector[i].GoToBegin();
	}

	for (int k=0; k<1; k++)
	{
		for (int i=0; i<3; i++) { LambdaIterVector[i].GoToBegin(); }
		
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		//for (int j=0; j<3; j++) // loop over type: line, blob, sheet
		//{
			S->FillBuffer(0.0);
			IteratorTypeFloat4WithIndex S_iter(S,region);

			hessian_iter.GoToBegin();
			tagged_iter.GoToBegin();
			S_iter.GoToBegin();
			
			int count=0;

			while (!hessian_iter.IsAtEnd()) 
			{

				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByValueDesc);

				for (int i=0; i<3; i++)
				{
					LambdaIterVector[i].Set( lambda[i] );
				}

				if (tagged_iter.Get() == 1)	// only tagged regions
				{
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

					val = S_sheet_dark(l, alpha, gamma);

					S_iter.Set(val);

				}

				count++;

				++hessian_iter;
				++S_iter;
				++tagged_iter;

				for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }

			}

			std::string type;
			type = "dark_sheet";
			//switch (j)
			//{
			//case 0: type="line";break;
			//case 1: type="blob";break;
			//case 2: type="sheet";break;
			//}

			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_type_normalizeTRUE" << type << ".hdr";

			WriteITK(S,ss.str());

			// Write eigenvals
			for (int i=0; i<3; i++)
			{		
				std::stringstream ss;
				ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
				//WriteITK(LambdaImageVector[i], ss.str());
			}
		//}
	}
}
*/

void ComputeFrangiHessian(ImageType::Pointer input_aniso)
{
	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	WriteITK(input,"input_isotropic.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	ImageType::Pointer output = AllocateNewImage(region);
	output->FillBuffer(0);
	output->SetSpacing(iso_spacing);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);
	
	// Compute hessian across sigma scales
	double sigma[2] = {0.56, 1.4};

	for (int k=0; k<2; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();

		ObjectnessFilterType::Pointer ofilter = ObjectnessFilterType::New();
		ofilter->SetInput(hessianFilter->GetOutput());
		ofilter->SetBrightObject(false);
		ofilter->SetObjectDimension(2); //plate
		ofilter->Update();

		IteratorTypeFloat4WithIndex object_iter(ofilter->GetOutput(),region);

		for (input_iter.GoToBegin(), object_iter.GoToBegin(), output_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++object_iter, ++output_iter)
		{
			if (input_iter.Get() > 200)
			{
				output_iter.Set( object_iter.Get() );
			}
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_dark_plate.hdr";
		WriteITK( ResampleImage( output, aniso_spacing ) ,ss.str());
	}
}

ImageType::Pointer HessianResponse(ImageType::Pointer input_aniso, double alpha, double beta, double gamma, double eta)
{	
	// Set intensity and lambda3 thresholds
	float intensity_threshold_lower = 150;
	float intensity_threshold_upper = 1000;
	float lambda3_threshold = 0;

	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	WriteITK(input,"input_isotropic.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	ImageType::Pointer output = AllocateNewImage(region);
	output->SetSpacing(iso_spacing);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);

	// Compute hessian across sigma scales
	double sigma[2] = {0.56, 1.4};

	for (int k=0; k<1; k++)
	{
		// Reset buffer to see each scale independently
		output->FillBuffer(0);

		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		
		//// Create helper images and iters to visualize algorithm
		std::vector< ImageType::Pointer > LambdaImageVector(3);
		std::vector< ImageType::Pointer > FImageVector(6); // A, B, C1, C2, Cup, Rut

		std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);
		std::vector< itk::ImageRegionIterator<ImageType> > FIterVector(6);

		ImageType::Pointer eigen_norm = AllocateNewImage(region);
		eigen_norm->SetSpacing( iso_spacing );
		eigen_norm->FillBuffer(0);
		IteratorTypeFloat4WithIndex eigen_norm_iter(eigen_norm,region);

		for (int i=0; i<3; i++)
		{
			LambdaImageVector[i] = AllocateNewImage(region);
			LambdaImageVector[i]->SetSpacing( iso_spacing );
			LambdaImageVector[i]->FillBuffer(0);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		for (int i=0; i<6; i++)
		{
			FImageVector[i] = AllocateNewImage(region);
			FImageVector[i]->SetSpacing( iso_spacing );
			FImageVector[i]->FillBuffer(0);
			FIterVector[i] = itk::ImageRegionIterator<ImageType>( FImageVector[i], region );
			FIterVector[i].GoToBegin();
		}

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output_iter.GoToBegin();
		eigen_norm_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			// Get eigenvalues
			EigenValueArrayType lambda;

			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

			double Lambda1 = lambda[0];
			double Lambda2 = lambda[1];
			double Lambda3 = lambda[2];

			// Set submerged and non-submerged conditions
			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper && Lambda3 >= lambda3_threshold ) 
			{				

				eigen_norm_iter.Set( sqrt( Lambda1*Lambda1 + Lambda2*Lambda2 + Lambda3*Lambda3 ) );

				double Lambda3Abs = vnl_math_abs( Lambda3 );

				for (int i=0; i<3; i++)
				{
					LambdaIterVector[i].Set( vnl_math_abs( lambda[i] ) );
				}

				// A, B, C1, C2, Cup, Rut
				FIterVector[0].Set( fA(lambda, alpha) );
				FIterVector[1].Set( fB(lambda, beta, gamma) );
				FIterVector[2].Set( Lambda3Abs * fC(lambda[0], lambda[1], eta) );
				FIterVector[3].Set( Lambda3Abs * fC(lambda[1], lambda[2], eta) );
				FIterVector[4].Set( Lambda3Abs * fCup(lambda, eta) );
				FIterVector[5].Set( Lambda3Abs * fRut(lambda, alpha, beta, gamma) );
				
				double val = vnl_math_max( Lambda3Abs * fRut(lambda, alpha, beta, gamma), Lambda3Abs * fCup(lambda, eta) );

				if ( val > output_iter.Get() )
				{
					output_iter.Set( val );		
				}
			}

			if (count++ % 500000 == 0)
			{
				std::cout << count << std::endl;
			}	

			++hessian_iter;
			++input_iter;
			++output_iter;
			++eigen_norm_iter;

			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
			for (int i=0; i<6; i++) { ++FIterVector[i]; }
		}

		

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_abs_" << i << ".hdr";
			WriteITK(  LambdaImageVector[i] , ss.str());
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_eigen_norm.hdr";
		WriteITK( eigen_norm, ss.str() );

		for (int i=0; i<6; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_F_";
			
			// A, B, C1, C2, Cup, Rut
			switch (i)
			{
				case 0: ss<<"a_alpha_"<<alpha; break;
				case 1: ss<<"b_beta_"<<beta<<"_gamma_"<<gamma; break;
				case 2: ss<<"c1_eta_"<<eta; break;
				case 3: ss<<"c2_eta_"<<eta; break;
				case 4: ss<<"cup_eta_"<<eta; break;
				case 5: ss<<"rut_alpha_"<<alpha<<"_beta_"<<beta; break;
			}

			ss << ".hdr";
			
			WriteITK( FImageVector[i], ss.str());
		}
	}

	//WriteITK(output, "Hessian_iso.hdr");

	ImageType::Pointer output_aniso = ResampleImage(output, aniso_spacing);
	
	std::stringstream ss;
	ss << "Hessian.hdr";
	//WriteITK(output_aniso, ss.str());

	return output_aniso;
}

ImageType::Pointer SatoResponse(ImageType::Pointer input_aniso, ByteImageType::Pointer chamfer_colon, double alpha, double gamma)
{	
	// Set intensity and lambda thresholds
	float intensity_threshold_lower = 150;
	float intensity_threshold_upper = 1000;
	float lambda1_threshold = 0;

	//// Make image isotropic
	//ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	//
	//ImageType::SpacingType iso_spacing = aniso_spacing;
	//iso_spacing[1] = iso_spacing[0];
	//iso_spacing[2] = iso_spacing[0];

	//ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	//WriteITK(input,"input_isotropic.hdr");

	ImageType::Pointer input = input_aniso;

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	ImageType::Pointer output = AllocateNewImage(region);
	output->SetSpacing( input->GetSpacing() );

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);

	// Compute hessian across 2 sigma scales
	double sigma[2] = {0.56,1.4};

	for (int k=0; k<2; k++)
	{
		// Reset output buffer
		//output->FillBuffer(0);

		//// Keep track of weight functions
		//ImageType::Pointer w1 = AllocateNewImage(region);
		//ImageType::Pointer w2 = AllocateNewImage(region);

		//w1->SetSpacing( iso_spacing );
		//w2->SetSpacing( iso_spacing );
		//
		//w1->FillBuffer(0);
		//w2->FillBuffer(0);

		//IteratorTypeFloat4WithIndex w1_iter(w1,region);
		//IteratorTypeFloat4WithIndex w2_iter(w2,region);

		//w1_iter.GoToBegin();
		//w2_iter.GoToBegin();

		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);
		
		//// Create helper images and iters to visualize algorithm
		//std::vector< ImageType::Pointer > LambdaImageVector(3);
		//std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

		//for (int i=0; i<3; i++)
		//{
		//	LambdaImageVector[i] = AllocateNewImage(region);
		//	LambdaImageVector[i]->SetSpacing( iso_spacing );
		//	LambdaImageVector[i]->FillBuffer(0);
		//	LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
		//	LambdaIterVector[i].GoToBegin();
		//}

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output_iter.GoToBegin();
		chamfer_colon_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{
			// Set submerged and non-submerged conditions

			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper && chamfer_colon_iter.Get() == 1)
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
			/*++w1_iter;
			++w2_iter;*/

			/*for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }*/
		}

		/*for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << ".hdr";
			WriteITK(  LambdaImageVector[i], ss.str());
		}

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_w1.hdr";
		WriteITK( w1, ss.str() );

		ss.str("");
		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_w2.hdr";
		WriteITK( w2, ss.str() );*/

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_alpha_" << alpha << "_gamma_" << gamma << "_sato_output.hdr";
		WriteITK(output,ss.str());
		
		/*RescaleIntensityFilterType::Pointer rescaler = RescaleIntensityFilterType::New();
		rescaler->SetInput(output);
		rescaler->SetOutputMaximum(1);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();

		WriteITK( rescaler->GetOutput(), ss.str() );*/

	}

	//ProcessHessian(output);

	//ImageType::Pointer output_aniso = ResampleImage(output, aniso_spacing);
	
	//std::stringstream ss;
	//ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".hdr";
	////WriteITK(output_aniso, ss.str());

	return output;
}

ImageType::Pointer SatoModifiedResponse(ImageType::Pointer input_aniso, ByteImageType::Pointer chamfer_colon, double sigma, double gamma)
{	
	// Set intensity and lambda thresholds
	float intensity_threshold_lower = 150;
	float intensity_threshold_upper = 1000;
	float lambda1_threshold = 0;

	//// Make image isotropic
	//ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	//
	//ImageType::SpacingType iso_spacing = aniso_spacing;
	//iso_spacing[1] = iso_spacing[0];
	//iso_spacing[2] = iso_spacing[0];

	//ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	//WriteITK(input,"input_isotropic.hdr");

	ImageType::Pointer input = input_aniso;

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Create response image
	ImageType::Pointer output = AllocateNewImage(region);
	output->FillBuffer(0);
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation(input);

	// Create eigen norm image
	ImageType::Pointer eigen_norm = AllocateNewImage(region);
	eigen_norm->FillBuffer(0);
	eigen_norm->CopyInformation(input);

	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);
	IteratorTypeFloat4WithIndex eigen_norm_iter(eigen_norm,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);

	// Compute smoothed Hessian
	HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
	hessianFilter->SetInput(input);
	hessianFilter->SetNormalizeAcrossScale(true);
	hessianFilter->SetSigma( sigma );
	hessianFilter->Update();
	itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);
	
	// Create helper images and iters to visualize algorithm
	std::vector< ImageType::Pointer > LambdaImageVector(3);
	std::vector< itk::ImageRegionIterator<ImageType> > LambdaIterVector(3);

	for (int i=0; i<3; i++)
	{
		LambdaImageVector[i] = AllocateNewImage(region);
		LambdaImageVector[i]->SetSpacing( input->GetSpacing() );
		LambdaImageVector[i]->FillBuffer(0);
		LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
		LambdaIterVector[i].GoToBegin();
	}

	// Find max Lambda3
	double Lambda3_max = 0;

	hessian_iter.GoToBegin();
	input_iter.GoToBegin();
	chamfer_colon_iter.GoToBegin();

	while (!hessian_iter.IsAtEnd()) 
	{
		if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper && chamfer_colon_iter.Get() == 1)
		{				
			EigenValueArrayType lambda;

			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

			if ( lambda[2] > Lambda3_max )
			{
				Lambda3_max = lambda[2];
			}
		}

		++hessian_iter;
		++input_iter;
		++chamfer_colon_iter;
	}


	hessian_iter.GoToBegin();
	input_iter.GoToBegin();
	output_iter.GoToBegin();
	eigen_norm_iter.GoToBegin();
	chamfer_colon_iter.GoToBegin();

	int count=0;

	while (!hessian_iter.IsAtEnd()) 
	{
		if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper && chamfer_colon_iter.Get() == 1)
		{				
			EigenValueArrayType lambda;

			hessian_iter.Get().ComputeEigenValues(lambda);
			std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

			double Lambda1 = lambda[0];
			double Lambda2 = lambda[1];
			double Lambda3 = lambda[2];

			for (int i=0; i<3; i++)
			{
				LambdaIterVector[i].Set( lambda[i] );
			}

			eigen_norm_iter.Set( sqrt( Lambda1*Lambda1 + Lambda2*Lambda2 + Lambda3*Lambda3 ) );
			
			if ( Lambda3 >= 0)
			{

				double val = ( Lambda3 / Lambda3_max ) *pow( 1 - vnl_math_abs( Lambda1 / Lambda3), gamma) * pow( 1 - vnl_math_abs( Lambda2 / Lambda3), gamma );

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
		++eigen_norm_iter;
		++chamfer_colon_iter;

		for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }
	}

	for (int i=0; i<3; i++)
	{
		std::stringstream ss;
		ss << "sigma_" << sigma << "_lambda_" << i+1 << ".hdr";
		WriteITK(  LambdaImageVector[i], ss.str());
	}

	std::stringstream ss;
	ss << "sigma_" << sigma << "_eigen_norm.hdr";
	WriteITK(eigen_norm,ss.str());

	ss.str("");
	ss << "sigma_" << sigma << "_gamma_" << gamma << "_sato_output.hdr";
	WriteITK(output,ss.str());
		
	/*RescaleIntensityFilterType::Pointer rescaler = RescaleIntensityFilterType::New();
	rescaler->SetInput(output);
	rescaler->SetOutputMaximum(1);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();

	WriteITK( rescaler->GetOutput(), ss.str() );*/

	//ProcessHessian(output);

	//ImageType::Pointer output_aniso = ResampleImage(output, aniso_spacing);
	
	//std::stringstream ss;
	//ss << "Sato_Hessian_alpha_" << alpha << "_gamma_" << gamma  << ".hdr";
	////WriteITK(output_aniso, ss.str());

	//return output_aniso;

	return output;
}

void ProcessHessian(ImageType::Pointer hessian, VoxelTypeImage::Pointer voxel_type)
{
	// Clean hessian image using size, connectedness, etc

	WriteITK(hessian,"hessian_pre_process.hdr");

	ImageType::RegionType region = hessian->GetLargestPossibleRegion();

	// Diffuse image with coherence-enhancing diffusion
	typedef itk::AnisotropicCoherenceEnhancingDiffusionImageFilter<ImageType,ImageType> DiffusionFilterType;
	DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
	diffusionFilter->SetInput( hessian );
	diffusionFilter->Update();
	hessian = diffusionFilter->GetOutput();

	WriteITK(hessian,"hessian_coherence_diffusion.hdr");

	// Rescale hessian to [0,1]
	RescaleIntensityFilterType::Pointer rescaler = RescaleIntensityFilterType::New();
	rescaler->SetInput(hessian);
	rescaler->SetOutputMaximum(1);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	hessian=rescaler->GetOutput();
	WriteITK(hessian,"hessian_rescaled.hdr");

	// Threshold to create binary hessian mask
	ByteImageType::Pointer hessian_mask = AllocateNewByteImage(region);

	IteratorTypeByteWithIndex hessian_mask_iter(hessian_mask,region);
	IteratorTypeFloat4WithIndex hessian_iter(hessian,region);
	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);
	
	for (hessian_iter.GoToBegin(), hessian_mask_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !hessian_iter.IsAtEnd(); ++hessian_iter, ++hessian_mask_iter, ++voxel_type_iter)
	{
		if (hessian_iter.Get() > 0.25 || voxel_type_iter.Get() == Tissue)
		{
			hessian_mask_iter.Set(1);
		} else {
			hessian_mask_iter.Set(0);
		}
	}

	WriteITK(hessian_mask,"hessian_threshold_mask.hdr");

	// Run connected component
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput(hessian_mask);
	ccFilter->FullyConnectedOn();
	ccFilter->Update();
	IntImageType::Pointer label_image = ccFilter->GetOutput();
	ccFilter.~SmartPointer();

	WriteITK(label_image,"hessian_cc.hdr");

	// Keep only the largest component (tissue) by size
	typedef itk::LabelShapeKeepNObjectsImageFilter< IntImageType > LabelShapeKeepNObjectsImageFilterType;
	LabelShapeKeepNObjectsImageFilterType::Pointer labelFilter = LabelShapeKeepNObjectsImageFilterType::New();
	labelFilter->SetAttribute("Size");
	labelFilter->SetBackgroundValue(0);
	labelFilter->SetNumberOfObjects(1);
	labelFilter->SetInput( label_image );
	labelFilter->Update();
	label_image = labelFilter->GetOutput();

	WriteITK(label_image,"hessian_label_size.hdr");

	// Update voxel edge with hessian information
	IteratorTypeIntWithIndex label_image_iter(label_image,region);

	for (voxel_type_iter.GoToBegin(), label_image_iter.GoToBegin(); !voxel_type_iter.IsAtEnd(); ++voxel_type_iter, ++label_image_iter)
	{
		if ( label_image_iter.Get() == 1 )
		{
			voxel_type_iter.Set( Tissue );
		}
	}

	WriteITK(voxel_type, "voxel_type_hessian.hdr");

	//// Remove objects based on size
	//typedef itk::LabelShapeOpeningImageFilter< IntImageType > LabelShapeOpeningImageFilterType;
	//LabelShapeOpeningImageFilterType::Pointer labelFilter = LabelShapeOpeningImageFilterType::New();

	//labelFilter->SetInput( label_image );
	//labelFilter->SetAttribute("Size");
	//labelFilter->SetBackgroundValue(0);
	//labelFilter->SetLambda(1000);
	////labelFilter->SetReverseOrdering(true);
	//labelFilter->Update();
	//label_image=labelFilter->GetOutput();

	//WriteITK(label_image,"hessian_label_size.hdr");


	//LabelImageToShapeLabelMapFilterType::Pointer converter = LabelImageToShapeLabelMapFilterType::New();
	//converter->SetInput( cc );
	//converter->Update();
	//LabelMapType::Pointer labelMap = converter->GetOutput();
	//converter.~SmartPointer();

	//std::vector<LabelObjectType::Pointer> labelVector = labelMap->GetLabelObjects();
	//sort(labelVector.begin(), labelVector.end(), compareSizeOnBorder);
	//int externalAirLabel = labelVector[0]->GetLabel();
	//
	////labelVector.erase( labelVector.begin() );
	////sort(labelVector.begin(), labelVector.end(), compareSize);
	////int colonLabel = labelVector[0]->GetLabel();

	//labelVector.clear();
	//labelMap.~SmartPointer();

	//IteratorTypeIntWithIndex cc_iter( cc, fullRegion );
	//chamfer_colon_iter = IteratorTypeByteWithIndex(chamfer_colon,fullRegion);
	//
	//for (chamfer_colon_iter.GoToBegin(), cc_iter.GoToBegin();
 //       !chamfer_colon_iter.IsAtEnd() && !cc_iter.IsAtEnd();  
 //       ++chamfer_colon_iter, ++cc_iter) 
	//{
	//	if (cc_iter.Get() == externalAirLabel)
	//	{
	//		chamfer_colon_iter.Set(0);
	//	}
	//}

	//cc.~SmartPointer();


}



void EnhanceVoxelType(ImageType::Pointer input, ImageType::Pointer hessian, VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	
	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeVoxelType voxel_type_iter(voxel_type,region);
	IteratorTypeByteWithIndex chamfer_colon_iter(chamfer_colon,region);
	IteratorTypeFloat4WithIndex hessian_iter(hessian,region);
	
	/*
	// Threshold hessian based on max value
	float hmax=0;
	for (hessian_iter.GoToBegin(); !hessian_iter.IsAtEnd(); ++hessian_iter)
	{
		if ( hessian_iter.Get() > hmax )
		{
			hmax = hessian_iter.Get();
		}
	}
	
	for (hessian_iter.GoToBegin(); !hessian_iter.IsAtEnd(); ++hessian_iter)
	{
		if ( hessian_iter.Get() < 0.25*hmax )
		{
			hessian_iter.Set( 0 );
		}
	}
	

	// Apply gaussian filter to hessian
	DiscreteGaussianFilterType::Pointer dgFilter = DiscreteGaussianFilterType::New();
	dgFilter->SetInput(hessian);

	std::cout << "Smoothing hessian by " << hessian->GetSpacing()[0] << "mm" << std::endl;
	dgFilter->SetVariance(hessian->GetSpacing()[0]*hessian->GetSpacing()[0]);
	dgFilter->Update();
	hessian=dgFilter->GetOutput();
	hessian_iter=IteratorTypeFloat4WithIndex(hessian,region);

	WriteITK(hessian,"hessian_smoothed.hdr");

	*/

	// Make binary hessian mask
	ByteImageType::Pointer hessian_mask = AllocateNewByteImage(region);
	hessian_mask->FillBuffer(0);
	IteratorTypeByteWithIndex hessian_mask_iter(hessian_mask,region);
	
	for (chamfer_colon_iter.GoToBegin(), hessian_iter.GoToBegin(), hessian_mask_iter.GoToBegin(), voxel_type_iter.GoToBegin(); !chamfer_colon_iter.IsAtEnd(); ++chamfer_colon_iter, ++hessian_iter, ++hessian_mask_iter, ++voxel_type_iter)
	{
		if (chamfer_colon_iter.Get() == 1 && ( hessian_iter.Get() > 0 || voxel_type_iter.Get() == Tissue ) )
		{
			hessian_mask_iter.Set(1);
		}
	}

	WriteITK(hessian_mask, "hessian_tissue_mask.hdr");

	// Run connected component
	ConnectedComponentFilterType::Pointer ccFilter = ConnectedComponentFilterType::New();
	ccFilter->SetInput(hessian_mask);
	ccFilter->Update();
	IntImageType::Pointer cc = ccFilter->GetOutput();
	IteratorTypeIntWithIndex cc_iter(cc,region);
	
	ccFilter.~SmartPointer();
	WriteITK(cc,"hessian_tissue_mask_cc.hdr");

	// Convert to Label map
	LabelImageToShapeLabelMapFilterType::Pointer converter = LabelImageToShapeLabelMapFilterType::New();
	converter->SetInput(cc);
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();
	converter.~SmartPointer();

	// Sort based on size
	std::vector<LabelObjectType::Pointer> labelVector = labelMap->GetLabelObjects();
	sort(labelVector.begin(), labelVector.end(), compareSize);

	int tissueLabel = labelVector[0]->GetLabel();

	for (hessian_iter.GoToBegin(), cc_iter.GoToBegin(); !hessian_iter.IsAtEnd(); ++hessian_iter, ++cc_iter)
	{
		if ( cc_iter.Get() != tissueLabel )
		{
			hessian_iter.Set( 0 );
		}
	}

	WriteITK(hessian, "hessian_post_cc.hdr");

	/*

	// Find small sized components
	



	std::vector<unsigned long> smallLabelVector;
	
	for (unsigned int i=0; i < labelVector.size() ; i++)
	{
		if (labelVector[i]->GetSize() < 5)
		{
			smallLabelVector.push_back( labelVector[i]->GetLabel() );
		}
	}
	
	// Delete small labels from hessian image
	for (cc_iter.GoToBegin(), hessian_iter.GoToBegin(); !cc_iter.IsAtEnd(); ++cc_iter, ++hessian_iter)
	{
		if (cc_iter.Get() > 0)
		{
			bool isSmallComp = false;
			for (unsigned int i=0; i < smallLabelVector.size() ; i++)
			{
				if (cc_iter.Get() == smallLabelVector[i])
				{
					isSmallComp = true;
					break;
				}
			}
			
			if (isSmallComp)
			{
				cc_iter.Set(0);
				hessian_iter.Set(0);
			}
		}
	}

	

	WriteITK(hessian, "hessian_post_cc.hdr");
	
	
	for (voxel_type_iter.GoToBegin(), hessian_iter.GoToBegin(), chamfer_colon_iter.GoToBegin(); !voxel_type_iter.IsAtEnd(); ++voxel_type_iter, ++hessian_iter, ++chamfer_colon_iter)
	{
		if (chamfer_colon_iter.Get() == 1 && ( voxel_type_iter.Get() == Unclassified || voxel_type_iter.Get() == Stool) && hessian_iter.Get() > 0)
		{
			voxel_type_iter.Set( Tissue );
		}
	}
	
	WriteITK(voxel_type, "hessian_enhanced.hdr");
	*/
}

/*
ImageType::Pointer SatoHessianEdgeEnhancingDiffusion(ImageType::Pointer input_aniso)
{
	float intensity_threshold_lower = 150;
	float intensity_threshold_upper = 700;
	float lambda1_threshold = 100;

	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	WriteITK(input,"input_isotropic.hdr");

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	ImageType::Pointer output = AllocateNewImage(region);
	output->FillBuffer(0);
	
	IteratorTypeFloat4WithIndex input_iter(input,region);
	IteratorTypeFloat4WithIndex output_iter(output,region);

	// Compute response across scales
	float sigma[2] = {0.56, 1.4};

	for (int k=0; k<2; k++)
	{
		// Diffusion filtering
		typedef itk::AnisotropicEdgeEnhancementDiffusionImageFilter<ImageType,ImageType> DiffusionFilterType;
		DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
		diffusionFilter->SetInput( input );
		diffusionFilter->SetSigma( sigma[k] );
		diffusionFilter->Update();
		ImageType::Pointer diffusion = diffusionFilter->GetOutput();
		IteratorTypeFloat4WithIndex diffusion_iter(diffusion,region);

		std::stringstream ss;
		ss << "diffusion_" << sigma[k] << ".hdr";
		WriteITK(diffusion,ss.str());

		// Compute hessian on diffused input with NO gaussian smoothing (variance = 0)
		typedef itk::DiscreteHessianGaussianImageFunction<ImageType> HessianFunctionType;
		HessianFunctionType::Pointer hessianFunction = HessianFunctionType::New();
		hessianFunction->SetInputImage( diffusion );
		hessianFunction->NormalizeAcrossScaleOn();
		hessianFunction->Initialize();

		for (input_iter.GoToBegin(), output_iter.GoToBegin(), diffusion_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter, ++diffusion_iter)
		{
			if ( input_iter.Get() >= intensity_threshold_lower && input_iter.Get() <= intensity_threshold_upper )
			{
				// Get hessian matrix
				SymmetricSecondRankTensorType hessianMatrix = hessianFunction->EvaluateAtIndex( input_iter.GetIndex() );
				
				// Get eigenvalues and sort
				SymmetricSecondRankTensorType::EigenValuesArrayType lambda;
				hessianMatrix.ComputeEigenValues( lambda );

				std::sort(lambda.Begin(), lambda.End(), OrderByValueDesc);

				double Lambda1 = lambda[0];
				double Lambda2 = lambda[1];
				double Lambda3 = lambda[2];

				if ( Lambda1 >= lambda1_threshold )
				{
					float l[3] = { lambda[0], lambda[1], lambda[2] };
					
					double alpha = 0.25;
					double gamma = 0.5;

					double val = S_sheet_dark(l, alpha, gamma);

					if ( val > output_iter.Get() )
					{
						output_iter.Set( val );
					}
				}
			}
		}
	}

	return output;
}

*/

void GetHessian(ImageType::Pointer input_aniso)
{
	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	WriteITK(input,"input_isotropic.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	IteratorTypeFloat4WithIndex input_iter(input,region);

	// Compute hessian across 2 sigma scales
	double sigma[2] = {0.56, 1.4};

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
			LambdaImageVector[i]->SetSpacing( input->GetSpacing() );
			LambdaImageVector[i]->FillBuffer(0);
			LambdaIterVector[i] = itk::ImageRegionIterator<ImageType>( LambdaImageVector[i], region);
			LambdaIterVector[i].GoToBegin();
		}

		// Allocate eigen norm image
		ImageType::Pointer eigen_norm = AllocateNewImage(region);
		eigen_norm->SetSpacing( input->GetSpacing() );
		eigen_norm->FillBuffer(0);
		IteratorTypeFloat4WithIndex eigen_norm_iter(eigen_norm,region);

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		eigen_norm_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{			
			if (input_iter.Get() >= 150 && input_iter.Get() <= 1000)
			{
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

				if ( lambda[2] >= 0 )
				{
					for (int i=0; i<3; i++)
					{
						LambdaIterVector[i].Set( lambda[i] );
					}

					eigen_norm_iter.Set( sqrt( lambda[0]*lambda[0] + lambda[1]*lambda[1] + lambda[2]*lambda[2] ) );

				}

			}
			
			if (++count % 500000 == 0)
			{
				std::cout << count << std::endl;
			}	

			++hessian_iter;
			++input_iter;
		

			++eigen_norm_iter;
			for (int i=0; i<3; i++) { ++LambdaIterVector[i]; }


		} // end hessian loop

		for (int i=0; i<3; i++)
		{
			std::stringstream ss;
			ss << "sigma_" << sigma[k] << "_lambda_" << i << "_ORDER_MAG_REGION_input_150_1000_lambda1_0.hdr";
			WriteITK( LambdaImageVector[i], ss.str());
		}

		//// Rescale eigen norm
		//RescaleIntensityFilterType::Pointer rescaler = RescaleIntensityFilterType::New();
		//rescaler->SetInput(eigen_norm);
		//rescaler->SetOutputMaximum(1);
		//rescaler->SetOutputMinimum(0);
		//rescaler->Update();
		//eigen_norm=rescaler->GetOutput();

		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_eigen_norm_ORDER_MAG_REGION_input_150_700_lambda1_0.hdr";
		WriteITK( eigen_norm, ss.str() );
	}
}

void UnderstandHessian(ImageType::Pointer input_aniso)
{
	// Make image isotropic
	ImageType::SpacingType aniso_spacing = input_aniso->GetSpacing();
	
	ImageType::SpacingType iso_spacing = aniso_spacing;
	iso_spacing[1] = iso_spacing[0];
	iso_spacing[2] = iso_spacing[0];

	ImageType::Pointer input = ResampleImage(input_aniso, iso_spacing);
	WriteITK(input,"input_isotropic.hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Allocate output across scales
	ImageType::Pointer output = AllocateNewImage(region);
	output->FillBuffer(0);

	IteratorTypeFloat4WithIndex output_iter(output,region);
	IteratorTypeFloat4WithIndex input_iter(input,region);
	
	double gamma = 2.0;


	// Compute hessian across 2 sigma scales
	double sigma[2] = {0.56, 1.4};

	for (int k=0; k<2; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);
		
		// Response at each scale N
		ImageType::Pointer N = AllocateNewImage(region);
		N->FillBuffer(0);
		IteratorTypeFloat4WithIndex N_iter(N,region);

		/*double lambda3_max = 0;

		//Find lambda3 max
		for (hessian_iter.GoToBegin(), input_iter.GoToBegin(); !hessian_iter.IsAtEnd(); ++hessian_iter, ++input_iter)
		{
			if (input_iter.Get() >= 150 && input_iter.Get() <= 700)
			{
				EigenValueArrayType lambda;
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

				if ( lambda[2] > lambda3_max )
				{
					lambda3_max = lambda[2];
				}
			}
		}*/

		hessian_iter.GoToBegin();
		input_iter.GoToBegin();
		output_iter.GoToBegin();
		N_iter.GoToBegin();

		int count=0;

		while (!hessian_iter.IsAtEnd()) 
		{			
			if (input_iter.Get() >= 150 && input_iter.Get() <= 1000)
			{
				EigenValueArrayType lambda;

				// Get eigenvalues
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

				if ( lambda[2] >= 0 )
				{
					double Ra = abs( lambda[0] / lambda[2] );
					double Rb = abs( lambda[1] / lambda[2] );
					

					double val = pow( 1-Ra , gamma ) * pow( 1-Rb, gamma );

					N_iter.Set( val );

					if ( val > output_iter.Get() )
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
			++N_iter;
		} // end hessian loop


		std::stringstream ss;
		ss << "sigma_" << sigma[k] << "_gamma_" << gamma << "_N.hdr";
		WriteITK(N,ss.str());

	}

	std::stringstream ss;
	ss << "gamma_" << gamma << "_hessian.hdr";
	WriteITK(output,ss.str());
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