float PolyDist(float X, vnl_real_polynomial poly, float x, float y) 
{
    return sqrtf(vnl_math_sqr(X-x)+vnl_math_sqr(poly.evaluate(X)-y));
}

float PolyMinDist(vnl_real_polynomial poly, float x, float y)
{
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

float Stool_Air_ComputeSmax(float intensity[], float gradient_magnitude[], int size) {
	float intensity2[5];
	for (int i=0;i<size;i++) {
		intensity2[i]=intensity[i]+1000;
	}

	return ComputeSmax(intensity2,gradient_magnitude, size)-1000;
}

float ComputeSmaxFit(PixelType intensity[], float gradient_magnitude[], PixelType Smax) 
{
	float R=0;
	
	for (int i=0;i<5;i++) 
	{
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



ImageType::Pointer ComputeNeighborhoodSmax(ImageType::Pointer &input, VoxelImageType::Pointer &v, ByteImageType::Pointer &mask)
{
	ImageType::Pointer smax = ImageType::New();
	smax->SetSpacing( input->GetSpacing() );
	smax->SetDirection( input->GetDirection() );
	smax->SetRegions( REGION );
	smax->CopyInformation( input );
	smax->Allocate();
	smax->FillBuffer(0);

	IteratorType smax_iter(smax,REGION);
	IteratorType input_iter(input,REGION);
	ByteIteratorType mask_iter(mask,REGION);

	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorVoxelType;
	NeighborhoodIteratorVoxelType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 4;
	radius[1] = 4;

	NeighborhoodIteratorVoxelType nit(radius, v, v->GetLargestPossibleRegion());

	VoxelImageType::SizeType size = nit.GetSize();
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
			
				VoxelImageType::IndexType idx = nit.GetIndex(i);
				
				if ( REGION.IsInside( idx ) )
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

bool checkBounds(ImageType::RegionType &region, ContinuousIndexType &index)
{
	ImageType::IndexType idx1 = region.GetIndex();
	ImageType::SizeType size = region.GetSize();
	
	ImageType::IndexType idx2;

	for (int i=0; i<3; i++)
	{
		idx2[i] = idx1[i] + size[i] - 1;

		if ( index[i] < idx1[i] || index[i] > idx2[i] )
			return false;
	}

	return true;
}

FloatImageType::Pointer ComputePartials(ImageType::Pointer &input, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, ImageType::Pointer &smax)
{
	FloatImageType::Pointer partial = FloatImageType::New();
	partial->SetRegions( REGION );
	partial->SetSpacing( input->GetSpacing() );
	partial->Allocate();
	partial->FillBuffer(0);
	
	FloatIteratorType partial_iter(partial,REGION);
	IteratorType input_iter(input,REGION);
	VoxelIteratorType vmap_iter(vmap,REGION);
	ByteIteratorType colon_iter(colon,REGION);
	IteratorType smax_iter(smax,REGION);

	for (partial_iter.GoToBegin(), input_iter.GoToBegin(), vmap_iter.GoToBegin(), colon_iter.GoToBegin(), smax_iter.GoToBegin(); !partial_iter.IsAtEnd();
		++partial_iter, ++input_iter, ++vmap_iter, ++colon_iter, ++smax_iter)
	{

		if ( colon_iter.Get() == 255 )
		{
			float value = 0;

			float I = (float) input_iter.Get() + BACKGROUND;
			float S = (float) smax_iter.Get() + BACKGROUND;

			switch ( vmap_iter.Get() ) 
			{
				case Air:

					value=0;

					break;

				case Tissue:

					value=1;

					break;

				case Stool:

					value=0;

					break;

				case TissueAir:

					value=1+(I/1000);
					
					if (value <= 0) { value = 0; }
					if (value >= 1) { value = 1; }

					break;

				case TissueStool:

					float value2;
					value2=I/S;
				
					if (value2 <= 0) { value2 = 0; }
					if (value2 >= 1) { value2 = 1; }

					value = 1-value2;

					break; 

				case StoolAir:

					value=0;

					break;
			}   
			
			partial_iter.Set(value);
		}
	}

	return partial;
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
FloatImageType::Pointer QuadraticRegression(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, FloatImageType::Pointer &gradient_magnitude)
{
	ImageType::SpacingType spacing = input->GetSpacing();

	// get gradient vector
	/*
	typedef itk::DiscreteGaussianImageFilter<ImageType, FloatImageType> DiscreteGaussianImageFilterType;
	DiscreteGaussianImageFilterType::Pointer smoother = DiscreteGaussianImageFilterType::New();
	smoother->SetInput( input );
	smoother->SetVariance( spacing[0]*spacing[0] );
	smoother->SetMaximumKernelWidth( 3 );
	smoother->Update();

	Write(smoother->GetOutput(),"smoothed.nii");
	*/

	typedef itk::GradientImageFilter<ImageType> GradientImageFilterType;
	GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->SetUseImageDirection(false);
	gradientFilter->Update();
	VectorImageType::Pointer gradient = gradientFilter->GetOutput();
	VectorIteratorType gradient_iter(gradient,REGION);

	// normalize
	for (gradient_iter.GoToBegin(); !gradient_iter.IsAtEnd(); ++gradient_iter)
	{
		VectorType g = gradient_iter.Get();

		float norm = g.GetNorm();

		if ( norm > 0 )
		{
			for (int i=0; i<3; i++)
				g[i] /= norm;
		}

		// tolerance
		for (int i=0; i<3; i++)
			g[i] = abs(g[i]) > 0.001 ? g[i] : 0; 

		gradient_iter.Set(g);
	}

	// interp grad, use 2x gradient magnitude
	typedef itk::MultiplyByConstantImageFilter<FloatImageType, int, FloatImageType> MultiplyByConstantImageFilterType;
	MultiplyByConstantImageFilterType::Pointer multiplier = MultiplyByConstantImageFilterType::New();
	multiplier->SetInput( gradient_magnitude );
	multiplier->SetConstant( 2 );
	multiplier->Update();

	Write(multiplier->GetOutput(),"gradient2x.nii");

	typedef itk::BSplineInterpolateImageFunction<FloatImageType> InterpolatorFloatType;
	InterpolatorFloatType::Pointer grad_interp = InterpolatorFloatType::New();
	grad_interp->SetSplineOrder(3);
	grad_interp->SetInputImage( multiplier->GetOutput() );

	// shift input to air -1024 and cast to short
	typedef itk::CastImageFilter<ImageType,ShortImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput( input );

	typedef itk::AddConstantToImageFilter< ShortImageType, short, ShortImageType > AddConstantToImageFilterType;
	AddConstantToImageFilterType::Pointer adder = AddConstantToImageFilterType::New();
	adder->SetInput( caster->GetOutput() );
	adder->SetConstant( BACKGROUND );
	adder->Update();

	ShortImageType::Pointer input_short = adder->GetOutput();
	Write(input_short,"input_short.nii");

	// interp input
	typedef itk::BSplineInterpolateImageFunction<ShortImageType> InterpolatorShortType;
	InterpolatorShortType::Pointer input_interp = InterpolatorShortType::New();
	input_interp->SetSplineOrder(3);
	input_interp->SetInputImage( input_short );

	InterpolatorShortType::ContinuousIndexType startIndex = input_interp->GetStartIndex();
	
	std::cout << "Start QR Index: " << startIndex[0] << " " << startIndex[1] << " " << startIndex[2] << std::endl;

	InterpolatorShortType::ContinuousIndexType endIndex = input_interp->GetEndIndex();
	
	std::cout << "End QR Index: " << endIndex[0] << " " << endIndex[1] << " " << endIndex[2] << std::endl;

	// get local stool max
	ImageType::Pointer smax = ComputeNeighborhoodSmax(input,vmap,colon);
	Write(smax,"smax.nii");

	// ve
	IteratorType smax_iter(smax,REGION);
	ByteIteratorType colon_iter(colon,REGION);
	VoxelIteratorType vmap_iter(vmap,REGION);

	float d[3][5];
	d[0][0] = -1.5; d[0][1] = -1.0; d[0][2] = 0; d[0][3] = 1.0; d[0][4] = 1.5;
	d[1][0] = -1.0; d[1][1] = -0.5; d[1][2] = 0; d[1][3] = 0.5; d[1][4] = 1.0;
	d[2][0] = -0.6; d[2][1] = -0.3; d[2][2] = 0; d[2][3] = 0.3; d[2][4] = 0.6;

	int count=0;

	std::ofstream file;
	file.open("qr.txt");

	for (smax_iter.GoToBegin(), gradient_iter.GoToBegin(), colon_iter.GoToBegin(), vmap_iter.GoToBegin(); !smax_iter.IsAtEnd(); 
		++smax_iter, ++gradient_iter, ++colon_iter, ++vmap_iter)
	{

		if (colon_iter.Get() == 255 && vmap_iter.Get() == Unclassified)
		{
			if (count++ % 1000 == 0)
				std::cout << count << std::endl;

			ImageType::IndexType idx = gradient_iter.GetIndex();

			VoxelType v;

			//if (smax_iter.Get() > 0)
			//{

				float dist = 0;

				for (int i=0; i<3; i++) // iterate through each set of distances
				{
					float in[5];
					float gr[5];
					
					for (int j=0; j<5; j++) // iterate through each distance
					{
						ContinuousIndexType odx = idx;

						for (int k=0; k<3; k++)
							odx[k] += gradient_iter.Get()[k]*d[i][j];
						
						//if ( !checkBounds(REGION,odx) )
						//	odx = idx;
						
						in[j] = input_interp->EvaluateAtContinuousIndex(odx);
						gr[j] = grad_interp->EvaluateAtContinuousIndex(odx);
					}

					float distanceTS=0;
					float distanceSA=0;
					float distanceTA=0;

					float localDist=0;
					VoxelType localV=Unclassified;
					
					//distanceTS=AverageTissueStoolDist(smax_iter.Get()+BACKGROUND,in,gr);
					//distanceSA=AverageStoolAirDist(smax_iter.Get()+BACKGROUND,in,gr);

					float TS_smax = ComputeSmax(in,gr,5);
					float SA_smax = Stool_Air_ComputeSmax(in,gr,5);

					distanceTS=AverageTissueStoolDist( TS_smax, in, gr );
					distanceSA=AverageStoolAirDist( SA_smax , in, gr );
					distanceTA=AverageTissueAirDist(in,gr);

					if (distanceSA<=distanceTS && distanceSA<=distanceTA) 
					{
						localDist = distanceSA;
						localV = StoolAir;

					} else if (distanceTS<=distanceTA && distanceTS<=distanceSA) {
						localDist = distanceTS;
						localV = TissueStool;

					} else {
						localDist = distanceTA;
						localV = TissueAir;
					}

					if ( i == 0 )
					{
						dist = localDist;
						v = localV;
					} else {
						if ( localDist < dist )
						{
							dist = localDist;
							v = localV;
						}
					}

					/*std::cout << "iteration: " << i << std::endl;
					for (int k=0; k<5;k++)
						std::cout << in[k] << " ";
					std::cout << std::endl;
					for (int k=0; k<5;k++)
						std::cout << gr[k] << " ";
					std::cout << std::endl;
					*/
					/*std::cout << std::endl;
					std::cout << "Distances (TS,SA,TA)" << std::endl;
					std::cout << distanceTS << " " << distanceSA << " " << distanceTA << std::endl;*/

					ImageType::SpacingType spacing = input->GetSpacing();

					file << "\n";
					file << idx[0]*spacing[0] << " " << idx[1]*spacing[1] << " " << idx[2]*spacing[2] << "\n";
					for (int k=0; k<5; k++)
						file << in[k] << " ";
					file << "\n";

					for (int k=0; k<5; k++)
						file << gr[k] << " ";
					file << "\n";

					file << distanceTS << " " << distanceSA << " " << distanceTA << "\n";




				}
			//} else {
			//	v = TissueAir;
			//}

			vmap_iter.Set(v);

		}		
	}

	file.close();
	
	Write(vmap,"qr.nii");

	return ComputePartials(input,vmap,colon,smax);
}

void LocalBoundary(VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon)
{
	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorType;
	
	VoxelImageType::SizeType radius;
	radius.Fill(0);
	radius[0]=2;
	radius[1]=2;

	NeighborhoodIteratorType nit(radius,vmap,REGION);
	ByteIteratorType colon_iter(colon,REGION);

	int n = nit.Size();

	std::cout << n << std::endl;

	unsigned int voxelCount[3] = {0,0,0}; // air/tissue/stool

	for (nit.GoToBegin(), colon_iter.GoToBegin(); !nit.IsAtEnd(); ++nit, ++colon_iter)
	{
		if (colon_iter.Get() == 255)
		{
			VoxelType centerVoxel = nit.GetCenterPixel();

			if (centerVoxel == Unclassified)
			{
				// Reset voxel count
				for(int i=0; i<3; i++)
					voxelCount[i] = 0;

				// Count neighboring solid classes
				for (int i=0; i<n; i++)
				{
					VoxelType v = nit.GetPixel(i);

					if (v == Air)
					{
						voxelCount[0]++;
					} else if (v == Tissue) {
						voxelCount[1]++;
					} else if (v == Stool) {
						voxelCount[2]++;
					}
				}
				
				// Set boundary based on neighbors
				if (voxelCount[0] == 0 && voxelCount[1] != 0 && voxelCount[2] != 0)
				{
					nit.SetCenterPixel(TissueStool);
				} else if (voxelCount[0] != 0 && voxelCount[1] == 0 && voxelCount[2] != 0) {
					nit.SetCenterPixel(StoolAir);
				} else if (voxelCount[0] != 0 && voxelCount[1] != 0 && voxelCount[2] == 0) {
					nit.SetCenterPixel(TissueAir);
				}

				/*if (voxelCount[0] > voxelCount[2] && voxelCount[1] > voxelCount[2])
				{
					nit.SetCenterPixel(TissueAir);
				} else if (voxelCount[1] > voxelCount[0] && voxelCount[2] > voxelCount[0]) {
					nit.SetCenterPixel(TissueStool);
				} else if (voxelCount[0] > voxelCount[1] && voxelCount[2] > voxelCount[1]) {
					nit.SetCenterPixel(StoolAir);
				}*/
			}
		}
		
	}

	Write(vmap,"local_vmap.nii");
}	

