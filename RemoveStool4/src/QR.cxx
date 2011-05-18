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

	IteratorType smaxIt(smax,REGION);
	IteratorType inputIt(input,REGION);
	ByteIteratorType maskIt(mask,REGION);

	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorVoxelType;
	NeighborhoodIteratorVoxelType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 4;
	radius[1] = 4;

	NeighborhoodIteratorVoxelType nit(radius, v, v->GetLargestPossibleRegion());

	VoxelImageType::SizeType size = nit.GetSize();
	int n = size[0]*size[1]*size[2] - 1;

	nit.GoToBegin();
	maskIt.GoToBegin();
	smaxIt.GoToBegin();
	inputIt.GoToBegin();

	PixelType max=0;

	while (!nit.IsAtEnd())
	{
		if (maskIt.Get() == 255 && nit.GetCenterPixel()==Unclassified)
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

			smaxIt.Set( max );
		
		} else {
			smaxIt.Set( 0 );
		}

		++nit;
		++maskIt;
		++smaxIt;
		++inputIt;
	}

	return smax;
}

bool checkBounds(ImageType::RegionType &region, ImageType::IndexType &index)
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

ArrayImageType::Pointer ComputePartials(ImageType::Pointer &input, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, FloatImageType::Pointer &smax)
{
	ArrayImageType::Pointer partial = ArrayImageType::New();
	partial->SetRegions( REGION );
	partial->SetSpacing( input->GetSpacing() );
	partial->Allocate();
	
	float zeros[3] = {0,0,0};
	partial->FillBuffer(ArrayType(zeros));
	
	ArrayIteratorType partialIt(partial,REGION);

	IteratorType inputIt(input,REGION);
	VoxelIteratorType vmapIt(vmap,REGION);
	ByteIteratorType colonIt(colon,REGION);
	FloatIteratorType smaxIt(smax,REGION);

	for (partialIt.GoToBegin(), inputIt.GoToBegin(), vmapIt.GoToBegin(), colonIt.GoToBegin(), smaxIt.GoToBegin(); !partialIt.IsAtEnd();
		++partialIt, ++inputIt, ++vmapIt, ++colonIt, ++smaxIt)
	{

		if ( colonIt.Get() == 255 )
		{
			float value[3] = {0,0,0};

			float I = (float) inputIt.Get();
			float S = (float) smaxIt.Get();

			switch ( vmapIt.Get() ) 
			{
				case Air:

					value[0]=1;					
					value[1]=0;
					value[2]=0;

					break;

				case Tissue:

					value[0]=0;
					value[1]=1;					
					value[2]=0;

					break;

				case Stool:

					value[0]=0;
					value[1]=0;
					value[2]=1;

					break;

				case TissueAir:

					value[1]=1+(I/1000);

					if (value[1] <= 0) { value[1] = 0; }
					if (value[1] >= 1) { value[1] = 1; }

					value[0]=1-value[1];
					value[2]=0;

					break;

				case TissueStool:

					value[0]=0;
					value[2]=I/S;

					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[1]=1-value[2];

					break; 

				case StoolAir:

					value[2]=(I+1000)/(S+1000);

					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[0]=1-value[2];
					value[1]=0;

					break;
			}   

			ArrayType data(value);
			partialIt.Set( data );
		}
	}

	Write(partial,"partial.nii");

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
ArrayImageType::Pointer QuadraticRegression(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, FloatImageType::Pointer &gradientMagnitude, PixelType tissueStoolThreshold)
{
	ImageType::SpacingType spacing = input->GetSpacing();

	typedef itk::GradientImageFilter<ImageType> GradientImageFilterType;
	GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->SetUseImageDirection(false);
	gradientFilter->Update();
	VectorImageType::Pointer gradient = gradientFilter->GetOutput();
	VectorIteratorType gradientIt(gradient,REGION);

	// normalize
	for (gradientIt.GoToBegin(); !gradientIt.IsAtEnd(); ++gradientIt)
	{
		VectorType g = gradientIt.Get();

		float norm = g.GetNorm();

		if ( norm > 0 )
		{
			for (int i=0; i<3; i++)
				g[i] /= norm;
		}

		// tolerance
		for (int i=0; i<3; i++)
			g[i] = abs(g[i]) > 0.001 ? g[i] : 0; 

		gradientIt.Set(g);
	}

	// interp grad, use 2x gradient magnitude
	typedef itk::MultiplyByConstantImageFilter<FloatImageType, int, FloatImageType> MultiplyByConstantImageFilterType;
	MultiplyByConstantImageFilterType::Pointer multiplier = MultiplyByConstantImageFilterType::New();
	multiplier->SetInput( gradientMagnitude );
	multiplier->SetConstant( 2 );
	multiplier->Update();
	Write(multiplier->GetOutput(),"gradient2x.nii");

	typedef itk::BSplineInterpolateImageFunction<FloatImageType> InterpolatorFloatType;
	InterpolatorFloatType::Pointer gradInterp = InterpolatorFloatType::New();
	gradInterp->SetSplineOrder(3);
	gradInterp->SetInputImage( multiplier->GetOutput() );

	// interp input
	typedef itk::BSplineInterpolateImageFunction<ImageType> InterpolatorType;
	InterpolatorType::Pointer inputInterp = InterpolatorType::New();
	inputInterp->SetSplineOrder(3);
	inputInterp->SetInputImage( input );

	// assign voxel boundaries to transition types
	ByteIteratorType colonIt(colon,REGION);
	VoxelIteratorType vmapIt(vmap,REGION);

	float d[3][5];
	d[0][0] = -1.5; d[0][1] = -1.0; d[0][2] = 0; d[0][3] = 1.0; d[0][4] = 1.5;
	d[1][0] = -1.0; d[1][1] = -0.5; d[1][2] = 0; d[1][3] = 0.5; d[1][4] = 1.0;
	d[2][0] = -0.6; d[2][1] = -0.3; d[2][2] = 0; d[2][3] = 0.3; d[2][4] = 0.6;

	// save smax to image
	FloatImageType::Pointer smaxImage = FloatImageType::New();
	smaxImage->SetRegions(REGION);
	smaxImage->SetSpacing(input->GetSpacing());
	smaxImage->Allocate();
	smaxImage->FillBuffer(0);
	FloatIteratorType smaxIt(smaxImage,REGION);

	int count=0;

	for (smaxIt.GoToBegin(), gradientIt.GoToBegin(), colonIt.GoToBegin(), vmapIt.GoToBegin(); !smaxIt.IsAtEnd(); 
		++smaxIt, ++gradientIt, ++colonIt, ++vmapIt)
	{

		if (colonIt.Get() == 255 && vmapIt.Get() == Unclassified)
		{
			if (count++ % 1000 == 0)
				std::cout << count << std::endl;

			ImageType::IndexType idx = gradientIt.GetIndex();

			VoxelType v;

			float dist = 0;

			for (int i=0; i<3; i++) // iterate through each set of distances
			{
				float in[5];
				float gr[5];
				
				for (int j=0; j<5; j++) // iterate through each distance
				{
					ContinuousIndexType odx = idx;

					for (int k=0; k<3; k++)
						odx[k] += gradientIt.Get()[k]*d[i][j];
					
					if ( !checkBounds(REGION,odx) )
						odx = idx;
					
					in[j] = inputInterp->EvaluateAtContinuousIndex(odx);
					gr[j] = gradInterp->EvaluateAtContinuousIndex(odx);
				}

				float distanceTS=0;
				float distanceSA=0;
				float distanceTA=0;

				float localDist=0;
				VoxelType localV=Unclassified;

				// Compute smax and validate by setting bounds of parabola
				
				float smaxSA = Stool_Air_ComputeSmax(in,gr,5);
				//float smaxTS = ComputeSmax(in,gr,5);
				float smaxTS=smaxSA; // enforce same smax

				smaxIt.Set(smaxSA);

				bool useTS = true;
				bool useSA = true;

				if (smaxTS < tissueStoolThreshold || smaxTS > 1500)
				{
					useTS = false;
				} else {
					distanceTS=AverageTissueStoolDist( smaxTS, in, gr );
				}

				if (smaxSA < tissueStoolThreshold || smaxSA > 1500)
				{
					useSA = false;
				} else {
					distanceSA=AverageStoolAirDist( smaxSA, in, gr );
				}

				distanceTA=AverageTissueAirDist(in,gr);				
				
				if (useSA && distanceSA<=distanceTS && distanceSA<=distanceTA) 
				{
					localDist = distanceSA;
					localV = StoolAir;

				} else if (useTS && distanceTS<=distanceTA && distanceTS<=distanceSA) {
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
			}

			vmapIt.Set(v);

		}		
	}
	
	Write(smaxImage,"smax.nii");
	Write(vmap,"qr.nii");

	return ComputePartials(input,vmap,colon,smaxImage);
}

/*
void LocalBoundary(VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon)
{
	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorType;
	
	VoxelImageType::SizeType radius;
	radius.Fill(0);
	radius[0]=2;
	radius[1]=2;

	NeighborhoodIteratorType nit(radius,vmap,REGION);
	ByteIteratorType colonIt(colon,REGION);

	int n = nit.Size();

	std::cout << n << std::endl;

	unsigned int voxelCount[3] = {0,0,0}; // air/tissue/stool

	for (nit.GoToBegin(), colonIt.GoToBegin(); !nit.IsAtEnd(); ++nit, ++colonIt)
	{
		if (colonIt.Get() == 255)
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

				//if (voxelCount[0] > voxelCount[2] && voxelCount[1] > voxelCount[2])
				//{
				//	nit.SetCenterPixel(TissueAir);
				//} else if (voxelCount[1] > voxelCount[0] && voxelCount[2] > voxelCount[0]) {
				//	nit.SetCenterPixel(TissueStool);
				//} else if (voxelCount[0] > voxelCount[1] && voxelCount[2] > voxelCount[1]) {
				//	nit.SetCenterPixel(StoolAir);
				//}
			}
		}
		
	}

	Write(vmap,"local_vmap.nii");
}	
*/