float PolyDist(float X, vnl_real_polynomial poly, float x, float y) 
{
    return sqrtf(vnl_math_sqr(X-x)+vnl_math_sqr(poly.evaluate(X)-y));
}

double cubeRoot(double x)
{
    if (x < 0)
        return -pow(-x, 0.333333333333333);
    else
        return pow(x, 0.333333333333333);
}

// use cardano's formula
std::vector<double> solveCubic(double a, double b, double c, double d)
{
	std::vector<double> out;

    // find the discriminant
    double f, g, h;
    f = (3 * c / a - pow(b, 2) / pow(a, 2)) / 3;
    g = (2 * pow(b, 3) / pow(a, 3) - 9 * b * c / pow(a, 2) + 27 * d / a) / 27;
    h = pow(g, 2) / 4 + pow(f, 3) / 27;
    // evaluate discriminant
    if (f == 0 && g == 0 && h == 0)
    {
        // 3 equal roots
        double x;
        // when f, g, and h all equal 0 the roots can be found by the following line
        x = -cubeRoot(d / a);
        // print solutions
        /*cout
            << "x = " << endl
            << " " << x << endl
            << " " << x << endl
            << " " << x << endl << endl;*/

		out.push_back(x);
    }
    else if (h <= 0)
    {
        // 3 real roots
        double q, i, j, k, l, m, n, p;
        // complicated maths making use of the method
        i = pow(pow(g, 2) / 4 - h, 0.5);
        j = cubeRoot(i);
        k = acos(-(g / (2 * i)));
        m = cos(k / 3);
        n = 1.73205080756888 * sin(k / 3);
        p = -(b / (3 * a));
        // print solutions
        /*cout
            << "x = " << endl
            << " " << 2 * j * m + p << endl
            << " " << -j * (m + n) + p << endl
            << " " << -j * (m - n) + p << endl << endl;*/

		out.push_back(2 * j * m + p);
		out.push_back(-j * (m + n) + p);
		out.push_back(-j * (m - n) + p);
    }
    else if (h > 0)
    {
        // 1 real root and 2 complex roots
        double r, s, t, u, p;
        // complicated maths making use of the method
        r = -(g / 2) + pow(h, 0.5);
        s = cubeRoot(r);
        t = -(g / 2) - pow(h, 0.5);
        u = cubeRoot(t);
        p = -(b / (3 * a));
        // print solutions
        /*cout
            << "x = " << endl
            << " " << (s + u) + p << endl
            << " " << -(s + u) / 2 + p << " +" << (s - u) * ROOTTHREE / 2 << "i" << endl
            << " " << -(s + u) / 2 + p << " " << -(s - u) * ROOTTHREE / 2 << "i" << endl << endl;*/
		out.push_back((s + u) + p);
    }

	return out;
}

float PolyMinDist(vnl_real_polynomial poly, float x, float y)
{
	std::vector<double> rRoots = solveCubic(vnl_math_sqr(poly[0])*2,
											poly[0]*3*poly[1],
											(2*poly[0]*poly[2]-2*poly[0]*y+vnl_math_sqr(poly[1])+1),
											 poly[1]*poly[2]-poly[1]*y-x
											);

	int size = rRoots.size();

	float distance = PolyDist(rRoots[0],poly,x,y);
    
	for(int i=1;i<size;i++) {
        float temp_dist=PolyDist(rRoots[i],poly,x,y);
        if (distance>temp_dist) {
            distance=temp_dist;
        }
    }

    return distance;
}


//float PolyMinDist(vnl_real_polynomial poly, float x, float y)
//{
//    double coefficients[4] ={
//		vnl_math_sqr(poly[0])*2,
//		poly[0]*3*poly[1],
//		(2*poly[0]*poly[2]-2*poly[0]*y+vnl_math_sqr(poly[1])+1),
//        poly[1]*poly[2]-poly[1]*y-x
//    };
//	
//    vnl_real_polynomial polyRoot(coefficients,4);
//    vnl_rpoly_roots rooter(polyRoot);
//    //rooter.compute();
//    vnl_vector<double> rRoots = rooter.realroots(0.0005);
//    int size = rRoots.size();
//
//	float distance = PolyDist(rRoots(0),poly,x,y);
//    for(int i=1;i<size;i++) {
//        float temp_dist=PolyDist(rRoots(i),poly,x,y);
//        if (distance>temp_dist) {
//            distance=temp_dist;
//        }
//    }
//	//std::cerr<<distance<<std::endl;
//    return distance;
//}

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

					// check divide by 0
					if ( S == 0 )
					{
						value[2] = 0;
					} else {
						value[2]=I/S;
					}

					if (value[2] <= 0) { value[2] = 0; }
					if (value[2] >= 1) { value[2] = 1; }

					value[1]=1-value[2];

					break; 

				case StoolAir:

					if ( (S+1000) == 0 )
					{
						value[2] = 0;
					} else {
						value[2]=(I+1000)/(S+1000);
					}

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
FloatImageType::Pointer QuadraticRegression(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, FloatImageType::Pointer &gradientMagnitude, PixelType tissueStoolThreshold)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// save computed smax
	FloatImageType::Pointer smaxImage = AllocateFloatImage(input);
	FloatIteratorType smaxIt(smaxImage,region);
	
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

	// distances (in mm?)
	float d[3][5];
	d[0][0] = -1.5; d[0][1] = -1.0; d[0][2] = 0; d[0][3] = 1.0; d[0][4] = 1.5;
	d[1][0] = -1.0; d[1][1] = -0.5; d[1][2] = 0; d[1][3] = 0.5; d[1][4] = 1.0;
	d[2][0] = -0.6; d[2][1] = -0.3; d[2][2] = 0; d[2][3] = 0.3; d[2][4] = 0.6;

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
						odx[k] += gradientIt.Get()[k] * ( d[i][j] / spacing[k] );
					
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
				
				//float smaxSA = (float) smaxIt.Get();
				float smaxTS=smaxSA; // enforce same smax

				smaxIt.Set(smaxSA);

				bool useTS = true;
				bool useSA = true;

				if (smaxTS < tissueStoolThreshold || smaxTS > 1500)
				{
					useTS = false;
					
					smaxIt.Set(0);

				} else {
					distanceTS=AverageTissueStoolDist( smaxTS, in, gr );
				}

				if (smaxSA < tissueStoolThreshold || smaxSA > 1500)
				{
					useSA = false;

					smaxIt.Set(0);

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

	/*for (smaxIt.GoToBegin(), smax2It.GoToBegin(); !smaxIt.IsAtEnd(); ++smaxIt, ++smax2It)
	{
		if (smaxIt.Get() < 0)
			smaxIt.Set(0);

		if (smaxIt.Get() > 1200)
			smaxIt.Set(1300);

		if (smax2It.Get() < 0)
			smax2It.Set(0);

		if (smax2It.Get() > 1200)
			smax2It.Set(1300);
	}*/
	
	Write(smaxImage,"smax.nii");
	//Write(smaxImage2,"smax2.nii");
	Write(vmap,"qr.nii");

	return smaxImage;
}

/*
* Find tissue that was erroneusly marked SA
* by the QR using standard deviation and
* assign it a tissue a partial based on intensity
*/

void FixATT(ImageType::Pointer &input, ArrayImageType::Pointer &partial, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, FloatImageType::Pointer &smax, PixelType tissueStoolThreshold)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// get stool air mask
	ByteImageType::Pointer sam = BinaryThreshold(vmap,StoolAir);
	//Write(sam,"sam.nii");

	// get sd of input in stool air mask
	FloatImageType::Pointer sdF = StandardDeviation(input,sam,1);

	// cast
	typedef itk::CastImageFilter<FloatImageType,ImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(sdF);
	caster->Update();
	ImageType::Pointer sd = caster->GetOutput();
	Write(sd,"sd.nii");

	// get nonzero min, max of sd
	IteratorType sdIt(sd,region);
	float min = itk::NumericTraits<PixelType>::max();
	float max = 0;

	for (sdIt.GoToBegin(); !sdIt.IsAtEnd(); ++sdIt)
	{
		float s = sdIt.Get();

		if (s > 0)
		{
			if (s < min)
				min = s;
			if (s > max)
				max = s;
		}
	}

	std::cout << "sd min: " << min << std::endl;
	std::cout << "sd max: " << max << std::endl;

	// otsu threshold sd
	typedef itk::OtsuThresholdImageCalculatorModified<ImageType> OtsuType;
	OtsuType::Pointer otsu = OtsuType::New();
	otsu->SetImage(sd);
	otsu->SetMinMax(true);
	otsu->SetHistogramMin(min);
	otsu->SetHistogramMax(max);
	//otsu->SetPrintHistogram("sdHistogram.csv");
	otsu->Compute();
	float ot = otsu->GetThreshold();

	std::cout << "sd otsu: " << ot << std::endl;

	ByteImageType::Pointer sdmask = AllocateByteImage(input);
	ByteIteratorType sdmaskIt(sdmask,region);

	for (sdIt.GoToBegin(), sdmaskIt.GoToBegin(); !sdIt.IsAtEnd(); ++sdIt, ++sdmaskIt)
	{
		if (sdIt.Get() > min && sdIt.Get() < ot)
		{
			sdmaskIt.Set(255);
		}
	}

	//Write(sdmask,"sdmask.nii");

	// get input in tissue intensity range mask
	ByteImageType::Pointer timask = BinaryThreshold(input,-300,tissueStoolThreshold);
	//Write(timask,"timask.nii"); 

	// mask sd with tissue intensity mask
	sdmask = Mask(sdmask,timask);
	//Write(sdmask,"sdmask2.nii");

	// ------------------------------------------------------------------------
	// for each sd component, count number of times it touches tissue
	// ------------------------------------------------------------------------

	// get connected components		
	typedef itk::ConnectedComponentImageFilter<ByteImageType, LabelImageType> ConnectedType;
	ConnectedType::Pointer connecter = ConnectedType::New();
	connecter->SetInput(sdmask);
	connecter->SetBackgroundValue(0);
	
	typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType> RelabelType;
	RelabelType::Pointer relabeler = RelabelType::New();
	relabeler->SetInput(connecter->GetOutput());
	relabeler->SetMinimumObjectSize(3);
	relabeler->Update();

	unsigned long originalNumOfObjects = relabeler->GetOriginalNumberOfObjects();
	unsigned long numOfObjects = relabeler->GetNumberOfObjects();

	std::cout << "Original # of objects: " << originalNumOfObjects << std::endl;
	std::cout << "Number of objects: " << numOfObjects << std::endl;

	LabelImageType::Pointer cc = relabeler->GetOutput();
	LabelIteratorType ccIt(cc,region);

	Write(cc,"ccRelabel.nii");

	std::vector<unsigned long> countVector;
	countVector.resize(numOfObjects+1);

	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorType;

	VoxelImageType::SizeType radius;
	radius.Fill(1);

	NeighborhoodIteratorType vIt(radius,vmap,region);

	for (vIt.GoToBegin(), ccIt.GoToBegin(); !vIt.IsAtEnd(); ++vIt, ++ccIt)
	{
		if ( ccIt.Get() > 0 )
		{
			bool tissueNeighbor = false;

			for (int i=0; i<vIt.Size(); i++)
			{
				if ( i != (vIt.Size()-1)/2 )
				{
					if ( vIt.GetPixel(i) == Tissue )
					{
						tissueNeighbor = true;
						break;
					}
				}
			}

			if (tissueNeighbor)
				countVector[ ccIt.Get() ]++;
		}
	}

	sdmaskIt = ByteIteratorType(sdmask,region);

	for (ccIt.GoToBegin(), sdmaskIt.GoToBegin(); !ccIt.IsAtEnd(); ++ccIt, ++sdmaskIt)
	{
		if ( countVector[ ccIt.Get() ] < 10 )
		{
			sdmaskIt.Set( 0 );	
		}
	}

	//Write(sdmask,"sdmaskPost.nii");

	// compute partial of tissue

	IteratorType inputIt(input,region);
	FloatIteratorType smaxIt(smax,region);
	ArrayIteratorType partialIt(partial,region);

	for (inputIt.GoToBegin(), smaxIt.GoToBegin(), partialIt.GoToBegin(), sdmaskIt.GoToBegin(); !inputIt.IsAtEnd();
		 ++inputIt, ++smaxIt, ++partialIt, ++sdmaskIt)
	{
		if (sdmaskIt.Get() != 0)
		{
			float I = (float) inputIt.Get();
			float S = (float) smaxIt.Get();

			ArrayType p;

			if ( S > 0 )
			{
				p[0] = 0;

				p[1] = 1 - (I/S);

				// bounds check
				p[1] = (p[1] > 1) ? 1 : p[1];
				p[1] = (p[1] < 0) ? 0 : p[1];

				p[2] = 1 - p[1];

			} else {
				p[0] = 0;
				p[1] = 1;
				p[2] = 0;
			}

			partialIt.Set( p );
		}
	}

	Write(partial,"partialATTFix.nii");
}










//void FixATT(ImageType::Pointer &input, ArrayImageType::Pointer &partial, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon)
//{
//	ImageType::RegionType region = input->GetLargestPossibleRegion();
//
//	// make stool air binary mask
//	ByteImageType::Pointer sa = ByteImageType::New();
//	sa->SetRegions(REGION);
//	sa->SetSpacing(input->GetSpacing());
//	sa->CopyInformation(input);
//	sa->Allocate();
//	sa->FillBuffer(0);
//	ByteIteratorType saIt(sa,REGION);
//	VoxelIteratorType vmapIt(vmap,REGION);
//
//	for (saIt.GoToBegin(), vmapIt.GoToBegin(); !saIt.IsAtEnd(); ++saIt, ++vmapIt)
//	{
//		if (vmapIt.Get() == StoolAir)
//		{
//			saIt.Set(255);
//		}
//	}
//
//	Write(sa,"sa.nii");
//
//	// find tissue automatically
//	// compute sd within stool air mask smoothed
//	/*float sigma[5] = {0.1, 0.25, 0.5, 1, 2};
//
//	for (int i=0; i<5; i++)
//	{*/
//		//// smooth input
//		//typedef itk::DiscreteGaussianImageFilter<ImageType,FloatImageType> GaussianType;
//		//GaussianType::Pointer smoother = GaussianType::New();
//		//smoother->SetInput(input);
//		//smoother->SetVariance( sigma[i]*sigma[i] );
//		//smoother->Update();
//		//FloatImageType::Pointer smoothedInput = smoother->GetOutput();
//
//		// get sd within SA mask
//		FloatImageType::Pointer sasd = StandardDeviation ( input, sa, 1 );
//		std::stringstream ss;
//		//ss << "stoolAirSD_rad_1_sigma_" << sigma[i] << ".nii";
//		Write(sasd,"sasd.nii");
//
//
//		// tissue can only be within strict intensity range
//		typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> 
//
//
//		// otsu
//
//		// get non-zero min and max of image
//		FloatIteratorType sasdIt(sasd,region);
//		float min = itk::NumericTraits<float>::max();
//		float max = 0;
//
//		for (sasdIt.GoToBegin(); !sasdIt.IsAtEnd(); ++sasdIt)
//		{
//			float val = sasdIt.Get();
//
//			if ( val > 0 )
//			{
//				if (val < min)
//					min = val;
//
//				if (val > max)
//					val = max;
//			}
//		}
//
//		std::cout << "min: " << min << std::endl;
//		std::cout << "max: " << max << std::endl;
//
//		// otsu threshold
//		typedef itk::OtsuThresholdImageCalculatorModified<FloatImageType> OtsuType;
//		OtsuType::Pointer otsuCalc = OtsuType::New();
//		otsuCalc->SetImage(sasd);
//		otsuCalc->SetMinMax(true);
//		otsuCalc->SetHistogramMin(min);
//		otsuCalc->SetHistogramMax(max);
//		otsuCalc->SetPrintHistogram("otsu.csv");
//		otsuCalc->Compute();
//		float ot = otsuCalc->GetThreshold();
//
//		std::cout << "Otsu: " << ot << std::endl;
//
//		// threshold
//		typedef itk::BinaryThresholdImageFilter<FloatImageType,ByteImageType> ThresholdType;
//		ThresholdType::Pointer thresholder = ThresholdType::New();
//		thresholder->SetInput(sasd);
//		thresholder->SetLowerThreshold( ot );
//		thresholder->SetInsideValue(0);
//		thresholder->SetOutsideValue(255);
//		thresholder->Update();
//
//		ss.str("");
//		ss << "sasdOtsu_sigma_" << sigma[i] << ".nii";
//		Write(thresholder->GetOutput(),ss.str());		
//
//		std::cout << std::endl;
//	}	
//	//FloatImageType::Pointer sasd = StandardDeviation(input,sa,1);
//	//Write(sasd,"stoolAirSD_rad_1.nii");
//
//	//// get standard deviation histogram
//	//typedef itk::Statistics::ScalarImageToHistogramGenerator< FloatImageType > HistogramGeneratorType;
//	//HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
//	//histogramGenerator->SetInput(sasd);
//	//histogramGenerator->SetNumberOfBins( 256 );
//	//histogramGenerator->SetMarginalScale( 10.0 );
//
//	//histogramGenerator->SetHistogramMin(  4.5 );
//	//histogramGenerator->SetHistogramMax( 1500.5 );
//	//histogramGenerator->Compute();
//
//	//typedef HistogramGeneratorType::HistogramType  HistogramType;
//
//	//const HistogramType * histogram = histogramGenerator->GetOutput();
//
//	//// write histogram to txt file
//	//const unsigned int histogramSize = histogram->Size();
//	//std::cout << "Histogram size " << histogramSize << std::endl;
//
//	//std::ofstream sdfile;
//	//sdfile.open("sdfile.csv");
//
//	//unsigned int bin;
//	//for( bin=1; bin < histogramSize; bin++ )
//	//{
//	//	sdfile << bin + 5 << "," << histogram->GetFrequency(bin,0) << "\n";
//	//}
//	//sdfile.close();
//
//	//// convert histogram to 1D image
//	//ImageType1D::Pointer histImage = ImageType1D::New();
//	//
//	//ImageType1D::IndexType index;
//	//index[0] = 0;
//	//
//	//ImageType1D::SizeType size;
//	//size[0] = 256;
//	//
//	//ImageType1D::RegionType region1D;
//	//region1D.SetIndex(index);
//	//region1D.SetSize(size);
//
//	//histImage->SetRegions(region1D);
//	//histImage->Allocate();
//	//histImage->FillBuffer(0);
//	//IteratorType1D histImageIt(histImage,region1D);
//
//	//histImageIt.GoToBegin();
//
//	//HistogramType::ConstIterator itr = histogram->Begin();
// //   HistogramType::ConstIterator end = histogram->End();
// //
//	//unsigned int binNumber = 0;
//	//while( itr != end && !histImageIt.IsAtEnd() )
//	//{
//	//	histImageIt.Set( itr.GetFrequency() );
//	//	
//	//	++itr;
//	//	++histImageIt;
//	//	++binNumber;
//	//}
//
//	//Write(histImage,"histImage.nii");
//
//	//// smooth 1D image
//	//typedef itk::DiscreteGaussianImageFilter<ImageType1D,ImageType1D> SmootherType;
//	//SmootherType::Pointer smoother = SmootherType::New();
//	//smoother->SetInput(histImage);
//	//smoother->SetVariance(1);
//
//	//// convert histogram image back into a histogram
//	//typedef itk::Statistics::ImageToHistogramFilter<ImageType1D> ImageToHistogramFilterType;
//	//ImageToHistogramFilterType::Pointer image2histFilter = ImageToHistogramFilterType::New();
//	//image2histFilter->SetInput(smoother->GetOutput());
//	//image2histFilter->Update();
//	//
//	//const HistogramType * smoothHistogram = image2histFilter->GetOutput();
//
//	//// write histogram to txt file
//	//const unsigned int histogramSize2 = smoothHistogram->Size();
//	//std::cout << "Histogram size " << histogramSize2 << std::endl;
//
//	//std::ofstream sdfile2;
//	//sdfile2.open("sdfile2.csv");
//
//	//unsigned int bin2;
//	//for( bin2=1; bin2 < histogramSize; bin2++ )
//	//{
//	//	sdfile2 << bin2 + 5 << "," << histogram->GetFrequency(bin2,0) << "\n";
//	//}
//	//sdfile2.close();
//}


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
