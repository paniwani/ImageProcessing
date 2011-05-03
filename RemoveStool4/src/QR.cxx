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

float ComputeSmax(PixelType intensity[], float gradient_magnitude[], int size) 
{
	double Smax =0;
	double subVariable=0;
	double square =0;
	
	for (int i=0;i<size;i++) 
	{
		square = intensity[i]*intensity[i];
		Smax+=square*gradient_magnitude[i];
		subVariable+=square*intensity[i];
	}
	
	Smax-=4*subVariable;
	subVariable=0;
	
	for (int i=0;i<size;i++) 
	{
		square = intensity[i]*intensity[i];
		square *= square;
		subVariable+=square;
	}

	return -4*subVariable/Smax;
}

float Stool_Air_ComputeSmax(PixelType intensity[], float gradient_magnitude[], int size) 
{
	PixelType intensity2[5];
	for (int i=0;i<size;i++) 
	{
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
float AverageTissueAirDist(PixelType intensity[], float gradient_magnitude[]) 
{
	double coefficients[3]={-4.0/1000,-4.0,0};

	vnl_real_polynomial poly(coefficients,3);
    double average=0;
    for(int i=0;i<5;i++) 
	{
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }
    return average/5;
}
float AverageTissueStoolDist(PixelType Smax, PixelType intensity[], float gradient_magnitude[]) 
{
    double coefficients[3] ={-4.0/Smax,4.0,0};
    vnl_real_polynomial poly(coefficients,3);
    double average=0;
    
	for(int i=0;i<5;i++) 
	{
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }

    return average/5;
}

float AverageStoolAirDist(PixelType Smax, PixelType intensity[], float gradient_magnitude[]) 
{
    double coefficients[3] ={-4/(Smax+1000),-4*(1000-Smax)/(Smax+1000),4000*Smax/(Smax+1000)};
    vnl_real_polynomial poly(coefficients,3);
    double average=0;
    
	for(int i=0;i<5;i++) 
	{
        average+=PolyMinDist(poly,intensity[i],gradient_magnitude[i]);
    }

    return average/5;
}



ImageType::Pointer ComputeNeighborhoodSmax(ImageType::Pointer &input, VoxelImageType::Pointer &v, ByteImageType::Pointer &mask)
{
	ImageType::Pointer smax = ImageType::New();
	smax->SetSpacing( input->GetSpacing() );
	smax->SetDirection( input->GetDirection() );
	smax->SetRegions( region );
	smax->CopyInformation( input );
	smax->Allocate();
	smax->FillBuffer(0);



	IteratorType smax_iter(smax,region);
	IteratorType input_iter(input,region);
	ByteIteratorType mask_iter(mask,region);

	typedef itk::NeighborhoodIterator<VoxelImageType> NeighborhoodIteratorVoxelType;
	NeighborhoodIteratorVoxelType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 2;
	radius[1] = 2;

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
				
				if ( region.IsInside( idx ) )
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