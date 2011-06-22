template <typename T1, typename T2>
typename T2::Pointer Rescale(typename T1::Pointer &im, typename T2::PixelType v1, typename T2::PixelType v2)
{
	typedef itk::RescaleIntensityImageFilter<T1,T2> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(im);
	rescaler->SetOutputMinimum(v1);
	rescaler->SetOutputMaximum(v2);
	rescaler->Update();
	return rescaler->GetOutput();
}

template <typename T1, typename T2>
typename T2::Pointer Cast(typename T1::Pointer &im)
{
	typedef itk::CastImageFilter<T1,T2> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(im);
	caster->Update();
	return caster->GetOutput();
}

template <typename T1, typename T2>
void Mask(typename T1::Pointer &im, typename T2::Pointer &m, typename T1::PixelType outsideValue=0)
{
	typedef itk::MaskImageFilter<T1,T2> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im);
	masker->SetInput2(m);
	masker->SetOutsideValue(outsideValue);
	masker->Update();
	im = masker->GetOutput();
}

template <typename T1, typename T2>
typename T1::Pointer AllocateImage(typename T2::Pointer &cim)
{
	T1::Pointer im = T1::New();
	im->SetRegions(cim->GetLargestPossibleRegion());
	im->CopyInformation(cim);
	im->Allocate();
	im->FillBuffer(0);
	return im;
}

template <typename T>
void BinaryDilate(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<T, T, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( im );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();
	im = dilater->GetOutput();
}

template <typename T>
typename T::Pointer BinaryErode(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter<T, T, StucturingElementType> BinaryErodeImageFilterType;
	BinaryErodeImageFilterType::Pointer eroder = BinaryErodeImageFilterType::New();
	eroder->SetInput( im );
	eroder->SetKernel( se );
	eroder->SetForegroundValue(255);
	eroder->SetBackgroundValue(0);
	eroder->Update();
	return eroder->GetOutput();
}

template <typename T>
typename T::Pointer BinaryClose(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryMorphologicalClosingImageFilter<T, T, StucturingElementType> BinaryMorphologicalClosingImageFilterType;
	BinaryMorphologicalClosingImageFilterType::Pointer closer = BinaryMorphologicalClosingImageFilterType::New();
	closer->SetInput( im );
	closer->SetKernel( se );
	closer->SetForegroundValue(255);
	closer->Update();
	return closer->GetOutput();
}

template <typename T>
typename T::Pointer BinaryOpen(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::BinaryBallStructuringElement<T::PixelType,Dimension> StucturingElementType;
	StucturingElementType se;

	T::SizeType rad;
	rad.Fill(0);
	rad[0] = radius; rad[1] = radius;

	if (Dimension == 3 && use3D)
		rad[2] = radius;

	se.SetRadius(rad);
	se.CreateStructuringElement();

	typedef itk::BinaryMorphologicalOpeningImageFilter<T, T, StucturingElementType> BinaryMorphologicalOpeningImageFilterType;
	BinaryMorphologicalOpeningImageFilterType::Pointer opener = BinaryMorphologicalOpeningImageFilterType::New();
	opener->SetInput( im );
	opener->SetKernel( se );
	opener->SetForegroundValue(255);
	opener->Update();
	return opener->GetOutput();
}

template <typename T>
typename T::RegionType BinaryCrop(typename T::Pointer &im, unsigned int pad=5, bool cropZ=false)
{
	T::RegionType region = im->GetLargestPossibleRegion();
	T::SizeType size = region.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0,minZ=0,maxZ=0;
	
	if (Dimension == 3 && cropZ)
	{
		minZ=size[2];
		maxZ=0;
	}

	typedef itk::ImageRegionIteratorWithIndex<T> IteratorType;
	IteratorType it(im,region);

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if ( it.Get() != 0 )
		{
			T::IndexType idx = it.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];

			if (Dimension == 3 && cropZ)
			{
				if (idx[2] < minZ)
					minZ = idx[2];
				if (idx[2] > maxZ)
					maxZ = idx[2];
			}
		}
	}

	T::IndexType edx;
	
	edx[0] = (minX-pad) > 0 ? minX-pad : 0;
	edx[1] = (minY-pad) > 0 ? minY-pad : 0;

	if (Dimension == 3)
	{
		if (cropZ)
		{
			edx[2] = (minZ-pad) > 0 ? minZ-pad : 0;
		} else {
			edx[2] = 0;
		}
	}

	T::SizeType esize;
	esize[0] = maxX-minX+2*pad+1 < size[0] ? maxX-minX+2*pad+1 : size[0];
	esize[1] = maxY-minY+2*pad+1 < size[1] ? maxY-minY+2*pad+1 : size[1];
	
	if (Dimension == 3)
	{
		if (cropZ)
		{
			esize[2] = maxZ-minZ+2*pad+1 < size[2] ? maxZ-minZ+2*pad+1 : size[2];
		} else {
			esize[2] = size[2];
		}
	}

	T::RegionType ror;
	ror.SetIndex( edx );
	ror.SetSize( esize );

	typedef itk::RegionOfInterestImageFilter<T,T> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( im );
	cropper->SetRegionOfInterest( ror );
	cropper->Update();
	im = cropper->GetOutput();

	return ror;
}

template <typename T>
typename T::Pointer CropByRegion(typename T::Pointer &im, typename T::RegionType region)
{
	typedef itk::RegionOfInterestImageFilter<T,T> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( im );
	cropper->SetRegionOfInterest( region );
	cropper->Update();
	return cropper->GetOutput();
}	

template <typename T>
typename T::Pointer Duplicate(typename T::Pointer &im)
{
	typedef itk::ImageDuplicator<T> DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(im);
	duplicator->Update();
	return duplicator->GetOutput();
} 

template <typename T>
typename T::Pointer Add(typename T::Pointer &im1, typename T::Pointer &im2)
{
	typedef itk::AddImageFilter<T> AdderType;
	AdderType::Pointer adder = AdderType::New();
	adder->SetInput1(im1);
	adder->SetInput2(im2);
	adder->Update();
	return adder->GetOutput();
}	

template <typename T>
void DivideByConstant(typename T::Pointer &im, typename T::PixelType c)
{
	typedef itk::DivideByConstantImageFilter<T,T::PixelType,T> DividerType;
	DividerType::Pointer divider = DividerType::New();
	divider->SetInput(im);
	divider->SetConstant(c);
	divider->Update();
	im = divider->GetOutput();
}

template <typename T>
typename T::Pointer Median(typename T::Pointer &im, unsigned int radius, bool use3D=false)
{
	typedef itk::MedianImageFilter<ImageType,ImageType> MedianType;
	MedianType::Pointer medianFilter = MedianType::New();
	medianFilter->SetInput(im);

	ImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = 1;
	rad[1] = 1;
	
	if (Dimension == 3 && use3D)
		rad[2] = radius;

	medianFilter->SetRadius(rad);
	medianFilter->Update();
	return medianFilter->GetOutput();
}