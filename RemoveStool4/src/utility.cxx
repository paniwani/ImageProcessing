ByteImageType::Pointer AllocateByteImage(ImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(VoxelImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(ByteImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ImageType::Pointer AllocateImage(ImageType::Pointer &in)
{
	ImageType::Pointer out = ImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer AllocateByteImage(FloatImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

FloatImageType::Pointer AllocateFloatImage(ImageType::Pointer &in)
{
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

FloatImageType::Pointer AllocateFloatImage(FloatImageType::Pointer &in)
{
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer BinaryThreshold(ImageType::Pointer &in, PixelType t1=itk::NumericTraits<PixelType>::NonpositiveMin(), PixelType t2=itk::NumericTraits<PixelType>::max())
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	IteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() > t1 && it.Get() < t2)
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(FloatImageType::Pointer &in, float t1=itk::NumericTraits<float>::NonpositiveMin(), float t2=itk::NumericTraits<float>::max())
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	FloatIteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() > t1 && it.Get() < t2)
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer BinaryThreshold(VoxelImageType::Pointer &vmap, VoxelType v)
{
	ByteImageType::Pointer out = AllocateByteImage(vmap);
	
	VoxelIteratorType it(vmap,vmap->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() == v )
		{
			oit.Set(255);
		}
	}
	
	return out;
}

ByteImageType::Pointer Mask(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<ByteImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

FloatImageType::Pointer Mask(FloatImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<FloatImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

ImageType::Pointer Mask(ImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::MaskImageFilter<ImageType,ByteImageType> MaskerType;
	MaskerType::Pointer masker = MaskerType::New();
	masker->SetInput1(im1);
	masker->SetInput2(im2);
	masker->Update();
	return masker->GetOutput();
}

ByteImageType::Pointer BinaryShapeKeeper(ByteImageType::Pointer &in, std::string attribute, unsigned int numOfObjects, bool reverse=false)
{
	typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType> KeeperType;
	KeeperType::Pointer keeper = KeeperType::New();
	keeper->SetInput(in);
	keeper->SetAttribute(attribute);
	keeper->SetNumberOfObjects(numOfObjects);
	keeper->SetReverseOrdering(reverse);
	keeper->SetForegroundValue(255);
	keeper->SetBackgroundValue(0);
	keeper->Update();
	return keeper->GetOutput();
}

ByteImageType::Pointer BinaryShapeOpen(ByteImageType::Pointer &in, std::string attribute, double lambda, bool reverse=false)
{
	typedef itk::BinaryShapeOpeningImageFilter<ByteImageType> OpenerType;
	OpenerType::Pointer opener = OpenerType::New();
	opener->SetInput(in);
	opener->SetAttribute(attribute);
	opener->SetLambda(lambda);
	opener->SetReverseOrdering(reverse);
	opener->SetForegroundValue(255);
	opener->SetBackgroundValue(0);
	opener->Update();
	return opener->GetOutput();
}

ByteImageType::Pointer BinaryOr(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::OrImageFilter<ByteImageType> OrType;
	OrType::Pointer or = OrType::New();
	or->SetInput1(im1);
	or->SetInput2(im2);
	or->Update();
	return or->GetOutput();
}

ByteImageType::Pointer BinaryInvert(ByteImageType::Pointer &im)
{	
	ByteImageType::Pointer out = AllocateByteImage(im);
	ByteIteratorType ot(out,out->GetLargestPossibleRegion());
	ByteIteratorType it(im,im->GetLargestPossibleRegion());

	for (it.GoToBegin(), ot.GoToBegin(); !it.IsAtEnd(); ++it, ++ot)
	{
		if (it.Get() == 255)
		{
			ot.Set(0);
		} else {
			ot.Set(255);
		}
	}

	return out;
}

ByteImageType::Pointer BinarySubtract(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{
	typedef itk::SubtractImageFilter<ByteImageType> SubtractType;
	SubtractType::Pointer subtracter = SubtractType::New();
	subtracter->SetInput1(im1);
	subtracter->SetInput2(im2);
	subtracter->Update();
	return subtracter->GetOutput();
}

template <typename T1, typename T2>
typename T2::Pointer Cast(typename T1::Pointer &in)
{
	typedef itk::CastImageFilter<T1,T2> CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput(in);
	caster->Update();
	return caster->GetOutput();
}

void Rescale(FloatImageType::Pointer &in, float min, float max)
{
	typedef itk::RescaleIntensityImageFilter<FloatImageType> RescaleType;
	RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput(in);
	rescaler->SetOutputMaximum(max);
	rescaler->SetOutputMinimum(min);
	rescaler->Update();
	in = rescaler->GetOutput();
}

ImageType::Pointer Crop(ImageType::Pointer input, ImageType::RegionType requestedRegion)
{
	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> CropType;
	CropType::Pointer cropper = CropType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(requestedRegion);
	cropper->Update();
	return cropper->GetOutput();
}

ByteImageType::Pointer Crop(ByteImageType::Pointer &input, ByteImageType::RegionType requestedRegion)
{
	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> CropType;
	CropType::Pointer cropper = CropType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(requestedRegion);
	cropper->Update();
	return cropper->GetOutput();
}

ImageType::Pointer Median(ImageType::Pointer &im, unsigned int rad)
{
	typedef itk::MedianImageFilter<ImageType,ImageType> MedianType;
	MedianType::Pointer medianFilter = MedianType::New();
	medianFilter->SetInput(im);

	ImageType::SizeType radius;
	radius.Fill(0);
	radius[0] = rad;
	radius[1] = rad;

	medianFilter->SetRadius(radius);
	medianFilter->Update();
	return medianFilter->GetOutput();
}

ByteImageType::Pointer Median(ByteImageType::Pointer &im, unsigned int rad)
{
	typedef itk::BinaryMedianImageFilter<ByteImageType,ByteImageType> MedianType;
	MedianType::Pointer medianFilter = MedianType::New();
	medianFilter->SetInput(im);
	medianFilter->SetBackgroundValue( 0 );
	medianFilter->SetForegroundValue( 255 );

	ByteImageType::SizeType radius;
	radius.Fill(0);
	radius[0] = rad;
	radius[1] = rad;

	medianFilter->SetRadius(radius);
	medianFilter->Update();
	return medianFilter->GetOutput();
}

ByteImageType::Pointer BinaryDilate(ByteImageType::Pointer &img, unsigned int radius)
{
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;

	se.SetRadius( rad );
	se.CreateStructuringElement();

	/*
	FlatStructuringElementType::RadiusType elementRadius;
	elementRadius.Fill(0);
	elementRadius[0] = radius;
	elementRadius[1] = radius;

	FlatStructuringElementType se = FlatStructuringElementType::Ball(elementRadius);
	*/

	typedef itk::BinaryDilateImageFilter<ByteImageType, ByteImageType, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( img );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();

	return dilater->GetOutput();
}

ByteImageType::Pointer BinaryErode(ByteImageType::Pointer &img, unsigned int radius)
{
	
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;

	se.SetRadius( rad );
	se.CreateStructuringElement();
	
	/*

	FlatStructuringElementType::RadiusType elementRadius;
	elementRadius.Fill(0);
	elementRadius[0] = radius;
	elementRadius[1] = radius;

	FlatStructuringElementType se = FlatStructuringElementType::Ball(elementRadius);

	*/

	typedef itk::BinaryErodeImageFilter<ByteImageType, ByteImageType, StructuringElementType> BinaryErodeImageFilterType;
	BinaryErodeImageFilterType::Pointer eroder = BinaryErodeImageFilterType::New();
	eroder->SetInput( img );
	eroder->SetKernel( se );
	eroder->SetForegroundValue(255);
	eroder->SetBackgroundValue(0);
	eroder->Update();

	return eroder->GetOutput();
}

ByteImageType::Pointer BinaryOpen(ByteImageType::Pointer &img, unsigned int radius)
{
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;

	se.SetRadius( rad );
	se.CreateStructuringElement();	

	typedef itk::BinaryMorphologicalOpeningImageFilter<ByteImageType, ByteImageType, StructuringElementType> BinaryMorphologicalOpeningImageFilterType;
	BinaryMorphologicalOpeningImageFilterType::Pointer opener = BinaryMorphologicalOpeningImageFilterType::New();
	opener->SetInput( img );
	opener->SetKernel( se );
	opener->SetForegroundValue(255);
	opener->Update();
	return opener->GetOutput();
}

/*
void BinaryFillHoles2D(ByteImageType::Pointer &im)
{
	// Invert image
	ByteIteratorType it(im,im->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() == 255)
		{
			it.Set(0);
		} else {
			it.Set(255);
		}
	}

	// Find largest background component in each 2D slice
	typedef itk::BinaryShapeKeepNObjectsImageFilter< ByteImageType2D > KeeperType2D;
	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, KeeperType2D > SliceKeeperType;

	KeeperType2D::Pointer keeper = KeeperType2D::New();
	keeper->SetForegroundValue(255);
	keeper->SetBackgroundValue(0);
	keeper->SetAttribute("Size");
	keeper->SetNumberOfObjects(1);
	
	SliceKeeperType::Pointer slicer = SliceKeeperType::New();
	slicer->SetInput(im);
	slicer->SetFilter(keeper);
	slicer->Update();

	ByteImageType::Pointer bkg = slicer->GetOutput();

	// Invert image back and fill holes
	ByteIteratorType bit(bkg,bkg->GetLargestPossibleRegion());

	for (it.GoToBegin(),bit.GoToBegin(); !it.IsAtEnd(); ++it,++bit)
	{
		if (it.Get() == 255)
		{
			it.Set(0);
		} else {
			it.Set(255);
		}

		if (bit.Get() == 0)
		{
			it.Set(255);
		}
	}
}
*/

void BinaryFillHoles2D(ByteImageType::Pointer &im)
{
	typedef itk::BinaryShapeOpeningImageFilter< ByteImageType2D > BinaryShapeOpeningImageFilter2D;
	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, BinaryShapeOpeningImageFilter2D > SliceBySliceImageFilterBackgroundType;

	BinaryShapeOpeningImageFilter2D::Pointer bkgFilter2D = BinaryShapeOpeningImageFilter2D::New();
	bkgFilter2D->SetAttribute("SizeOnBorder");
	bkgFilter2D->SetBackgroundValue(255);
	bkgFilter2D->SetForegroundValue(0);
	bkgFilter2D->SetLambda(1);
	bkgFilter2D->SetReverseOrdering(false);

	SliceBySliceImageFilterBackgroundType::Pointer bkgRemover = SliceBySliceImageFilterBackgroundType::New();
	bkgRemover->SetInput( im );
	bkgRemover->SetFilter( bkgFilter2D );
	bkgRemover->Update();
	im = bkgRemover->GetOutput();
}

void BinaryFillHoles(ByteImageType::Pointer &im)
{
	// Invert image
	ByteIteratorType it(im,im->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() == 255)
		{
			it.Set(0);
		} else {
			it.Set(255);
		}
	}

	// Find largest background component
	ByteImageType::Pointer bkg = BinaryShapeKeeper(im,"Size",1);

	// Invert image back and fill holes
	ByteIteratorType bit(bkg,bkg->GetLargestPossibleRegion());

	for (it.GoToBegin(),bit.GoToBegin(); !it.IsAtEnd(); ++it,++bit)
	{
		if (it.Get() == 255)
		{
			it.Set(0);
		} else {
			it.Set(255);
		}

		if (bit.Get() == 0)
		{
			it.Set(255);
		}
	}
}


