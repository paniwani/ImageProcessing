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

ByteImageType::Pointer BinaryKeeper(ByteImageType::Pointer &in, std::string attribute, unsigned int numOfObjects, bool reverse=false)
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

ByteImageType::Pointer BinaryOr(ByteImageType::Pointer &im1, ByteImageType::Pointer &im2)
{	
	typedef itk::OrImageFilter<ByteImageType> OrType;
	OrType::Pointer or = OrType::New();
	or->SetInput1(im1);
	or->SetInput2(im2);
	or->Update();
	return or->GetOutput();
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

ImageType::Pointer Crop(ImageType::Pointer &input, ImageType::RegionType requestedRegion)
{
	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> CropType;
	CropType::Pointer cropper = CropType::New();
	cropper->SetInput(input);
	cropper->SetRegionOfInterest(requestedRegion);
	cropper->Update();
	return cropper->GetOutput();
}




