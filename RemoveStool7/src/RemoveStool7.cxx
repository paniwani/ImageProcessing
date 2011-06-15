#include <itkImage.h> 						
#include <iostream> 	
#include <utils.h>
#include <itkColonSegmentationFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkFastBilateralImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>

ByteImageType2D::Pointer FindTissue(ImageType2D::Pointer input, FloatImageType2D::Pointer gm)
{
	// Get region
	ImageType2D::RegionType region = input->GetLargestPossibleRegion();
	
	WriteITK <FloatImageType2D> (gm,"gm.nii");

	IteratorType2D inputIt(input,region);
	FloatIteratorType2D gmIt(gm,region);

	ByteImageType2D::Pointer tissue = ByteImageType2D::New();
	tissue->SetRegions(region);
	tissue->CopyInformation(input);
	tissue->Allocate();
	tissue->FillBuffer(0);
	ByteIteratorType2D tissueIt(tissue,region);

	for (inputIt.GoToBegin(), tissueIt.GoToBegin(), gmIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++tissueIt, ++gmIt)
	{
		if (inputIt.Get() > -250 && inputIt.Get() < 180 && gmIt.Get() < 300)
			tissueIt.Set( 255 );
	}

	//WriteITK <ByteImageType2D> (tissue,"tissue.nii");

	return tissue;

	//// Remove small components from tissue
	//typedef itk::BinaryShapeOpeningImageFilter<ByteImageType2D> OpeningType;
	//OpeningType::Pointer opener = OpeningType::New();
	//opener->SetInput(tissue);
	//opener->SetForegroundValue(255);
	//opener->SetBackgroundValue(0);
	//opener->SetLambda(300);
	//opener->SetAttribute("Size");
	////opener->SetFullyConnected(true);
	//opener->Update();
	//tissue = opener->GetOutput();

	//WriteITK <ByteImageType2D> (tissue,"tissueOpen.nii");
}

 												
int main(int argc, char * argv[])				
{ 		
	// Load input
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("C:/ImageData/mr10_092_13p.i0344/dcm/mr10_092_13p_i0140.dcm");
	WriteITK <ImageType2D> (input,"input.nii");

	ImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType2D,ByteImageType2D> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType2D::Pointer colon = colonSegmenter->GetOutput();
	WriteITK <ByteImageType2D> (colon,"colon.nii");

	// Crop region
	ImageType2D::SizeType size = region.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0;
	short paddingXY = 5;
	
	IteratorType2D inputIt(input,region);
	ByteIteratorType2D colonIt(colon,region);

	for (inputIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt)
	{
		if ( colonIt.Get() != 0 )
		{
			ImageType2D::IndexType idx = inputIt.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];
		}
	}

	ImageType2D::IndexType edx;
	
	edx[0] = (minX-paddingXY) > 0 ? minX-paddingXY : 0;
	edx[1] = (minY-paddingXY) > 0 ? minY-paddingXY : 0;

	ImageType2D::SizeType esize;
	esize[0] = maxX-minX+2*paddingXY+1 < size[0] ? maxX-minX+2*paddingXY+1 : size[0];
	esize[1] = maxY-minY+2*paddingXY+1 < size[1] ? maxY-minY+2*paddingXY+1 : size[1];

	ImageType2D::RegionType extractRegion;
	extractRegion.SetIndex( edx );
	extractRegion.SetSize( esize );

	typedef itk::RegionOfInterestImageFilter<ImageType2D,ImageType2D> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( input );
	cropper->SetRegionOfInterest( extractRegion );
	cropper->Update();
	input = cropper->GetOutput();

	typedef itk::RegionOfInterestImageFilter<ByteImageType2D,ByteImageType2D> RegionOfInterestByteImageFilterType;
	RegionOfInterestByteImageFilterType::Pointer cropperByte = RegionOfInterestByteImageFilterType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();

	region = input->GetLargestPossibleRegion();

	WriteITK <ImageType2D> (input,"inputCropped.nii");
	WriteITK <ByteImageType2D> (colon,"colonCropped.nii");

	// Mask input with colon
	typedef itk::MaskImageFilter<ImageType2D,ByteImageType2D> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(input);
	masker->SetInput2(colon);
	masker->SetOutsideValue(-1024);
	masker->Update();
	input = masker->GetOutput();

	WriteITK <ImageType2D> (input,"inputMasked.nii");

	// Get gradient magnitude
	typedef itk::GradientMagnitudeImageFilter<ImageType2D,FloatImageType2D> GradientMagnitudeImageFilterType;
	GradientMagnitudeImageFilterType::Pointer gradientFilter = GradientMagnitudeImageFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->Update();
	FloatImageType2D::Pointer gm = gradientFilter->GetOutput();
	WriteITK <FloatImageType2D> (gm,"gm.nii");

	ByteImageType2D::Pointer tissue = FindTissue(input,gm);
	WriteITK <ByteImageType2D> (tissue,"tissue.nii");

	// Show gm < 600 only
	FloatIteratorType2D gmIt(gm,region);
	for (gmIt.GoToBegin(); !gmIt.IsAtEnd(); ++gmIt)
	{
		if (gmIt.Get() > 650)
			gmIt.Set(0);
	}

	WriteITK <FloatImageType2D> (gm,"gmThresholded.nii");

	// Smooth with bilateral filter
	typedef itk::FastBilateralImageFilter<ImageType2D,ImageType2D> FastBilateralImageFilterType;
	FastBilateralImageFilterType::Pointer bilateralFilter = FastBilateralImageFilterType::New();
	bilateralFilter->SetInput(input);
	bilateralFilter->SetRangeSigma(100);
	bilateralFilter->SetDomainSigma(5);
	bilateralFilter->Update();
	input = bilateralFilter->GetOutput();
	WriteITK <ImageType2D> (input,"inputBilateral.nii");

	gradientFilter = GradientMagnitudeImageFilterType::New();
	gradientFilter->SetInput(input);
	gradientFilter->Update();
	FloatImageType2D::Pointer gmSmooth = gradientFilter->GetOutput();
	WriteITK <FloatImageType2D> (gmSmooth,"gmSmooth.nii");

	ByteImageType2D::Pointer tissueSmooth = FindTissue(input,gmSmooth);
	WriteITK <ByteImageType2D> (tissueSmooth,"tissueSmooth.nii");

	system("pause"); 							
	return 0; 									
} 			


