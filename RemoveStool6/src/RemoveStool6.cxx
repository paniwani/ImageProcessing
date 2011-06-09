#include <itkImage.h> 						
#include <iostream> 		
#include <utils.h>
#include <itkBilateralImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkFastBilateralImageFilter.h>
//#include <itkColonSegmentationFilter.h>

ByteImageType::Pointer SegmentColon(ImageType::Pointer &input);
 												
int main(int argc, char * argv[])				
{ 	
	// get input
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",0,100);
	WriteITK <ImageType> (input,"input.nii");

	// 

	ByteImageType::Pointer colon = SegmentColon(input);
	
	// smooth
	/*double rangeSigma = 50.0;
	double domainSigma = 10.0;

	typedef itk::FastBilateralImageFilter<ImageType,ImageType> BilateralImageFilterType;
	BilateralImageFilterType::Pointer smoother = BilateralImageFilterType::New();
	smoother->SetInput(input);	
	smoother->SetRangeSigma(rangeSigma);
	smoother->SetDomainSigma(domainSigma);
	smoother->Update();
	ImageType::Pointer inputSmooth = smoother->GetOutput();
	WriteITK <ImageType> (inputSmooth,"inputSmooth.nii");

	WriteITK <ImageType> ( Subtract(inputSmooth,input) , "change.nii" );*/

	//system("pause"); 							
	return 0; 									
} 									

ByteImageType::Pointer SegmentColon(ImageType::Pointer &input)
{
	// get air
	PixelType min,max;
	GetMinMax(input,min,max);

	ByteImageType::Pointer air = BinaryThreshold(input,-600);
	WriteITK <ByteImageType> (air,"air.nii");
	
	UIntImageType::Pointer airMapU = CC(air);
	RelabelSize(airMapU,100);
	Threshold(airMapU,0,10);
	
	typedef itk::CastImageFilter<UIntImageType,ByteImageType> CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput(airMapU);
	caster->Update();
	ByteImageType::Pointer airMap = caster->GetOutput();

	airMapU.~SmartPointer();

	WriteITK <ByteImageType> (airMap,"airMap.nii");

	// remove bkg air
	// count number of times touches XY border
	RelabelBorder(airMap);
	WriteITK <ByteImageType> (airMap,"airMapBorder.nii");

	ByteIteratorType ait(airMap,airMap->GetLargestPossibleRegion());
	for (ait.GoToBegin(); !ait.IsAtEnd(); ++ait)
	{
		if (ait.Get() == 1)
			ait.Set(0);
	}

	WriteITK <ByteImageType> (airMap,"airMapNoBorder.nii");

	// remove lungs
	// find 2 components with lowest y centroid
	RelabelCentroidY(airMap);
	WriteITK <ByteImageType> (airMap,"airMapCY.nii");

	unsigned char mincy,maxcy;

	GetMinMax(airMap,mincy,maxcy);

	ait = ByteIteratorType(airMap,airMap->GetLargestPossibleRegion());

	for (ait.GoToBegin(); !ait.IsAtEnd(); ++ait)
	{
		if (ait.Get() >= maxcy-1)
		{
			ait.Set(0);
		} else if (ait.Get() > 0 && ait.Get() < maxcy-1) {
			ait.Set(255);
		}
	}

	air = airMap;

	WriteITK <ByteImageType> (airMap,"airMapNoLungs.nii");	

	// find tagged
	ByteImageType::Pointer tag = BinaryThreshold(input,250,true);
	WriteITK <ByteImageType> (tag,"tag.nii");

	UIntImageType::Pointer tagMapU = CC(tag);
	RelabelSize(tagMapU,5);

	WriteITK <UIntImageType> (tagMapU,"tag2.nii");

	return air;
}
