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

	system("pause"); 							
	return 0; 									
} 									

ByteImageType::Pointer SegmentColon(ImageType::Pointer &input)
{
	// get air
	PixelType min,max;
	GetMinMax(input,min,max);

	ByteImageType::Pointer air = BinaryThreshold(input,min,-600);
	
	// remove lungs
	UIntImageType::Pointer airMap = CC(air);
	RelabelSize(airMap,1000);
	Threshold(airMap,0,10);

	WriteITK <UIntImageType> (airMap,"airMap.nii");

	return air;
}
