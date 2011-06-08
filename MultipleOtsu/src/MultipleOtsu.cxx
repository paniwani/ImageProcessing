#include <itkImage.h> 						
#include <iostream> 		
#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <utils.h>
 												
int main(int argc, char * argv[])				
{ 						
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",130,200);
	WriteITK <ImageType> (input,"input.nii");
	
	typedef itk::OtsuMultipleThresholdsImageFilter<ImageType,ImageType> OtsuMultipleType;
	OtsuMultipleType::Pointer otsu = OtsuMultipleType::New();
	otsu->SetInput(input);
	otsu->SetNumberOfThresholds(2);
	otsu->Update();
	const OtsuMultipleType::OtsuCalculatorType::OutputType thresholds = otsu->GetThresholds();

	system("pause"); 							
	return 0; 									
} 												
