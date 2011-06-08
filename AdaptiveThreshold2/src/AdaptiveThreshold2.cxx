#include <itkImage.h> 						
#include <iostream> 		
#include <utils.h>
#include <itkAdaptiveOtsuThresholdImageFilter.h>
 												
int main(int argc, char * argv[])				
{ 						
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85,100);

	WriteITK <ImageType> (input,"input.nii");

	typedef itk::AdaptiveOtsuThresholdImageFilter<ImageType,ImageType> AOFilterType;
	AOFilterType::Pointer aoFilter = AOFilterType::New();
	aoFilter->SetInput(input);
	
	try {
		aoFilter->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error: " << err << std::endl;
		system("pause");
	} 

	WriteITK <ImageType> (aoFilter->GetOutput(), "out.nii");

	system("pause"); 							
	return 0; 									
} 												
