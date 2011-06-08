#include <itkImage.h> 						
#include <iostream> 	
#include <utils.h>
#include <itkImageRegionMultidimensionalSplitter.h>
 												
int main(int argc, char * argv[])				
{ 						
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/imagedata/mr10-uncleansed/mr10_092_13p.i0344/dcm",85,90);

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	typedef itk::ImageRegionMultidimensionalSplitter<3> SplitterType;
	SplitterType::Pointer splitter = SplitterType::New();
	
	unsigned int requestedSplits = 50;
	
	unsigned int divisions = splitter->GetNumberOfSplits(region,50);


	system("pause"); 							
	return 0; 									
} 												
