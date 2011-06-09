#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <itkTextureImageToImageFilter.h>
 												
int main(int argc, char * argv[])				
{ 		
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85);
	WriteITK <ImageType> (input,"input.nii");

	typedef itk::RescaleIntensityImageFilter<ImageType,ByteImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(input);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->SetReleaseDataFlag(true);
	
	typedef itk::TextureImageToImageFilter<ByteImageType> TextureFilterType;
	typedef TextureFilterType::OutputImageType OutputImageType;

	TextureFilterType::Pointer textureFilter = TextureFilterType::New();
	textureFilter->SetInput(rescaler->GetOutput());
	textureFilter->SetNeighborhoodRadius( 1 );

	std::vector<unsigned int> offsetScales(1,1);

	textureFilter->SetOffsetScales( offsetScales );
	textureFilter->SetNumberOfHistogramBins( 128 );
	textureFilter->SetNormalizeHistogram( false );
	textureFilter->Update();

	for (int i=0; i<8; i++)
	{
		std::stringstream ss;
		ss << "texture" << i << ".nii";
		WriteITK <OutputImageType> (textureFilter->GetOutput(i),ss.str());
	}

	system("pause"); 							
	return 0; 									
} 												
