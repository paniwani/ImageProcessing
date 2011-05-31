#include "itkImage.h"
#include <iostream>
#include <utils.h>
#include <itkTextureImageToImageFilter.h>

int main()
{
	ImageType2D::Pointer input = ReadDicom <ImageType2D> ("C:/ImageData/mr10-uncleansed/mr10_092_13p.i0344/dcm",85);

	WriteITK <ImageType2D> (input,"input.nii");

	/*typedef itk::TextureImageToImageFilter<ImageType2D,ImageType2D> TextureFilterType;
	TextureFilterType::Pointer tf = TextureFilterType::New();
	tf->SetInput(input);
	tf->Update();

	for (unsigned int i=0; i<8; i++)
	{
		std::stringstream ss;
		ss << "texture" << i << ".nii";
		WriteITK <ImageType2D> (tf->GetOutput(i),ss.str());
	}*/

	system("pause");
	return 0;
}