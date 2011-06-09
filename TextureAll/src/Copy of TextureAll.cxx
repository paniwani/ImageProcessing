#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <otbScalarImageToTexturesFilter.h>
#include <time.h>
#include <itkAddImageFilter.h>
#include <itkDivideByConstantImageFilter.h>
 											
FloatImageType2D::Pointer DivideByConstant(FloatImageType2D::Pointer &im1, float c)
{
	typedef itk::DivideByConstantImageFilter<FloatImageType2D,float,FloatImageType2D> DividerType;
	DividerType::Pointer divider = DividerType::New();
	divider->SetInput(im1);
	divider->SetConstant(c);
	divider->Update();
	return divider->GetOutput();
}

FloatImageType2D::Pointer Add(FloatImageType2D *im1, FloatImageType2D *im2)
{
	typedef itk::AddImageFilter<FloatImageType2D> AdderType;
	AdderType::Pointer adder = AdderType::New();
	adder->SetInput1(im1);
	adder->SetInput2(im2);
	adder->Update();
	return adder->GetOutput();
}

int main(int argc, char * argv[])				
{ 		
	// Start the clock
	clock_t init = clock();

	// Get image
	FloatImageType2D::Pointer input = ReadITK <FloatImageType2D> ("ADS40RoiSmall.png");
	WriteITK <FloatImageType2D> (input,"input.nii");

	// Rescale input
	typedef itk::RescaleIntensityImageFilter<FloatImageType2D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(input);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);

	// Setup texture params
	unsigned int radius = 1;

	std::vector<FloatImageType2D::OffsetType> offsets;
	offsets.resize(2);

	offsets[0][0] = 0; offsets[0][1] = 1;
	offsets[1][0] = 1; offsets[1][1] = 0;
	
	// Create texture filter
	typedef otb::ScalarImageToTexturesFilter<FloatImageType2D,FloatImageType2D> ScalarImageToTexturesFilterType;
	ScalarImageToTexturesFilterType::Pointer textureFilter = ScalarImageToTexturesFilterType::New();
	textureFilter->SetInput(rescaler->GetOutput());

	FloatImageType2D::SizeType rad;
	rad.Fill(radius);

	textureFilter->SetRadius(rad);
	textureFilter->SetInputImageMinimum(0);
	textureFilter->SetInputImageMaximum(255);
	textureFilter->SetNumberOfBinsPerAxis(8);

	textureFilter->SetNumberOfThreads(1);
	std::cout << "Number of threads: " << textureFilter->GetNumberOfThreads() << std::endl;

	// Allocate outputs
	std::vector<FloatImageType2D::Pointer> outputs;
	outputs.resize(8);

	for (unsigned int i=0; i<8; i++)
	{
		outputs[i] = FloatImageType2D::New();
		outputs[i]->SetRegions(input->GetLargestPossibleRegion());
		outputs[i]->SetSpacing(input->GetSpacing());
		outputs[i]->CopyInformation(input);
		outputs[i]->Allocate();
		outputs[i]->FillBuffer(0);
	}	

	// Sum texture outputs across all offsets
	for (unsigned int i=0; i<offsets.size(); i++)
	{
		textureFilter->SetOffset(offsets[i]);
		textureFilter->Update();

		for (unsigned int j=0; j<8; j++)
		{
			outputs[j] = Add(outputs[j],textureFilter->GetOutput(j));
		}
	}

	// Average outputs and rescale
	typedef itk::RescaleIntensityImageFilter<FloatImageType2D,ByteImageType2D> RescalerToByteType;
	RescalerToByteType::Pointer rescalerToByte = RescalerToByteType::New();
	rescalerToByte->SetOutputMaximum(255);
	rescalerToByte->SetOutputMinimum(0);	

	for (unsigned int j=0; j<8; j++)
	{
		outputs[j] = DivideByConstant(outputs[j],offsets.size());

		// rescale output
		rescalerToByte->SetInput(outputs[j]);
		
		std::stringstream ss;
		ss << "texture" << j << ".nii";

		WriteITK <ByteImageType2D> (rescalerToByte->GetOutput(),ss.str());
		
	}



	/*for (unsigned int i=0; i<8; i++)
	{
		std::stringstream ss;
		ss << "texture" << i << ".nii";

		typedef itk::RescaleIntensityImageFilter<FloatImageType2D,ByteImageType2D> RescalerType;
		RescalerType::Pointer rescaler = RescalerType::New();
		rescaler->SetInput(textureFilter->GetOutput(i));
		rescaler->SetOutputMaximum(255);
		rescaler->SetOutputMinimum(0);
		rescaler->Update();

		WriteITK <ByteImageType2D> (rescaler->GetOutput(),ss.str());
	}*/

	std::cout << "Time elapsed: " << ( (double) clock() - init ) / CLOCKS_PER_SEC << std::endl;

	system("pause"); 							
	return 0; 									
} 												
