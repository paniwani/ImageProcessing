#include <itkImage.h>
#include <utils.h>
#include <otbScalarImageToTexturesFilter.h>

int main()
{
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("c:/imagedata/Modified_mr10_092_13p.i0344_85_100_40_slice13.hdr");
	
	typedef otb::ScalarImageToTexturesFilter<ImageType2D,ImageType2D> ScalarImageToTexturesFilterType;
	ScalarImageToTexturesFilterType::Pointer textureFilter = ScalarImageToTexturesFilterType::New();

	textureFilter->SetInput(input);

	typedef ImageType2D::SizeType SizeType;
	SizeType sradius;
	sradius.Fill(1);

	textureFilter->SetRadius(sradius);

	typedef ImageType2D::OffsetType OffsetType;
	OffsetType offset;
	offset[0] =  1;
	offset[1] =  0;

	textureFilter->SetOffset(offset);

	//textureFilter->SetInputImageMinimum(0);
	//textureFilter->SetInputImageMaximum(255);

	try
	{
		textureFilter->Update();
	} catch ( itk::ExceptionObject & excp)
	{
		std::cerr << "Error: " << std::endl;
		std::cerr << excp << std::endl;
		system("pause");
		return EXIT_FAILURE;
	}

	WriteITK <ImageType2D> (textureFilter->GetOutput(0),"energy.nii");

	
	return 0;
}