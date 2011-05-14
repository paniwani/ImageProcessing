#include <itkImage.h>
#include <utils.h>
#include <otbScalarImageToTexturesFilter.h>

int main()
{
	ImageType2D::Pointer input = ReadITK <ImageType2D> ("c:/imagedata/Modified_mr10_092_13p.i0344_85_100_40_slice13.hdr");

	WriteITK <ImageType2D> (input,"input.nii");

	ImageType2D::Pointer mask = ImageType2D::New();
	mask->SetRegions(input->GetLargestPossibleRegion());
	mask->SetSpacing(input->GetSpacing());
	mask->Allocate();
	mask->FillBuffer(0);
	
	ImageType2D::IndexType maskIndex;
	maskIndex[0] = 100;
	maskIndex[1] = 100;

	ImageType2D::SizeType maskSize;
	maskSize[0] = 300;
	maskSize[1] = 300;

	ImageType2D::RegionType maskRegion;
	maskRegion.SetIndex(maskIndex);
	maskRegion.SetSize(maskSize);
	
	IteratorType2D mask_iter(mask,maskRegion);

	for (mask_iter.GoToBegin(); !mask_iter.IsAtEnd(); ++mask_iter)
	{
		mask_iter.Set(1);
	}
	
	typedef otb::ScalarImageToTexturesFilter<ImageType2D,ImageType2D,ImageType2D> ScalarImageToTexturesFilterType;
	ScalarImageToTexturesFilterType::Pointer textureFilter = ScalarImageToTexturesFilterType::New();

	textureFilter->SetInput(input);
	textureFilter->SetMaskImage(mask);

	typedef ImageType2D::SizeType SizeType;
	SizeType sradius;
	sradius.Fill(1);

	textureFilter->SetRadius(sradius);

	typedef ImageType2D::OffsetType OffsetType;
	OffsetType offset;
	offset[0] =  1;
	offset[1] =  0;

	textureFilter->SetOffset(offset);
	textureFilter->SetFeature("Entropy");

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

	WriteITK <ImageType2D> (textureFilter->GetOutput(0),"entropy.nii");

	
	return 0;
}