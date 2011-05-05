#include <itkImage.h>
#include <utils.h>
#include <itkGradientImageFilter.h>
#include <itkCovariantVector.h>

int main()
{
	// Create fake image
	ImageType2D::Pointer input = ImageType2D::New();

	ImageType2D::IndexType index;
	index[0] = 0;
	index[1] = 0;

	ImageType2D::SizeType size;
	size[0] = 3;
	size[1] = 3;

	ImageType2D::RegionType region;
	region.SetIndex( index );
	region.SetSize( size );

	ImageType2D::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	

	input->SetDirection( direction );
	input->SetRegions( region );
	input->Allocate();

	input->FillBuffer(0);

	ImageType2D::DirectionType d = input->GetDirection();
	std::cout << "Direction: " << d[0][0] << " " << d[1][1] << std::endl;

	IteratorType2D input_iter(input,region);
	
	int count = 0;

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		input_iter.Set( count++ );
	}

	input->SetPixel(index,100);

	WriteITK <ImageType2D> (input,"input.nii");

	typedef itk::GradientImageFilter<ImageType2D> GradientImageFilterType;
	GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
	gradientFilter->SetInput( input );
	//gradientFilter->SetUseImageSpacingOff();
	//gradientFilter->SetUseImageDirection(false);
	gradientFilter->Update();

	typedef itk::CovariantVector<float,2> VectorType;
	typedef itk::Image<VectorType,2> VectorImageType;
	typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;

	VectorImageType::Pointer grad = gradientFilter->GetOutput();

	VectorIteratorType grad_iter(grad,region);

	for (input_iter.GoToBegin(), grad_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		ImageType2D::IndexType idx = input_iter.GetIndex();

		VectorType g = grad_iter.Get();

		std::cout << idx[0] << "\t" << idx[1] << "\t" << input_iter.Get() << "\t" << g[0] << "\t" << g[1] << std::endl;
	}

	typedef itk::ImageFileWriter<VectorImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(grad);
	writer->SetFileName("grad.mhd");
	writer->Update();

	system("pause");
	return 0;
}
