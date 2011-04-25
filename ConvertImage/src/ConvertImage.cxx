/****************************************************************************************   
 Convert to display only tissue stool in colon (unsigned int)
****************************************************************************************/

#include <itkImage.h>
#include <utils.h>
#include <itkCastImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkColonSegmentationFilter.h>

int main()
{
	typedef itk::Image<unsigned short, 3> UShortImageType;

	ImageType::Pointer input = ReadDicom <ImageType> ("c:/imagedata/mr10-uncleansed/mr10_092_13p.i0344/dcm");

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	WriteITK <ImageType> (input,"input.nii");

	typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage(input);
	minMaxCalc->ComputeMinimum();
	PixelType min = minMaxCalc->GetMinimum();

	typedef itk::ColonSegmentationFilter<ImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colon_segmenter = ColonSegmentationFilterType::New();
	colon_segmenter->SetInput(input);
	colon_segmenter->SetRemoveBoneLung(true);
	colon_segmenter->SetForegroundValue(1);
	colon_segmenter->SetBackgroundValue(0);
	colon_segmenter->Update();
	ByteImageType::Pointer colon = colon_segmenter->GetOutput();

	WriteITK <ByteImageType> (colon,"colon.nii");
	
	IteratorType input_iter(input,region);
	ByteIteratorType colon_iter(colon,region);

	for(input_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++colon_iter)
	{
		if ( colon_iter.Get() == 255 )
		{
			if ( input_iter.Get() < -300 || input_iter.Get() > 1400 )
			{
				input_iter.Set( min );
			}
		} else {
			input_iter.Set( min );
		}
	}

	WriteITK <ImageType> (input,"shrunk.nii");

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		input_iter.Set( input_iter.Get() + abs(min) );
	}

	WriteITK <ImageType> (input,"shift.nii");

	typedef itk::CastImageFilter<ImageType, UShortImageType> CastImageFilterType;
	CastImageFilterType::Pointer caster = CastImageFilterType::New();
	caster->SetInput(input);
	caster->Update();
	
	WriteITK <UShortImageType> (caster->GetOutput(), "output.nii");

	return 0;
}
