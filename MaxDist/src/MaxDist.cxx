#include <itkImage.h> 						
#include <iostream>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkRandomImageSource.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileWriter.h>
 												
int main(int argc, char * argv[])				
{ 							
	typedef itk::Image<unsigned char,2> ImageType;
	typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

	// Make random binary image from [0,255]
	typedef itk::RandomImageSource<ImageType> RandomImageSourceType;
	RandomImageSourceType::Pointer rs = RandomImageSourceType::New();
	rs->SetMin(0);
	rs->SetMax(255);

	ImageType::SizeType size;
	size[0] = 256; size[1] = 256;

	rs->SetSize(size);
	rs->Update();

	// Threshold image
	typedef itk::BinaryThresholdImageFilter<ImageType,ImageType> ThesholdType;
	ThesholdType::Pointer thresholder = ThesholdType::New();
	thresholder->SetInput( rs->GetOutput() );
	thresholder->SetOutsideValue(255);
	thresholder->SetInsideValue(0);
	thresholder->SetLowerThreshold(128);
	thresholder->Update();

	ImageType::Pointer input = thresholder->GetOutput();
		
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Write input
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(input);
	writer->SetFileName("input.png");
	writer->Update();

	// Get max distance
	IteratorType it(input,region);
	IteratorType ot(input,region);

	float maxDist = 0;

	int count = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() == 255)
		{
		
			ImageType::IndexType idx1 = it.GetIndex();
			
			ot.GoToBegin();
			
			for (int i=0; i<count; i++)
			{
				++ot;
			}
			
			while (!ot.IsAtEnd())
			{
				
				if (ot.Get() == 255)
				{
					
					ImageType::IndexType idx2 = ot.GetIndex();
					
					float d = 0;
					
					for (int j=0; j<2; j++)
					{
						d += (idx2[j] - idx1[j])*(idx2[j] - idx1[j]);
					}
					
					d = sqrt(d);
					
					if (d > maxDist)
						maxDist = d;
				}
				
				++ot;
			}
		}
		
		count++;
	}

	std::cout << "Max distance: " << maxDist << std::endl;
	system("pause");
							
	return 0; 									
} 												
