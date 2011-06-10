#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <itkShapedNeighborhoodIterator.h>
 												
int main(int argc, char * argv[])				
{ 												
	// Load image
	ByteImageType2D::Pointer input = ReadITK <ByteImageType2D> ("ADS40RoiSmall.png");
	ByteImageType2D::RegionType region = input->GetLargestPossibleRegion();

	// Setup shape iterator
	int radius = 3;

	ByteImageType2D::SizeType rad;
	rad.Fill(radius);
	
	typedef itk::ShapedNeighborhoodIterator<ByteImageType2D> ShapedNeighborhoodIteratorType;
	ShapedNeighborhoodIteratorType it(rad,input,region);
 
	for (float y = -radius; y <= radius; y++)
	{
		for (float x = -radius; x <= radius; x++)
		{
			ShapedNeighborhoodIteratorType::OffsetType off;
			float disSquared = x*x + y*y;
			if (disSquared <= 10)
			{
				off[0] = static_cast<int>(x);
				off[1] = static_cast<int>(y);
				it.ActivateOffset(off);
			}
		}
	}

	// count number of elements
	unsigned int num = 0;
	ShapedNeighborhoodIteratorType::ConstIterator ci;
	for (ci = it.Begin(); ci != it.End(); ci++)
	{
		num++;
	}

	std::vector<int> dist;
	dist.resize(num);

	for (int i=0; i < num ; i++)
	{	
		ShapedNeighborhoodIteratorType::OffsetType off = it.GetOffset(i);
		dist[i] = off[0]*off[0] + off[1]*off[1];
	}



	std::cout << num << std::endl;

	//ShapedNeighborhoodIteratorType::OffsetType off = it.GetOffset(0);



	



	
	
	system("pause"); 							
	return 0; 									
} 												
