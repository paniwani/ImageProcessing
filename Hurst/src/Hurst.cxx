#include "itkImage.h"
#include <iostream>
#include <itkNeighborhoodIterator.h>


#include <utils.h>

struct pos {
	float distance;
	float bright;
	float dark;
};

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/wr1_019_10s.i0462/dcm");
	WriteITK <ImageType> (input, "input.nii");
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	typedef ikt::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 3;
	radius[1] = 3;

	NeighborhoodIteratorType it(radius, input, region);

	// Calculate distance for each position in neighborhood
	ImageType::SizeType size = it.GetSize();
	int num = size[0]*size[1];

	float distance[ num ];

	for (int i=0; i < num; i++)
	{
		ImageType::OffsetType offset = it->GetOffset( i );
		float dist = sqrt( offset[0]*offset[0] + offset[1]*offset[1] );

		distance[i] = dist;
	}

	// Create vector of distance/intensity struct

	std::vector< pos > pos_array;

	for (int i=0; i < (num-1)/2; i++)
	{
		bool new_dist = true;

		for (int j=0; j < pos_array.size ; j++)
		{
			if ( distance[i] == pos_array[i].distance )
			{
				new_dist = false;
				break;
			}
		}
		
		if (new_dist)
		{
			pos p;
			p.distance = distance[i];
			p.bright = ImageType::PixelType::Max;
			p.dark = 0;
			pos_array.push_back( p );
		}
	}

	// Iterator through image

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		for (int i=0; i < num; i++)
		{
			float val = it.GetPixel(i);
			
			for (int j=0; j < pos_array.size ; j++)
			{
				if ( distance[i] == pos_array[i].distance )
				{
					p = pos_array[i];

					if ( val > p.bright


					break;	
				}
			}
		}
			
		}
	}






	




	system("pause");
	return 0;
}