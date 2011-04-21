#include "itkImage.h"
#include <iostream>
#include <itkNeighborhoodIterator.h>
#include <utils.h>

struct pos {
	float distance;
	float bright;
	float dark;
	float range;
};

bool compDistance (pos a, pos b);

int main()
{
	ImageType::Pointer input = ReadDicom <ImageType> ("C:/ImageData/mr9/mr4_065_13p.i0375/dcm");
	WriteITK <ImageType> (input, "input.nii");
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(0);
	radius[0] = 1;
	radius[1] = 1;

	NeighborhoodIteratorType it(radius, input, region);

	// Calculate distance for each position in neighborhood
	ImageType::SizeType size = it.GetSize();
	int num = size[0]*size[1];

	float *distance;
	distance = new float[num];
	//float distance[ num ];

	for (int i=0; i < num; i++)
	{
		ImageType::OffsetType offset = it.GetOffset( i );
		float dist = sqrt( (float) offset[0]*offset[0] + (float) offset[1]*offset[1] );

		distance[i] = dist;
	}

	for (int i=0; i < num; i++)
	{
		std::cout << distance[i] << std::endl;
	}

	// Create vector of distance/intensity struct

	std::vector< pos > pos_array;

	for (int i=0; i < (num-1)/2; i++)
	{
		bool new_dist = true;

		for (int j=0; j < pos_array.size() ; j++)
		{
			if ( distance[i] == pos_array[j].distance )
			{
				new_dist = false;
				break;
			}
		}
		
		if (new_dist)
		{
			pos p;
			p.distance = distance[i];
			p.bright = itk::NumericTraits<PixelType>::NonpositiveMin();
			p.dark = itk::NumericTraits<PixelType>::max();
			p.range = 0;
			pos_array.push_back( p );
		}
	}

	// Sort based on distance
	sort( pos_array.begin(), pos_array.end(), compDistance );

	// Allocate output
	ImageType::Pointer output = ImageType::New();
	output->SetRegions( region );
	output->SetSpacing( input->GetSpacing() );
	output->SetDirection( input->GetDirection() );
	output->CopyInformation( input );
	output->Allocate();
	output->FillBuffer(0);

	IteratorType oit(output,region);

	// Iterator through image

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		// Reset values for struct vector
		for (int j=0; j < pos_array.size() ; j++)
		{
			pos_array[j].bright = itk::NumericTraits<PixelType>::NonpositiveMin();
			pos_array[j].dark = itk::NumericTraits<PixelType>::max();
			pos_array[j].range = 0;
		}			

		// Visit all pixels in neighborhood and get bright/dark values for each distance
		for (int i=0; i < num; i++)
		{
			float val = it.GetPixel(i);
			
			for (int j=0; j < pos_array.size() ; j++)
			{
				if ( distance[i] == pos_array[j].distance )
				{
					if ( val > pos_array[j].bright )
						pos_array[j].bright = val;

					if ( val < pos_array[j].dark )
						pos_array[j].dark = val;

					break;	
				}
			}			
		}

		// Compute linear regression of ln(distance) vs. ln(range)
		float s0,s1,s2,t0,t1;
		s1 = s2 = t0 = t1 = 0;
		s0 = 1;

		// Ensure that all ranges are not 0
		bool allRangeZero = true;
		for (int j=0; j < pos_array.size() ; j++)
		{
			if ( ( pos_array[j].bright - pos_array[j].dark ) != 0 )
				allRangeZero = false;
		}

		float slope = 0;

		if ( !allRangeZero )
		{
			for (int j=0; j < pos_array.size() ; j++)
			{

				if ( ( pos_array[j].bright - pos_array[j].dark ) != 0 )
				{
					s0++;
					s1 += log( pos_array[j].distance );
					s2 += log( pos_array[j].distance ) * log( pos_array[j].distance );
					t0 += log( pos_array[j].bright - pos_array[j].dark );
					t1 += log( pos_array[j].distance) * log( pos_array[j].bright - pos_array[j].dark );
				}
			}

			slope = ( s0*t1 - s1*t0 ) / (s0*s2 - s1*s1);
		}

		oit.Set( slope );
	}

	std::stringstream ss;
	ss << "hurst_" << radius[0] << "x" << radius[1] << ".nii";

	WriteITK < ImageType > (output, ss.str() );

	delete [] distance;

	system("pause");
	return 0;
}

bool compDistance (pos a, pos b) { return (a.distance < b.distance); }