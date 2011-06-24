const unsigned int Dimension = 3;
#include <itkImage.h> 						
#include <iostream> 			
#include <utils2.h>
#include <itkMaximumProjectionImageFilter.h>
#include <itkShapedNeighborhoodIterator.h>
#include <itkNeighborhood.h>
#include <algorithm>
#include <numeric>
#include <vector>

bool compare (PixelType i,PixelType j);
 												
int main(int argc, char * argv[])				
{ 		
	if( argc < 3 ) 
    { 
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " dicomDirectory radius" << std::endl;
		return EXIT_FAILURE;
    }

	ImageType::Pointer input = ReadDicom <ImageType> ( argv[1] );
	WriteITK <ImageType> (input,"input.nii");

	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// get global mip 
	/*typedef itk::MaximumProjectionImageFilter<ImageType,ImageType> MIPFilterType;
	MIPFilterType::Pointer mipF = MIPFilterType::New();
	mipF->SetInput(input);
	mipF->Update();
	WriteITK <ImageType> (mipF->GetOutput(),"globalMip.nii");*/

	// get local mip
	unsigned int Radius = atoi( argv[2] );
	
	unsigned int projectionDimension = 2;
	
	typedef itk::ShapedNeighborhoodIterator<ImageType> ShapedNeighborhoodIteratorType;
	typedef ImageType::OffsetType OffsetType;

	// make radius
	ImageType::SizeType radiusSize;
	radiusSize.Fill(Radius);

	ShapedNeighborhoodIteratorType it(radiusSize,input,region);
	it.ClearActiveList();
	
	typedef itk::Neighborhood<PixelType,Dimension> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(Radius);

	// activate offsets that move only in the projection dimension
	for (unsigned int i=0; i<hood.Size(); i++)
	{
		OffsetType offset = hood.GetOffset(i);

		bool activate = true;

		for (int j=0; j<Dimension; j++)
		{
			if ( j != projectionDimension ) 
			{
				if ( offset[j] != 0 )
				{
					activate = false;
					break;
				}
			}
		}

		if (activate)
		{
			it.ActivateOffset(offset);
		}
	}

	unsigned int numElements = it.GetActiveIndexListSize();

	std::cout << "Number of active offsets: " << numElements << std::endl;
	

	// Get max, mean, and median intensity
	ImageType::Pointer maxOut = AllocateImage<ImageType,ImageType> (input);
	ImageType::Pointer minOut = AllocateImage<ImageType,ImageType> (input);
	ImageType::Pointer meanOut = AllocateImage<ImageType,ImageType> (input);
	ImageType::Pointer medianOut = AllocateImage<ImageType,ImageType> (input);	
	
	IteratorType maxIt(maxOut,region);
	IteratorType minIt(minOut,region);
	IteratorType meanIt(meanOut,region);
	IteratorType medianIt(medianOut,region);

	for (it.GoToBegin(), maxIt.GoToBegin(), minIt.GoToBegin(), meanIt.GoToBegin(), medianIt.GoToBegin(); !it.IsAtEnd(); ++it, ++maxIt, ++minIt, ++meanIt, ++medianIt)
	{
		ShapedNeighborhoodIteratorType::ConstIterator ci;

		// Insert values into vector
		std::vector<PixelType> vec;
		//vec.resize( numElements );
		
		for (ci = it.Begin(); ci != it.End(); ci++)
		{
			vec.push_back(ci.Get());
		}

		// get stats
		sort(vec.begin(),vec.end(),compare);

		maxIt.Set( *max_element(vec.begin(),vec.end()) );
		minIt.Set( *min_element(vec.begin(),vec.end()) );
		meanIt.Set( (PixelType) ( (float) std::accumulate(vec.begin(),vec.end(),0) / (float) vec.size() ) );
		medianIt.Set( vec[ Radius ] );
	}

	std::stringstream ss;
	ss << "localMax_radius_" << Radius << ".nii";
	WriteITK <ImageType> (maxOut,ss.str());

	ss.str("");
	ss << "localMin_radius_" << Radius << ".nii";
	WriteITK <ImageType> (minOut,ss.str());

	ss.str("");
	ss << "localMean_radius_" << Radius << ".nii";
	WriteITK <ImageType> (meanOut,ss.str());

	ss.str("");
	ss << "localMedian_radius_" << Radius << ".nii";
	WriteITK <ImageType> (medianOut,ss.str());
					
	return 0; 									
} 												

bool compare (PixelType i,PixelType j) { return (i<j); }
