#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <itkNeighborhood.h>
#include <itkNeighborhoodIterator.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkRelabelLabelMapFilter.h>
//#include <itkShapedNeighborhoodIterator.h>

float ComputeLogSlope(std::vector<float> x, std::vector<float> y);

int main(int argc, char * argv[])				
{ 			
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << "RescaledRange inputImage radius";
		return EXIT_FAILURE;
	}

	// Load image
	ByteImageType2D::Pointer input = ReadITK <ByteImageType2D> (argv[1]);
	ByteImageType2D::RegionType region = input->GetLargestPossibleRegion();

	int radius = atoi(argv[2]);

	// Make map image
	ByteImageType2D::Pointer map = ByteImageType2D::New();
	
	ByteImageType2D::RegionType mapRegion;
	
	ByteImageType2D::IndexType mapIndex;
	mapIndex[0] = 0; mapIndex[1] = 0;

	ByteImageType2D::SizeType mapSize;
	mapSize[0] = 2*radius+1; mapSize[1] = 2*radius+1;

	mapRegion.SetIndex(mapIndex);
	mapRegion.SetSize(mapSize);

	map->SetRegions(mapRegion);

	map->Allocate();
	map->FillBuffer(0);

	ByteIteratorType2D mapIt(map,mapRegion);

	// Setup neighborhood and fill map with distance squared
	typedef itk::Neighborhood<unsigned char,2> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(radius);

	typedef ByteImageType2D::OffsetType OffsetType;
	
	int count=0;

	for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
	{
		OffsetType off = hood.GetOffset(count);
		
		unsigned char d = off[0]*off[0] + off[1]*off[1];

		// make octagonal shape
		/*if ( d <= 10 )
		{
			mapIt.Set(d);
		} else {
			mapIt.Set(0);
		}*/

		mapIt.Set(d);

		count++;
	}

	//WriteITK <ByteImageType2D> (map,"map.nii");

	// Convert map to label map to get number of distance classes
	typedef itk::LabelImageToLabelMapFilter<ByteImageType2D> LabelImageToLabelMapFilterType;
	typedef LabelImageToLabelMapFilterType::OutputImageType LabelMapType;
	typedef LabelImageToLabelMapFilterType::LabelObjectType LabelObjectType;
	typedef LabelObjectType::LabelType LabelType;
	
	LabelImageToLabelMapFilterType::Pointer converter = LabelImageToLabelMapFilterType::New();
	converter->SetInput(map);
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();
	std::vector< LabelType > labelVector = labelMap->GetLabels();
	
	unsigned int numClasses = labelVector.size();
	
	// Make distance vector
	std::vector<float> dv;
	dv.resize(numClasses);
	for (int i=0; i<numClasses; i++)
	{
		dv[i] = sqrt( (float) labelVector[i] );
	}

	// Relabel map consecutively to correspond to each distance class	
	typedef itk::RelabelLabelMapFilter<LabelMapType> RelabelLabelMapFilterType;
	RelabelLabelMapFilterType::Pointer relabeler = RelabelLabelMapFilterType::New();
	relabeler->SetInput(labelMap);
	
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,ByteImageType2D> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToImageFilter->SetInput(relabeler->GetOutput());
	labelMapToImageFilter->Update();

	map = labelMapToImageFilter->GetOutput();
	mapIt = ByteIteratorType2D(map,mapRegion);

	labelMap.~SmartPointer();

	//WriteITK <ByteImageType2D> (map,"map2.nii");

	// Allocate range vector to store min and max for each class
	struct srange {
		unsigned char min;
		unsigned char max;
	};
	
	std::vector<srange> rangeVector;
	rangeVector.resize(numClasses);

	for (int i=0; i<numClasses; i++)
	{
		struct srange r;
		r.min = itk::NumericTraits<unsigned char>::max();
		r.max = itk::NumericTraits<unsigned char>::NonpositiveMin();
		rangeVector.push_back(r);
	}

	//std::cout << "Number of classes: " << numClasses << std::endl;
	//std::cout << "rangeVector.size(): " << rangeVector.size() << std::endl;

	// Allocate vector to store final range values
	std::vector<float> rv;
	rv.resize(numClasses);

	//std::cout << "rv.size() " << rv.size() << std::endl;

	// Allocate output image
	FloatImageType2D::Pointer out = FloatImageType2D::New();
	out->SetRegions(region);
	out->CopyInformation(input);
	out->Allocate();
	out->FillBuffer(0);
	FloatIteratorType2D oit(out,region);

	// Iterate through image
	ByteImageType2D::SizeType rad;
	rad.Fill(radius);
	
	typedef itk::NeighborhoodIterator<ByteImageType2D> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(rad,input,region);

	for (nit.GoToBegin(), oit.GoToBegin(); !nit.IsAtEnd(); ++nit, ++oit)
	{

		// Reset range vector
		for (int i=0; i<numClasses; i++)
		{
			rangeVector[i].min = itk::NumericTraits<unsigned char>::max();
			rangeVector[i].max = itk::NumericTraits<unsigned char>::NonpositiveMin();

			rv[i] = 0;
		}

		// Iterate neighborhood and map image together
		int j=0;

		for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
		{
			if ( mapIt.Get() > 0 )
			{
				unsigned char val = nit.GetPixel(j);
				unsigned char m = mapIt.Get() - 1;

				if ( val < rangeVector[m].min )
					rangeVector[m].min = val;

				if ( val > rangeVector[m].max )
					rangeVector[m].max = val;
			}

			j++;
		}

		// Compute range for each distance class
		for (int i=0; i<numClasses; i++)
		{
			rv[i] = (float) rangeVector[i].max - (float) rangeVector[i].min;
		}

		// Get slope
		float slope = ComputeLogSlope(dv,rv);

		oit.Set(slope);
	}

	std::stringstream ss;
	ss << "rr" << radius << ".nii";
	WriteITK <FloatImageType2D> (out,ss.str());
							
	return 0; 									
} 							

float ComputeLogSlope(std::vector<float> x, std::vector<float> y)
{
/*
To find the
"least squares regression line" , y = Mx + B
for the points (x(i),y(i)) for i = 1,N

you would do the following :

s0 = N+1

s1 = sum_of(x(i)) for i = 1,N

s2 = sum_of(x(i)*x(i)) for i = 1,N

t0 = sum_of(y(i)) for i = 1,N

t1 = sum_of(x(i)*y(i)) for i = 1,N


then 

M = ( s0*t1 - s1*t0 ) / (s0*s2 - s1*s1)

B = ( s2*t0 - s1*t1 ) / (s0*s2 - s1*s1)

and the regression line is given by 

y = Mx + B*/
	
	unsigned int N = x.size();
	float s0,s1,s2,t0,t1;
	
	s0 = N+1;
	s1 = 0; s2=0; t0=0; t1=0;

	for (int i=0; i<N; i++)
	{
		if (y[i] > 0) // ensure non-zero range is used to compute log-slope
		{
			float lx = log(x[i]);
			float ly = log(y[i]);

			s1 += lx;
			s2 += lx*lx;
			t0 += ly;
			t1 += lx*ly;
		} else {
			s0 -= -1;
		}
	}

	// Don't divide by zero
	float num = s0*t1 - s1*t0;
	float den = s0*s2 - s1*s1;

	if (den == 0)
	{
		return 0;
	} else {
		return num/den;
	}

}
