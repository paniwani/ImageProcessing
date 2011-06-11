#include <itkImage.h> 						
#include <iostream> 
#include <utils.h>
#include <itkNeighborhood.h>
#include <itkNeighborhoodIterator.h>
#include <itkLabelImageToLabelMapFilter.h>
//#include <itkShapedNeighborhoodIterator.h>

float ComputeSlope(std::vector<float> x, std::vector<float> y);

int main(int argc, char * argv[])				
{ 												
	// Load image
	ByteImageType2D::Pointer input = ReadITK <ByteImageType2D> ("att.png");
	ByteImageType2D::RegionType region = input->GetLargestPossibleRegion();

	int radius = 3;

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

		if ( d <= 10 )
		{
			mapIt.Set(d);
		} else {
			mapIt.Set(0);
		}
		count++;
	}

	WriteITK <ByteImageType2D> (map,"map.nii");

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

	labelMap.~SmartPointer();
	
	// Make distance vector
	std::vector<float> dv;
	dv.resize(numClasses);
	for (int i=0; i<numClasses; i++)
	{
		dv[i] = log( sqrt( (float) labelVector[i] ) );
	}

	// Rescale map
	typedef itk::RescaleIntensityImageFilter<ByteImageType2D,ByteImageType2D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(map);
	rescaler->SetOutputMaximum(numClasses);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	map = rescaler->GetOutput();
	mapIt = ByteIteratorType2D(map,mapRegion);

	WriteITK <ByteImageType2D> (map,"map2.nii");

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

	// Allocate vector to store final range values
	std::vector<float> rv;
	rv.resize(numClasses);

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
		for (int i=0; i<rangeVector.size(); i++)
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
		for (int i=0; i<rangeVector.size(); i++)
		{
			rv[i] = log( (float) rangeVector[i].max - (float) rangeVector[i].min );
		}

		// Get slope
		float slope = ComputeSlope(dv,rv);

		oit.Set(slope);
	}

	WriteITK <FloatImageType2D> (out,"out.nii");
	
	system("pause"); 							
	return 0; 									
} 							

float ComputeSlope(std::vector<float> x, std::vector<float> y)
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

	// check that range is not all 0
	bool isZero = true;
	for (int i=0; i<N; i++)
	{
		if (y[i] != 0)
		{
			isZero = false;
			break;
		}
	}

	if (!isZero)
	{
		for (int i=0; i<N; i++)
		{
			s1 += x[i];
			s2 += x[i]*x[i];
			t0 += y[i];
			t1 += x[i]*y[i];
		}

		return ( s0*t1 - s1*t0 ) / ( s0*s2 - s1*s1 );
	} else {
		return 0;
	}
}
