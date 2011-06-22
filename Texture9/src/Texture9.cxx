const unsigned int Dimension = 2;
#include <itkImage.h> 						
#include <iostream>
#include <utils2.h>
#include <itkColonSegmentationFilter.h>
#include <itkConfidenceConnectedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkNoiseImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBlackTopHatImageFilter.h>
#include <itkConvolutionImageFilter.h>
#include <otbScalarImageToTexturesFilter.h>
#include <otbScalarImageToTexturesMaskFilter.h>
#include <itkScalarImageToHistogramGenerator.h>

ByteImageType::Pointer Range(ByteImageType::Pointer &input, unsigned int radius);
FloatImageType::Pointer LocalEntropy(ByteImageType::Pointer &input, unsigned int radius, unsigned int numBins=256);
 												
int main(int argc, char * argv[])				
{ 					
	if( argc < 2 )
    {
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " inputImage radius" << std::endl;
		system("pause");
		return 1;
    }

	// Load image
	ImageType::Pointer inputOriginal = ReadITK <ImageType> (argv[1]);
	
	// Convert to unsigned char [0-255]
	ByteImageType::Pointer input = Rescale <ImageType,ByteImageType> (inputOriginal,0,255);
	
	ByteImageType::RegionType region = input->GetLargestPossibleRegion();
	WriteITK <ByteImageType> (input,"input.nii");

	unsigned int Radius = atoi( argv[2] );

	std::stringstream ss;
	
	/*
	// Segment colon
	typedef itk::ColonSegmentationFilter<ByteImageType,ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput(input);
	colonSegmenter->SetRemoveBoneLung(false);
	colonSegmenter->Update();
	ByteImageType::Pointer colon = colonSegmenter->GetOutput();

	ByteImageType::RegionType colonRegion = BinaryCrop <ByteImageType> (colon);
	input = CropByRegion <ByteImageType> (input,colonRegion);
	Mask <ImageType,ByteImageType> (input,colon,-1025);
	WriteITK <ImageType> (input,"input.nii");
	WriteITK <ByteImageType> (colon,"colon.nii");

	region = input->GetLargestPossibleRegion();
	ByteImageType::SizeType size = region.GetSize();
	*/

	//// Get gradient magnitude
	//typedef itk::GradientMagnitudeImageFilter<ByteImageType,FloatImageType> GradientMagnitudeFilterType;
	//GradientMagnitudeFilterType::Pointer gmFilter = GradientMagnitudeFilterType::New();
	//gmFilter->SetInput(input);
	//gmFilter->Update();
	//WriteITK <FloatImageType> (gmFilter->GetOutput(),"gm.nii");

	// Get standard deviation
	typedef itk::NoiseImageFilter<ByteImageType,FloatImageType> NoiseFilterType;
	NoiseFilterType::Pointer noiseFilter = NoiseFilterType::New();
	noiseFilter->SetInput(input);
	
	ByteImageType::SizeType radiusSize;
	radiusSize.Fill(Radius);

	noiseFilter->SetRadius(radiusSize);
	noiseFilter->Update();
	
	ss.str("");
	ss << "noise_radius_" << Radius << ".nii";
	WriteITK <FloatImageType> (noiseFilter->GetOutput(),ss.str());

	// Get range
	ByteImageType::Pointer range = Range(input,Radius);
	ss.str("");
	ss << "range_radius_" << Radius << ".nii";
	WriteITK <ByteImageType> (range,ss.str());

	// Get entropy
	unsigned int numBins = 256;

	FloatImageType::Pointer entropy = LocalEntropy( input, Radius, numBins);
	ss.str("");
	ss << "entropy_radius_" << Radius << "_numBins_" << numBins << ".nii";
	WriteITK <FloatImageType> (entropy,ss.str());

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Haralick Textures
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//// Setup 8 outputs
	//std::vector<FloatImageType::Pointer> outVector;
	//outVector.resize(8);

	//for (int i=0; i<8; i++)
	//{
	//	outVector[i] = FloatImageType::New();
	//	outVector[i]->SetRegions(region);
	//	outVector[i]->CopyInformation(input);
	//	outVector[i]->Allocate();
	//	outVector[i]->FillBuffer(0);
	//}

	//// Setup offsets
	//typedef ByteImageType::OffsetType OffsetType;

	//std::vector<OffsetType> offsetVector;
	//offsetVector.resize(2);
	//offsetVector[0][0] = 0; offsetVector[0][1] = 1;
	//offsetVector[1][0] = 1; offsetVector[1][1] = 0;

	//// Setup texture filter
	//typedef otb::ScalarImageToTexturesFilter<FloatImageType,FloatImageType> TextureFilterType;
	//TextureFilterType::Pointer textureFilter = TextureFilterType::New();
	//textureFilter->SetInput( Rescale <ByteImageType,FloatImageType> (input,0,255) );
	//textureFilter->SetInputImageMinimum(0);
	//textureFilter->SetInputImageMaximum(255);
	//textureFilter->SetNumberOfThreads(1);
	////textureFilter->SetMaskImage( Cast <ByteImageType,FloatImageType> (colon) );

	//ByteImageType::SizeType radius;
	//radius.Fill(Radius);
	//textureFilter->SetRadius(radius);

	//for (int offcount = 0; offcount < 2; offcount++)
	//{
	//	// Get textures for a particular offset
	//	textureFilter->SetOffset(offsetVector[offcount]);
	//	textureFilter->Update();
	//	
	//	// Sum textures across offsets
	//	for (int i=0; i<8; i++)
	//	{
	//		FloatImageType::Pointer texture = textureFilter->GetOutput(i);		
	//		outVector[i] = Add <FloatImageType> (outVector[i],texture);

	//		std::stringstream ss;
	//		ss << "texture" << i << "_offset_" << offcount << ".nii";
	//		//WriteITK <FloatImageType> (texture,ss.str());

	//	}
	//}

	//// Average offsets
	//for (int i=0; i<8; i++)
	//{
	//	DivideByConstant <FloatImageType> (outVector[i],offsetVector.size());

	//	std::stringstream ss;
	//	ss << "textureAveragedOffset" << i << "_radius_" << Radius << ".nii";
	//	WriteITK <ByteImageType> ( Rescale <FloatImageType,ByteImageType> (outVector[i],0,255) ,ss.str());
	//}



	return 0;
}

ByteImageType::Pointer Range(ByteImageType::Pointer &input, unsigned int radius)
{
	// get region
	ByteImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ByteImageType::SizeType rad;
	rad.Fill(radius);

	// setup output image
	ByteImageType::Pointer out = AllocateImage <ByteImageType,ByteImageType> (input);
	ByteIteratorType outIt(out,region);
	
	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt)
	{
		PixelType max = itk::NumericTraits<PixelType>::NonpositiveMin();
		PixelType min = itk::NumericTraits<PixelType>::max();
		
		for (int i=0; i<nIt.Size(); i++)
		{
			PixelType val = nIt.GetPixel(i);

			if (val > max)
				max = val;

			if (val < min)
				min = val;
		}

		outIt.Set(max-min);
	}

	return out;
}

FloatImageType::Pointer LocalEntropy(ByteImageType::Pointer &input, unsigned int radius, unsigned int numBins)
{
	// note: use unsigned char pixel type

	// get region
	ByteImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ByteImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// allocate local copy image
	ByteImageType::Pointer localImage = ByteImageType::New();
	ByteImageType::RegionType localRegion;
	
	ByteImageType::IndexType localIndex;
	ByteImageType::SizeType localSize;

	for (int i=0; i<Dimension; i++)
	{
		localIndex[i] = 0;
		localSize[i] = 2*radius + 1;
	}

	localRegion.SetIndex(localIndex);
	localRegion.SetSize(localSize);

	localImage->SetRegions(localRegion);
	localImage->Allocate();
	localImage->FillBuffer(0);

	// setup image to histogram filter
	typedef itk::Statistics::ScalarImageToHistogramGenerator< ByteImageType > HistogramGeneratorType;
	typedef HistogramGeneratorType::HistogramType HistogramType;

	HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
	histogramGenerator->SetNumberOfBins(numBins);
	histogramGenerator->SetMarginalScale(10.0); //?

	// iterate images
	FloatIteratorType outIt(out,region);
	
	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt)
	{
		// copy neighborhood values into local image
		localImage->FillBuffer(0);
		ByteIteratorType lIt(localImage,localImage->GetLargestPossibleRegion());
		
		lIt.GoToBegin();
		int i=0;

		while (!lIt.IsAtEnd())
		{
			lIt.Set( nIt.GetPixel(i) );
			i++;
			++lIt;
		}

		// convert local image to histogram
		histogramGenerator->SetInput( localImage );
		histogramGenerator->Compute();
		
		const HistogramType * histogram = histogramGenerator->GetOutput();

		// get entropy from histogram
		HistogramType::ConstIterator itr = histogram->Begin();
		HistogramType::ConstIterator end = histogram->End();

		double Sum = histogram->GetTotalFrequency();

		double Entropy = 0.0;

		while( itr != end )
		{
			const double probability = itr.GetFrequency() / Sum;

			if( probability > 0.99 / Sum )
			{
				Entropy += - probability * vcl_log( probability ) / vcl_log( 2.0 );
			}
			++itr;
		}

		// save entropy to output
		outIt.Set( (float) Entropy );
	}

	return out;
}

FloatImageType::Pointer RescaledRange(ByteImageType::Pointer &input, unsigned int radius)
{

	ByteImageType::RegionType region = input->GetLargestPossibleRegion();

	// Make map image
	ByteImageType::Pointer map = ByteImageType::New();
	
	ByteImageType::RegionType mapRegion;
	
	ByteImageType::IndexType mapIndex;
	ByteImageType::SizeType mapSize;

	for (int i=0; i<Dimension; i++)
	{
		mapIndex[i] = 0;
		mapSize[i] = 2*radius+1;
	}

	mapRegion.SetIndex(mapIndex);
	mapRegion.SetSize(mapSize);

	map->SetRegions(mapRegion);

	map->Allocate();
	map->FillBuffer(0);

	ByteIteratorType mapIt(map,mapRegion);

	// Setup neighborhood and fill map with distance squared
	typedef itk::Neighborhood<BytePixelType,Dimension> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(radius);

	typedef ByteImageType::OffsetType OffsetType;
	
	int count=0;

	for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
	{
		OffsetType off = hood.GetOffset(count);
		
		PixelType d = 0;

		for (int i=0; i<3; i++)
		{
			d += off[i]*off[i];
		}

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

	WriteITK <ByteImageType> (map,"map.nii");

	// Convert map to label map to get number of distance classes
	typedef itk::LabelImageToLabelMapFilter<ByteImageType> LabelImageToLabelMapFilterType;
	typedef LabelImageToLabelMapFilterType::OutputImageType LabelMapType;
	typedef LabelImageToLabelMapFilterType::LabelObjectType LabelObjectType;
	typedef LabelObjectType::LabelType LabelType;
	
	LabelImageToLabelMapFilterType::Pointer converter = LabelImageToLabelMapFilterType::New();
	converter->SetInput(map);
	converter->Update();
	LabelMapType::Pointer labelMap = converter->GetOutput();
	std::vector< LabelType > labelVector = labelMap->GetLabels();
	
	// Remove bkg 0 label
	if (labelVector[0] == 0)
		labelVector.erase(labelVector.begin());

	unsigned int numClasses = labelVector.size(); //ignore bkg label
	
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
	
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,ByteImageType> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToImageFilter->SetInput(relabeler->GetOutput());
	labelMapToImageFilter->Update();

	map = labelMapToImageFilter->GetOutput();
	mapIt = IteratorType(map,mapRegion);

	labelMap.~SmartPointer();

	Write(map,"map2.nii");

	// Allocate vectors to store min and max for each class
	std::vector<BytePixelType> minVector;
	std::vector<BytePixelType> maxVector;

	minVector.resize(numClasses);
	maxVector.resize(numClasses);

	for (int i=0; i<numClasses; i++)
	{
		minVector[i] = ( itk::NumericTraits<BytePixelType>::max() );
		maxVector[i] = ( itk::NumericTraits<BytePixelType>::NonpositiveMin() );
	}

	//std::cout << "Number of classes: " << numClasses << std::endl;
	//std::cout << "rangeVector.size(): " << rangeVector.size() << std::endl;

	// Allocate vector to store final range values
	std::vector<float> rv;
	rv.resize(numClasses);

	//std::cout << "rv.size() " << rv.size() << std::endl;

	// Allocate output image
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetRegions(region);
	out->CopyInformation(input);
	out->Allocate();
	out->FillBuffer(0);
	FloatIteratorType oit(out,region);

	// Iterate through image
	ByteImageType::SizeType rad;
	rad.Fill(radius);
	
	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(rad,input,region);

	for (nit.GoToBegin(), oit.GoToBegin(); !nit.IsAtEnd(); ++nit, ++oit)
	{

		// Reset min and max vectors
		for (int i=0; i<numClasses; i++)
		{
			minVector[i] = itk::NumericTraits<BytePixelType>::max();
			maxVector[i] = itk::NumericTraits<BytePixelType>::NonpositiveMin();

			rv[i] = 0;
		}

		// Iterate neighborhood and map image together
		int j=0;

		for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
		{
			if ( mapIt.Get() > 0 )
			{
				BytePixelType val = nit.GetPixel(j);
				BytePixelType m = mapIt.Get() - 1;

				if ( val < minVector[m] )
					minVector[m] = val;

				if ( val > maxVector[m] )
					maxVector[m] = val;
			}

			j++;
		}

		// Compute range for each distance class
		for (int i=0; i<numClasses; i++)
		{
			rv[i] = (float) maxVector[i] - (float) minVector[i];
		}

		// Get slope
		float slope = ComputeLogSlope(dv,rv);

		oit.Set(slope);
	}

	/*std::stringstream ss;
	ss << "rr" << radius << ".nii";
	WriteITK <FloatImageType2D> (out,ss.str());*/
							
	return out; 									
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
