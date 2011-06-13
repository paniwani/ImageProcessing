#include "RemoveStool4.h"
#include "utility.cxx"
#include "io.cxx"
#include "scattercorrection.cxx"
#include "QR.cxx"
#include "EM.cxx"
#include "HessianFunctions2.cxx"

void Test(ImageType::Pointer &inputOriginal, ImageType::Pointer &input, ByteImageType::Pointer &colon)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	//// get 2d grad mag
	//typedef itk::GradientMagnitudeImageFilter<ImageType2D,FloatImageType2D> GradientMagnitudeImageFilterType2D;
	//typedef itk::SliceBySliceImageFilter<ImageType,FloatImageType,GradientMagnitudeImageFilterType2D> SliceBySliceImageFilterType;

	//GradientMagnitudeImageFilterType2D::Pointer gf = GradientMagnitudeImageFilterType2D::New();

	//SliceBySliceImageFilterType::Pointer sgf = SliceBySliceImageFilterType::New();
	//sgf->SetInput(input);
	//sgf->SetFilter(gf);
	//sgf->Update();
	//FloatImageType::Pointer gm = sgf->GetOutput();

	//Write(gm,"gm2D.nii");

	//// allocate vmap
	//VoxelImageType::Pointer vmap = VoxelImageType::New();
	//vmap->SetRegions(region);
	//vmap->Allocate();
	//vmap->FillBuffer(Unclassified);

	//ApplyThresholdRules(input, region, gm, vmap, colon, 180);

	//Write(vmap,"vmap.nii");

	//// Do scatter correction

	//ImageType::Pointer inputScatter = ScatterCorrection(inputOriginal,colon);
	//Write(inputScatter,"inputScatter.nii");

	//// Show change in scatter correction
	//typedef itk::SubtractImageFilter<ImageType,ImageType> SubtracterType;
	//SubtracterType::Pointer subtracter = SubtracterType::New();
	//subtracter->SetInput1(input);
	//subtracter->SetInput2(inputScatter);
	//subtracter->Update();
	//Write(subtracter->GetOutput(),"scatterChange.nii");

	//// Get vmap with updated input
	//vmap->FillBuffer(Unclassified);
	//ApplyThresholdRules(input, region, gm, vmap, colon, 180);

	typedef itk::FastBilateralImageFilter<ImageType,ImageType> FastBilateralImageFilterType;
	FastBilateralImageFilter::Pointer bilateralFilter = FastBilateralImageFilter::New();
	bilateralFilter->SetInput(input);
	bilateralFilter->SetDomainSigma(


}

int main(int argc, char * argv[])
{
	// Start clock
	clock_t init,final;
	init = clock();

	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}

	ImageType::Pointer				inputOriginal			= ImageType::New();
	ImageType::Pointer				input					= ImageType::New();
	ByteImageType::Pointer			colon					= ByteImageType::New();
	FloatImageType::Pointer			gradientMagnitude		= FloatImageType::New();
	VoxelImageType::Pointer			vmap					= VoxelImageType::New();

	debug.open("debug.txt");

	// Load images and segment colon
	Setup(argv[1],inputOriginal,input,colon,gradientMagnitude);


	Test(inputOriginal,input,colon);


	/*for (int i=1; i<10; i++)
	{
		std::stringstream ss;
		ss << "rr" << i << ".nii";

		FloatImageType::Pointer rr = RescaledRange(input,3);
		Write(rr,ss.str());
	}*/	

	/*for (int i=1; i<5; i++)
	{
		std::stringstream ss;
		ss << "sd" << i << ".nii";

		FloatImageType::Pointer sd = StandardDeviation(input,colon,i);
		Write(sd,ss.str());
	}

	for (int i=1; i<5; i++)
	{
		std::stringstream ss;
		ss << "range" << i << ".nii";

		ImageType::Pointer range = Range(input,colon,i);
		Write(range,ss.str());
	}*/



	//TextureTest(input,colon);
	//ConnectedTest(input,gradientMagnitude,colon);
	//AdaptiveThreshold(input,colon);
	//LocalThreshold(input,colon,gradientMagnitude,vmap);

	// Initial segmentation, save tissue stool threshold
	//PixelType tst = SingleMaterialClassification(input, gradientMagnitude, vmap, colon);

	//// Apply scatter correction
	//input = ScatterCorrection(inputOriginal,colon,vmap);

	//// Update vmap with new scatter input
	//ApplyThresholdRules(input,REGION,gradientMagnitude,vmap,colon,tst);

	//Write(vmap,"vmapScatter.nii");	

	//// Determine boundary types
	//FloatImageType::Pointer smax = QuadraticRegression(input,colon,vmap,gradientMagnitude,tst);

	//gradientMagnitude.~SmartPointer();

	//ArrayImageType::Pointer partial = ComputePartials(input,vmap,colon,smax);
	//
	//FixATT(input,partial,vmap,colon,smax,tst);

	//smax.~SmartPointer();

	//// Perform hessian analysis
	//HessianAnalysis(input,colon,vmap,partial,tst);	

	//// Expectation Maximization
	//EM(partial,colon,input);

	//// Perform subtraction
	//ImageType::Pointer carstonOutput = Subtraction(input,inputOriginal,colon,partial,vmap);

	//// End clock
	//final = clock() - init;

	//std::cout << (double) final / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;

	//system("pause");
	return 0;
}



void TextureTest(ImageType::Pointer &input, ByteImageType::Pointer &colon)
{
	ImageType::Pointer range = Range(input,colon,1);
	Write(range,"range.nii");

	FloatImageType::Pointer sd = StandardDeviation(input,colon,1);
	Write(sd,"sd.nii");	
}



void ConnectedTest(ImageType::Pointer &input, FloatImageType::Pointer &gradientMagnitude, ByteImageType::Pointer &colon)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	/*smooth input
	typedef itk::FastBilateralImageFilter<ImageType,ImageType> SmoothFilterType;
	SmoothFilterType::Pointer smoother = SmoothFilterType::New();
	smoother->SetInput(input);

	double rangeSigma = 50.0;
	double domainSigma = 10.0;

	smoother->SetRangeSigma(rangeSigma);
	smoother->SetDomainSigma(domainSigma);
	smoother->Update();
	input = smoother->GetOutput();
	WriteITK <ImageType> (inputSmooth,"inputSmooth.nii");*/

	// get conservative tissue estimate
	ByteImageType::Pointer tmap = AllocateByteImage(input);

	ByteIteratorType tIt(tmap,region);
	ByteIteratorType cIt(colon,region);
	IteratorType It(input,region);
	FloatIteratorType gIt(gradientMagnitude,region);

	for (It.GoToBegin(), tIt.GoToBegin(), gIt.GoToBegin(), cIt.GoToBegin(); !It.IsAtEnd(); ++It, ++tIt, ++gIt, ++cIt)
	{
		if (cIt.Get() != 0)
		{
			PixelType I = It.Get();
			float G = gIt.Get();

			if ( I >= -250 && I <= 150 && G <= 300)
			{
				tIt.Set(255);
			}
		}
	}

	Write(tmap,"tmap.nii");

	//// keep components with large (>100 pixels) size only
	//typedef itk::BinaryShapeOpeningImageFilter<ByteImageType> OpeningType;
	//OpeningType::Pointer open = OpeningType::New();
	//open->SetInput(tmap);
	//open->SetForegroundValue(255);
	//open->SetBackgroundValue(0);
	//open->SetAttribute("Size");
	//open->SetLambda(100);
	//open->Update();
	//ByteImageType::Pointer tmapBig = open->GetOutput();
	//Write(tmapBig,"tmapBig1.nii");

	//open->SetFullyConnected(true);
	//open->Update();
	//tmapBig = open->GetOutput();
	//Write(tmapBig,"tmapBig2.nii");

	// keep only those components that are fully surrounded by tissue

	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorType;
	
	ByteImageType::SizeType radius;
	radius.Fill(1);

	NeighborhoodIteratorType nIt(radius,tmap,region);

	ByteImageType::Pointer tmap2 = AllocateByteImage(input);
	ByteIteratorType t2It(tmap2,region);

	for (nIt.GoToBegin(), t2It.GoToBegin(), cIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++t2It, ++cIt)
	{
		if (cIt.Get() != 0)
		{
			if ( nIt.GetCenterPixel() != 0)
			{
				bool surrounded = true;

				for (unsigned int i=0; i<nIt.Size(); i++)
				{
					if ( i != ((nIt.Size()-1)/2) )
					{
						if ( nIt.GetPixel(i) == 0 )
						{
							surrounded = false;
							break;
						}
					}
				}

				if (surrounded)
				{
					t2It.Set(255);
				}
			}
		}
	}

	Write(tmap2,"tmap2.nii");

	typedef itk::ConfidenceConnectedImageFilter<ImageType,ByteImageType> ConfidenceConnectedType;
	ConfidenceConnectedType::Pointer connecter = ConfidenceConnectedType::New();
	connecter->SetInput(input);
	connecter->SetInitialNeighborhoodRadius(1);
	connecter->SetMultiplier(2);
	connecter->SetNumberOfIterations(1);
	connecter->SetReplaceValue(255);

	for (t2It.GoToBegin(); !t2It.IsAtEnd(); ++t2It)
	{
		if (t2It.Get() != 0)
		{
			connecter->AddSeed( t2It.GetIndex() );
		}
	}

	connecter->Update();

	std::cout << "mean: " << connecter->GetMean() << std::endl;
	std::cout << "std: " << sqrt(connecter->GetVariance()) << std::endl;

	Write(connecter->GetOutput(),"connected.nii");

}

void AdaptiveThreshold(ImageType::Pointer &input, ByteImageType::Pointer &colon)
{
	// get an adaptive threshold between tissue and stool
	typedef itk::AdaptiveOtsuThresholdImageFilter<ImageType,ImageType> AOFilterType;
	AOFilterType::Pointer otsuFilter = AOFilterType::New();
	otsuFilter->SetOutsideValue(0);
	otsuFilter->SetInsideValue(255);
	otsuFilter->SetInput(input);

	try {
		otsuFilter->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error: " << err << std::endl;
		system("pause");
	}
	
	ImageType::Pointer bth = otsuFilter->GetOutput();
	Write(bth,"bth.nii");
	
	ImageType::Pointer th = otsuFilter->GetThresholdImage();
	Write(th,"th.nii");
}







// divide image and auto threshold I and G in each region
void LocalThreshold(ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap)
{
	// get full region
	ImageType::RegionType fullRegion = input->GetLargestPossibleRegion();

	unsigned int requestedDivisions = 128;//itk::NumericTraits<unsigned int>::max();

	// divide image into hyper cubes
	typedef itk::ImageRegionMultidimensionalSplitter<3> SplitterType;
	SplitterType::Pointer splitter = SplitterType::New();

	unsigned int divisions = splitter->GetNumberOfSplits(fullRegion,requestedDivisions);

	std::cout << "requested number of divisions: " << requestedDivisions << std::endl;
	std::cout << "number of divisions: " << divisions << std::endl;

	// allocate vmap
	vmap->SetRegions( fullRegion );
	vmap->SetSpacing( input->GetSpacing() );
	vmap->SetDirection( input->GetDirection() );
	vmap->CopyInformation( input );
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);

	// setup otsu calculator
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(input);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1500);

	// get global otsu threshold
	otsuCalculator->Compute();
	PixelType globalOtsu = otsuCalculator->GetThreshold();

	// allocate otsu map image
	FloatImageType::Pointer otsuMap = AllocateFloatImage(input);
	
	for (unsigned int i=0; i<divisions; i++)
	{
		// get local region
		ImageType::RegionType localRegion = splitter->GetSplit(i,divisions,fullRegion);

		//// validate region
		//ImageType::IndexType localIndex = localRegion.GetIndex();
		//ImageType::SizeType localSize = localRegion.GetSize();

		//for (int j=0; j<3; j++)
		//{
		//	int over = ( localIndex[j] + localSize[j] ) - fullRegion.GetSize()[j];

		//	if (over > 0)
		//	{
		//		localSize[j] -= over;
		//	}
		//}
		//localRegion.SetSize(localSize);

		// set otsu region
		otsuCalculator->SetRegion(localRegion);
		std::stringstream ss;
		ss << i << "_intensity.csv";
		otsuCalculator->SetPrintHistogram(ss.str());
		otsuCalculator->Compute();
		
		PixelType otsu = otsuCalculator->GetThreshold();
		std::cout << "otsu " << i << ": " << otsu << std::endl;

		// set bounds on local tissue stool threshold
		if ( otsu < 100 || otsu > 700 )
		{
			otsu = globalOtsu;
		}

		// fill map with otsu value
		FloatIteratorType oit(otsuMap,localRegion);
		for (oit.GoToBegin(); !oit.IsAtEnd(); ++oit)
		{
			oit.Set(i);
		}

		ApplyThresholdRules(input,localRegion,gradientMagnitude,vmap,colon,otsu);

		/* write input in local region
		std::stringstream ss;
		ss << "region" << i << ".nii";

		Write( Crop(input,localRegion) , ss.str() );	*/

		/*std::stringstream ss;
		ss << "vmap" << i << ".nii";
		Write(vmap,ss.str());*/
	}

	std::stringstream ss;
	ss << "vmapLocal_ " << divisions << ".nii";
	Write(vmap,ss.str());

	ApplyThresholdRules(input,fullRegion,gradientMagnitude,vmap,colon,globalOtsu);
	Write(vmap,"vmapGlobal.nii");

	Write(otsuMap,"otsuMap.nii");
}

void LevelSet(ImageType::Pointer &input, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// get intial contour
	// get tissue
	ByteImageType::Pointer tissue = ByteImageType::New();
	tissue->SetRegions(region);
	tissue->SetSpacing(input->GetSpacing());
	tissue->CopyInformation(input);
	tissue->Allocate();
	tissue->FillBuffer(0);
	ByteIteratorType tissueIt(tissue,region);
	ByteIteratorType colonIt(colon,region);
	VoxelIteratorType vmapIt(vmap,region);

	for (vmapIt.GoToBegin(), colonIt.GoToBegin(), tissueIt.GoToBegin(); !vmapIt.IsAtEnd(); ++vmapIt, ++tissueIt, ++colonIt)
	{
		if (vmapIt.Get() == Tissue)
			tissueIt.Set(255);

		if (colonIt.Get() == 0)
			tissueIt.Set(255);
	}

	Write(tissue,"tissue.nii");

	//// get contours
	//typedef itk::BinaryContourImageFilter<ByteImageType,ByteImageType> ContourType;
	//ContourType::Pointer contourFilter = ContourType::New();
	//contourFilter->SetInput(tissue);
	//contourFilter->SetForegroundValue(255);
	//contourFilter->SetBackgroundValue(0);
	//contourFilter->Update();
	//ByteImageType::Pointer tissueContour = contourFilter->GetOutput();
	//Write(tissueContour,"tissueContour.nii");

	// get signed distance from tissue
	typedef itk::MorphologicalSignedDistanceTransformImageFilter<ByteImageType,FloatImageType> DistanceType;
	DistanceType::Pointer dFilter = DistanceType::New();
	dFilter->SetInput(tissue);
	dFilter->SetOutsideValue(255);
	dFilter->Update();
	Write(dFilter->GetOutput(),"distanceTissue.nii");

	// get zero crossing
	typedef itk::ZeroCrossingImageFilter<FloatImageType,ByteImageType> ZeroCrossingType;
	ZeroCrossingType::Pointer crossingFilter = ZeroCrossingType::New();
	crossingFilter->SetInput(dFilter->GetOutput());
	crossingFilter->Update();
	Write(crossingFilter->GetOutput(),"dzc.nii");

	// get edge potential map
	FloatImageType::Pointer ep = FloatImageType::New();
	ep->SetRegions(region);
	ep->SetSpacing(input->GetSpacing());
	ep->CopyInformation(input);
	ep->Allocate();
	ep->FillBuffer(0);
	FloatIteratorType epIt(ep,region);
	FloatIteratorType gIt(gradientMagnitude,region);

	for (gIt.GoToBegin(), epIt.GoToBegin(); !gIt.IsAtEnd(); ++gIt, ++epIt)
	{
		epIt.Set( 1 / ( 1 + gIt.Get() ) );
	}

	Write(ep,"edgePotential.nii");

	typedef itk::CurvesLevelSetImageFilter< FloatImageType, FloatImageType > CurvesFilterType;
	CurvesFilterType::Pointer curvesFilter = CurvesFilterType::New();
	curvesFilter->SetInput( dFilter->GetOutput() );
	curvesFilter->SetFeatureImage( ep );
	curvesFilter->SetPropagationScaling( 1 );
	curvesFilter->SetCurvatureScaling( 0.5 );
	curvesFilter->SetAdvectionScaling( 0.4 );
	curvesFilter->SetNumberOfIterations(50);
	curvesFilter->SetNumberOfThreads(2);
	curvesFilter->Update();

	Write(curvesFilter->GetOutput(),"curves.nii");

	crossingFilter = ZeroCrossingType::New();
	crossingFilter->SetInput(dFilter->GetOutput());
	crossingFilter->Update();
	Write(crossingFilter->GetOutput(),"lzc.nii");

}

void HeteroStoolRemoval(ImageType::Pointer &cOutput, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap)
{
	ImageType::RegionType region = cOutput->GetLargestPossibleRegion();

	// get air mask
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(cOutput->GetSpacing());
	air->Allocate();
	air->FillBuffer(0);

	IteratorType cOutputIt(cOutput,region);
	ByteIteratorType airIt(air,region);

	for (cOutputIt.GoToBegin(), airIt.GoToBegin(); !airIt.IsAtEnd(); ++cOutputIt, ++airIt)
	{
		if (cOutputIt.Get() < -600)
			airIt.Set(255);
	}

	Write(air,"air.nii");

	// remove bkg air
	typedef itk::BinaryShapeOpeningImageFilter< ByteImageType2D > BinaryShapeOpeningImageFilter2D;
	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, BinaryShapeOpeningImageFilter2D > SliceBySliceImageFilterBackgroundType;

	BinaryShapeOpeningImageFilter2D::Pointer bkgFilter2D = BinaryShapeOpeningImageFilter2D::New();
	bkgFilter2D->SetAttribute("SizeOnBorder");
	bkgFilter2D->SetBackgroundValue(0);
	bkgFilter2D->SetForegroundValue(255);
	bkgFilter2D->SetLambda(0);
	bkgFilter2D->SetReverseOrdering(true);
	
	SliceBySliceImageFilterBackgroundType::Pointer bkgRemover = SliceBySliceImageFilterBackgroundType::New();
	bkgRemover->SetInput( air );
	bkgRemover->SetFilter( bkgFilter2D );
	bkgRemover->Update();
	air = bkgRemover->GetOutput();

	Write(air,"air2.nii");

	// get distance map from air
	typedef itk::MorphologicalDistanceTransformImageFilter<ByteImageType,FloatImageType> DistanceType;
	DistanceType::Pointer distanceFilter = DistanceType::New();
	distanceFilter->SetInput(air);
	distanceFilter->SetOutsideValue(255);
	distanceFilter->Update();
	FloatImageType::Pointer airDist = distanceFilter->GetOutput();
	FloatIteratorType airDistIt(airDist,region);

	// set distance cutoff for full tissue
	float distCutoff = 3;

	for (airDistIt.GoToBegin(); !airDistIt.IsAtEnd(); ++airDistIt)
	{
		if (airDistIt.Get() > distCutoff)
			airDistIt.Set(distCutoff); 
	}	
	
	Write(airDist,"airDist.nii");

	// get sd
	FloatImageType::Pointer sd = StandardDeviation(cOutput,colon,2);
	FloatIteratorType sdIt(sd,region);
	Write(sd,"sd.nii");

	float sdCutoff = 350;

	// set full tissue if neighboring full tissue and low sd
	typedef itk::ImageDuplicator<FloatImageType> DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(airDist);
	duplicator->Update();
	FloatImageType::Pointer airDist2 = duplicator->GetOutput();
	FloatIteratorType airDist2It(airDist2,region);

	typedef itk::NeighborhoodIterator<FloatImageType> NeighborhoodIteratorType;
	FloatImageType::SizeType radius;
	radius.Fill(1);
	NeighborhoodIteratorType anIt(radius,airDist,region);

	for (anIt.GoToBegin(), sdIt.GoToBegin(), airDist2It.GoToBegin(); !anIt.IsAtEnd(); ++anIt, ++sdIt, ++airDist2It)
	{
		if (anIt.GetCenterPixel() < distCutoff && sdIt.Get() < sdCutoff)
		{
			for (unsigned int i=0; i<anIt.Size(); i++)
			{
				if (anIt.GetPixel(i) >= distCutoff)
				{
					airDist2It.Set(distCutoff);
					break;
				}
			}
		}
	}

	Write(airDist2,"airDist2.nii");

}

FloatImageType::Pointer StandardDeviation(ImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	FloatIteratorType outIt(out,region);
	ByteIteratorType maskIt(mask,region);
	
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(), maskIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt, ++maskIt)
	{
		// inside mask
		if (maskIt.Get() != 0)
		{
			float mean = 0;
			float std = 0;

			// get mean
			for (int i=0; i<nIt.Size(); i++)
			{
				mean += (float) nIt.GetPixel(i);
			}

			mean /= (float) nIt.Size();

			//std::cout << mean << std::endl;

			// get standard deviation
			for (int i=0; i<nIt.Size(); i++)
			{
				std += ( (float) nIt.GetPixel(i) - mean ) * ( (float) nIt.GetPixel(i) - mean );
			}

			std /= ( (float) nIt.Size()-1);
			std = sqrt(std);

			outIt.Set(std);
		}
	}

	return out;
}

FloatImageType::Pointer StandardDeviation(FloatImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius)
{
	// get region
	FloatImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	FloatImageType::Pointer out = FloatImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	FloatIteratorType outIt(out,region);
	ByteIteratorType maskIt(mask,region);
	
	typedef itk::NeighborhoodIterator<FloatImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(), maskIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt, ++maskIt)
	{
		// inside mask
		if (maskIt.Get() != 0)
		{
			float mean = 0;
			float std = 0;

			// get mean
			for (int i=0; i<nIt.Size(); i++)
			{
				mean += (float) nIt.GetPixel(i);
			}

			mean /= (float) nIt.Size();

			//std::cout << mean << std::endl;

			// get standard deviation
			for (int i=0; i<nIt.Size(); i++)
			{
				std += ( (float) nIt.GetPixel(i) - mean ) * ( (float) nIt.GetPixel(i) - mean );
			}

			std /= ( (float) nIt.Size()-1);
			std = sqrt(std);

			outIt.Set(std);
		}
	}

	return out;
}

ImageType::Pointer Subtraction(ImageType::Pointer &input, ImageType::Pointer &inputOriginal, ByteImageType::Pointer &colon, ArrayImageType::Pointer &partial, VoxelImageType::Pointer &vmap)
{
	// scale intensity by partial of tissue
	//ImageType::RegionType region = input->GetLargestPossibleRegion();

	//typedef itk::ImageDuplicator<ImageType> DuplicatorType;
	//DuplicatorType::Pointer duplicator = DuplicatorType::New();
	//duplicator->SetInputImage(input);
	//duplicator->Update();
	//ImageType::Pointer input2 = duplicator->GetOutput();

	//IteratorType inputIt(input,region);
	//IteratorType input2It(input2,region);
	//ByteIteratorType colonIt(colon,region);
	//ArrayIteratorType partialIt(partial,region);

	//for (inputIt.GoToBegin(), input2It.GoToBegin(), colonIt.GoToBegin(), partialIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++input2It, ++colonIt, ++partialIt)
	//{
	//	if (colonIt.Get() != 0)
	//	{
	//		input2It.Set( partialIt.Get()[1] * ( inputIt.Get() + 1000 ) - 1000 );	
	//	}
	//}

	//Write(input2,"carstonOutput.nii");

	//// plug cropped input back into original
	//DuplicatorType::Pointer duplicator = DuplicatorType::New();
	//duplicator->SetInputImage(input);
	//duplicator->Update();
	//ImageType::Pointer input2 = duplicator->GetOutput();

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// get pt binary mask
	ByteImageType::Pointer ptb = ByteImageType::New();
	ptb->SetRegions(region);
	ptb->CopyInformation(input);
	ptb->Allocate();
	ptb->FillBuffer(0);
	ByteIteratorType ptbIt(ptb,region);
	ArrayIteratorType partialIt(partial,region);

	for (ptbIt.GoToBegin(), partialIt.GoToBegin(); !ptbIt.IsAtEnd(); ++ptbIt, ++partialIt)
	{
		if (partialIt.Get()[1] > 0)
		{
			ptbIt.Set(255);
		}
	}

	Write(ptb,"partialTissueBinaryMask.nii");

	// keep only largest component
	typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType> KeeperType;
	KeeperType::Pointer keeper = KeeperType::New();
	keeper->SetInput(ptb);
	keeper->SetForegroundValue(255);
	keeper->SetBackgroundValue(0);
	keeper->SetNumberOfObjects(1);
	keeper->SetAttribute("Size");
	keeper->Update();
	ptb = keeper->GetOutput();
	Write(ptb,"partialTissueBinaryMaskCompd.nii");
	
	// update partial image
	ptbIt = ByteIteratorType(ptb,region);

	for (ptbIt.GoToBegin(), partialIt.GoToBegin(); !ptbIt.IsAtEnd(); ++ptbIt, ++partialIt)
	{
		if (ptbIt.Get() == 0 && partialIt.Get()[1] != 0)
		{
			ArrayType p = partialIt.Get();
			p[1] = 0;
			partialIt.Set(p);
		}
	}

	// scale input by partial of tissue
	IteratorType inputIt(input,region);
	ByteIteratorType colonIt(colon,region);

	for (inputIt.GoToBegin(), colonIt.GoToBegin(), partialIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt, ++partialIt)
	{
		if (colonIt.Get() != 0)
		{
			inputIt.Set( partialIt.Get()[1] * ( inputIt.Get() + 1000 ) - 1000 );	
		}
	}

	Write(input,"inputPtScaled.nii");

	// gaussian blur edges of air and tissue
	// get binary air mask
	ByteImageType::Pointer air = ByteImageType::New();
	air->SetRegions(region);
	air->SetSpacing(input->GetSpacing());
	air->Allocate();
	air->FillBuffer(0);
	ByteIteratorType airIt(air,region);

	for (airIt.GoToBegin(), inputIt.GoToBegin(); !airIt.IsAtEnd(); ++inputIt, ++airIt)
	{
		if (inputIt.Get() < -600)
		{
			airIt.Set(255);
		}
	}

	// remove bkg air and get edges on each 2d slice
	typedef itk::BinaryShapeOpeningImageFilter< ByteImageType2D > BinaryShapeOpeningImageFilter2D;
	typedef itk::BinaryContourImageFilter<ByteImageType2D,ByteImageType2D> BinaryEdgeFilterType2D;
	typedef itk::SliceBySliceImageFilter< ByteImageType, ByteImageType, BinaryShapeOpeningImageFilter2D, BinaryEdgeFilterType2D> SliceBySliceImageFilterType;

	BinaryShapeOpeningImageFilter2D::Pointer bkgFilter2D = BinaryShapeOpeningImageFilter2D::New();
	bkgFilter2D->SetAttribute("SizeOnBorder");
	bkgFilter2D->SetBackgroundValue(0);
	bkgFilter2D->SetForegroundValue(255);
	bkgFilter2D->SetLambda(0);
	bkgFilter2D->SetReverseOrdering(true);

	BinaryEdgeFilterType2D::Pointer edgeFinder2D = BinaryEdgeFilterType2D::New();
	edgeFinder2D->SetForegroundValue(255);
	edgeFinder2D->SetBackgroundValue(0);
	edgeFinder2D->SetInput(bkgFilter2D->GetOutput());
	
	SliceBySliceImageFilterType::Pointer slicer = SliceBySliceImageFilterType::New();
	slicer->SetInput( air );
	slicer->SetInputFilter( bkgFilter2D );
	slicer->SetOutputFilter( edgeFinder2D );
	slicer->Update();
	air = slicer->GetOutput();

	Write(air,"airEdges.nii");

	// dilate air edge
	Dilate(air,3);

	Write(air,"airEdgesDilated.nii");

	// smooth only near edge tissue interface
	typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> SmoothType;
	SmoothType::Pointer smoother = SmoothType::New();
	smoother->SetInput(input);
	
	/*ArrayType sigmaArray;
	sigmaArray[0] = 0.7;
	sigmaArray[1] = 0.7;
	sigmaArray[2] = 0;*/

	smoother->SetVariance(0.7*0.7);
	smoother->Update();
	ImageType::Pointer inputSmooth = smoother->GetOutput();

	Write(inputSmooth,"inputSmooth.nii");

	IteratorType inputSmoothIt(inputSmooth,region);

	airIt = ByteIteratorType(air,region);
	inputIt = IteratorType(input,region);

	for (inputIt.GoToBegin(), airIt.GoToBegin(), inputSmoothIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputSmoothIt, ++airIt)
	{
		if (airIt.Get() != 0)
		{
			inputIt.Set( inputSmoothIt.Get() );
		}
	}

	Write(input,"inputPostSmoothNearAir.nii");

	//// blur input within air edge mask
	//typedef itk::GaussianBlurImageFunction< ImageType > GFunctionType;
	//GFunctionType::Pointer gaussianFunction = GFunctionType::New();
	//gaussianFunction->SetInputImage( input );

	//GFunctionType::ErrorArrayType setError;
	//setError.Fill( 0.01 );
	//gaussianFunction->SetMaximumError( setError );

	//// blur only in 2D
	//ArrayType sigmaArray;
	//sigmaArray[0] = 0.7;
	//sigmaArray[1] = 0.7;
	//sigmaArray[2] = 0;

	//gaussianFunction->SetSigma( sigmaArray );

	//std::cout << "gaussian kernel: " << gaussianFunction->GetMaximumKernelWidth() << std::endl;
	//

	//for (airIt.GoToBegin(), inputIt.GoToBegin(); !airIt.IsAtEnd(); ++airIt, ++inputIt)
	//{
	//	if (airIt.Get() != 0)
	//	{
	//		inputIt.Set( (PixelType) gaussianFunction->EvaluateAtIndex(inputIt.GetIndex()) );
	//	}
	//}

	//Write(input,"inputBlur.nii");

	// plug input back into full original input
	IteratorType inputOriginalIt(inputOriginal,OLDREGION);

	for (inputIt.GoToBegin(), inputOriginalIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputOriginalIt, ++colonIt)
	{
		if (colonIt.Get() != 0)
		{
			inputOriginalIt.Set( inputIt.Get() );
		}
	}

	Write2(inputOriginal,"output.nii");

	return inputOriginal;
}

//void TextureAnalysis(ImageType::Pointer &input)
//{
//	// get gradient (no smoothing)
//	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeImageFilterType;
//	GradientMagnitudeImageFilterType::Pointer gmFilter = GradientMagnitudeImageFilterType::New();
//	gmFilter->SetInput(input);
//
//	// threshold gradient
//	typedef itk::BinaryThresholdImageFilter<FloatImageType,ByteImageType> BinaryThresholdImageFilterType;
//	BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
//	thresholdFilter->SetInsideValue(1);
//	thresholdFilter->SetOutsideValue(0);
//	thresholdFilter->SetLowerThreshold(400);
//	thresholdFilter->SetInput(gmFilter->GetOutput());
//	thresholdFilter->Update();
//
//	// write mask
//	typedef itk::MultiplyByConstantImageFilter<ByteImageType,unsigned char,ByteImageType> MultiplyFilterType;
//	MultiplyFilterType::Pointer multiplier = MultiplyFilterType::New();
//	multiplier->SetInput(thresholdFilter->GetOutput());
//	multiplier->SetConstant(255);
//	multiplier->Update();
//	Write(multiplier->GetOutput(),"mask.nii");
//
//	// get image minimum and maximum
//	typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
//	MinimumMaximumImageCalculatorType::Pointer imageCalc = MinimumMaximumImageCalculatorType::New();
//	imageCalc->SetImage(input);
//	imageCalc->ComputeMinimum();
//	imageCalc->ComputeMaximum();
//
//	typedef otb::ScalarImageToTexturesFilter2<ImageType,FloatImageType,ByteImageType> TextureFilterType;
//	TextureFilterType::Pointer textureFilter = TextureFilterType::New();
//	textureFilter->SetInput(input);
//	textureFilter->SetMaskImage(thresholdFilter->GetOutput());
//	textureFilter->SetInputImageMinimum(imageCalc->GetMinimum());
//	textureFilter->SetInputImageMaximum(imageCalc->GetMaximum());
//	
//	ImageType::OffsetType offset;
//	offset[0] = 1;
//	offset[1] = 0;
//	offset[2] = 0;
//
//	textureFilter->SetOffset(offset);
//
//	ImageType::SizeType radius;
//	radius[0] = 1;
//	radius[1] = 1;
//	radius[2] = 0;
//
//	textureFilter->SetRadius(radius);
//	textureFilter->Update();
//
//	for (int i=0; i<8; i++)
//	{
//		std::stringstream ss;
//		ss << "texture" << i << ".nii";
//		Write(textureFilter->GetOutput(i),ss.str());
//	}
//
//}

void DirectionalGradient(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap)
{
	// create binary tissue mask
	ByteImageType::Pointer mask = ByteImageType::New();
	mask->SetSpacing(input->GetSpacing());
	mask->SetRegions(REGION);
	mask->Allocate();
	mask->FillBuffer(0);
	
	ByteIteratorType maskIt(mask,REGION);
	VoxelIteratorType vmapIt(vmap,REGION);

	for (maskIt.GoToBegin(), vmapIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++vmapIt)
	{
		if (vmapIt.Get() == Air)
			maskIt.Set(255);
	}

	Write(mask,"airMask.nii");

	// mask dg adjacent to air only

	typedef itk::ImageDuplicator<ByteImageType> ImageDuplicatorType;
	ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
	duplicator->SetInputImage(mask);
	duplicator->Update();
	ByteImageType::Pointer maskDilated = duplicator->GetOutput();

	Dilate(maskDilated,3);
	Write(maskDilated,"airMaskDilated.nii");

	typedef itk::SubtractImageFilter<ByteImageType,ByteImageType,ByteImageType> SubtractImageFilterType;
	SubtractImageFilterType::Pointer subtracter = SubtractImageFilterType::New();
	subtracter->SetInput1(maskDilated);
	subtracter->SetInput2(mask);
	subtracter->Update();
	ByteImageType::Pointer airBorder = subtracter->GetOutput();

	// only see dg in unclassified air border
	ByteIteratorType airBorderIt(airBorder,REGION);

	for (airBorderIt.GoToBegin(), vmapIt.GoToBegin(); !airBorderIt.IsAtEnd(); ++airBorderIt, ++vmapIt)
	{
		if (vmapIt.Get() != Unclassified)
			airBorderIt.Set(0);
	}

	Write(airBorder,"airBorder.nii");

	// directional gradient
	typedef itk::DirectionalGradientImageFilter2<ImageType,ByteImageType,ImageType> DirectionalGradientImageFilterType;
	DirectionalGradientImageFilterType::Pointer dgFilter = DirectionalGradientImageFilterType::New();
	dgFilter->SetInput(input);
	//dgFilter->SetSigma(input->GetSpacing()[0]);
	dgFilter->SetScale(-1);
	dgFilter->SetMaskImage(mask);
	dgFilter->SetOutsideValue(255);
	dgFilter->Update();
	ImageType::Pointer dg = dgFilter->GetOutput();
	IteratorType dgIt(dg,REGION);

	// set -ve to 0
	for (dgIt.GoToBegin(); !dgIt.IsAtEnd(); ++dgIt)
	{
		if (dgIt.Get() < 0)
			dgIt.Set(0);
	}

	Write(dg,"dg.nii");
	
	typedef itk::MaskImageFilter<ImageType,ByteImageType,ImageType> MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1(dg);
	masker->SetInput2(airBorder);
	masker->Update();
	Write(masker->GetOutput(),"dgMasked.nii");
	
}

void Dilate(ByteImageType::Pointer &img, unsigned int radius)
{
	StructuringElementType se;
	
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;

	se.SetRadius( rad );
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<ByteImageType, ByteImageType, StructuringElementType> BinaryDilateImageFilterType;
	BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput( img );
	dilater->SetKernel( se );
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();

	img = dilater->GetOutput();
}


/*********************************************************
- Load image
- Shift input so that air is 0 HU
- Segment colon
- Mask input with colon
- Crop input and colon in XY plane
*********************************************************/
void Setup(std::string dataset, ImageType::Pointer  &inputOriginal, ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude)
{
	//----------------------------------------------
	// Load image
	//----------------------------------------------
	std::vector<std::string> datasetArray = explode( "\\", dataset );
	std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	if (truncateOn)
	{
		std::cout << "Truncating data" << std::endl;
		std::cout << "Slices: " << truncateArray[0] << " to " << truncateArray[1] << std::endl;
		
		std::stringstream ss;
		ss << dsname << "_" << truncateArray[0] << "_" << truncateArray[1];
		dsname = ss.str();
	}

	// Set writer prefix
	note = dsname;

	// Load dicom files
	if (!truncateOn)
	{
		inputOriginal = ReadDicom < ImageType > ( dataset );
	} else {
		inputOriginal = ReadDicom < ImageType > ( dataset, truncateArray[0], truncateArray[1] );
	}

	// Set global region
	REGION = inputOriginal->GetLargestPossibleRegion();	

	//Write2(inputOriginal,"inputOriginal.nii");

	// Get image minimum and set as global background
	typedef itk::MinimumMaximumImageCalculator<ImageType> MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minCalc = MinimumMaximumImageCalculatorType::New();
	minCalc->SetImage(inputOriginal);
	minCalc->ComputeMinimum();
	BACKGROUND = minCalc->GetMinimum();

	//----------------------------------------------
	// Segment colon
	//----------------------------------------------
	typedef itk::ColonSegmentationFilter< ImageType, ByteImageType > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( inputOriginal );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	//colonSegmenter->SetPrintImages(true);

	if ( truncateOn )
		colonSegmenter->SetRemoveBoneLung( false );

	colonSegmenter->Update();
	colon = colonSegmenter->GetOutput();

	//----------------------------------------------
	// Mask input with colon
	//----------------------------------------------
	typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( inputOriginal );
	masker->SetInput2( colon );
	masker->SetOutsideValue( BACKGROUND );
	masker->Update();
	input = masker->GetOutput();

	//----------------------------------------------
	// Calculate gradient magnitude
	//----------------------------------------------
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientMagnitudeFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientMagnitudeFilter->SetInput( inputOriginal );
	gradientMagnitudeFilter->SetSigma( inputOriginal->GetSpacing()[0] );
	gradientMagnitudeFilter->Update();
	gradientMagnitude = gradientMagnitudeFilter->GetOutput();

	/*typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType> GradientMagnitudeFilterType;
	GradientMagnitudeFilterType::Pointer gmFilter = GradientMagnitudeFilterType::New();
	gmFilter->SetInput( inputOriginal );
	gmFilter->Update();
	gradientMagnitude = gmFilter->GetOutput();*/

	//----------------------------------------------
	// Crop images in XY plane
	//----------------------------------------------

	ImageType::SizeType size = REGION.GetSize();

	long minX=size[0],minY=size[1],maxX=0,maxY=0;
	short paddingXY = 5;
	
	IteratorType inputIt(input,REGION);
	ByteIteratorType colonIt(colon,REGION);

	for (inputIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++colonIt)
	{
		if ( colonIt.Get() != 0 )
		{
			ImageType::IndexType idx = inputIt.GetIndex();

			if (idx[0] < minX)
				minX = idx[0];
			if (idx[0] > maxX)
				maxX = idx[0];
			if (idx[1] < minY)
				minY = idx[1];
			if (idx[1] > maxY)
				maxY = idx[1];
		}
	}

	ImageType::IndexType edx;
	
	edx[0] = (minX-paddingXY) > 0 ? minX-paddingXY : 0;
	edx[1] = (minY-paddingXY) > 0 ? minY-paddingXY : 0;
	edx[2] = 0;

	ImageType::SizeType esize;
	esize[0] = maxX-minX+2*paddingXY+1 < size[0] ? maxX-minX+2*paddingXY+1 : size[0];
	esize[1] = maxY-minY+2*paddingXY+1 < size[1] ? maxY-minY+2*paddingXY+1 : size[1];
	esize[2] = size[2];

	ImageType::RegionType extractRegion;
	extractRegion.SetIndex( edx );
	extractRegion.SetSize( esize );

	OLDREGION = extractRegion;

	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( input );
	cropper->SetRegionOfInterest( extractRegion );
	cropper->Update();
	input = cropper->GetOutput();

	typedef itk::RegionOfInterestImageFilter<FloatImageType,FloatImageType> RegionOfInterestImageFilterFloatType;
	RegionOfInterestImageFilterFloatType::Pointer cropperFloat = RegionOfInterestImageFilterFloatType::New();
	cropperFloat->SetInput( gradientMagnitude );
	cropperFloat->SetRegionOfInterest( extractRegion );
	cropperFloat->Update();
	gradientMagnitude = cropperFloat->GetOutput();

	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();

	// Set cropped region globally
	REGION = colon->GetLargestPossibleRegion();

	Write(input,"input.nii");
	//Write(colon,"colon.nii");
	Write(gradientMagnitude,"gradientMagnitudeSmoothed.nii");
}

/*********************************************************
- Allocate voxel map
- Compute otsu threshold separating tissue and stool
- Apply initial threshold rules
*********************************************************/
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon) 
{
	//----------------------------------------------
	// Allocate voxel map
	//----------------------------------------------
	vmap->SetRegions( input->GetLargestPossibleRegion() );
	vmap->SetSpacing( input->GetSpacing() );
	vmap->SetDirection( input->GetDirection() );
	vmap->CopyInformation( input );
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);

	//----------------------------------------------
	// Compute otsu threshold separating tissue and stool
	//----------------------------------------------

	//typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	//OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	//otsuCalculator->SetImage( input );
	//otsuCalculator->SetMinMax(true);
	//otsuCalculator->SetHistogramMin(-300);
	//otsuCalculator->SetHistogramMax(1500);
	////otsuCalculator->SetPrintHistogram(note+"_intensity.csv");
	//otsuCalculator->Compute();

	//PixelType tissueStoolThreshold = otsuCalculator->GetThreshold();
	//std::cout << "Tissue Stool Otsu Threshold: " << tissueStoolThreshold << std::endl;

	PixelType tissueStoolThreshold = 180;

	//----------------------------------------------
	// Apply initial threshold rules
	//----------------------------------------------
	ApplyThresholdRules( input, input->GetLargestPossibleRegion(), gradientMagnitude, vmap, colon, tissueStoolThreshold );

	Write(vmap,"vmap.nii");

	// output means of air tissue stool classes
	float mean[3] = {0,0,0};
	int count[3] = {0,0,0};
	
	VoxelIteratorType vmapIt(vmap,vmap->GetLargestPossibleRegion());
	IteratorType inputIt(input,input->GetLargestPossibleRegion());
	
	for (vmapIt.GoToBegin(), inputIt.GoToBegin(); !vmapIt.IsAtEnd(); ++vmapIt, ++inputIt)
	{
		VoxelType v = vmapIt.Get();
		PixelType I = (float) inputIt.Get();

		if (v == Air)
		{
			mean[0] += I;
			count[0]++;
		} else if (v == Tissue) {
			mean[1] += I;
			count[1]++;
		} else if (v == Stool) {
			mean[2] += I;
			count[2]++;
		}
	}

	std::cout << "Mean" << std::endl;
	std::cout << "Air: " << mean[0]/count[0] << std::endl;
	std::cout << "Tissue: " << mean[1]/count[1] << std::endl;
	std::cout << "Stool: " << mean[2]/count[2] << std::endl;

	return tissueStoolThreshold;
}

void ApplyThresholdRules( ImageType::Pointer &input, ImageType::RegionType localRegion, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissueStoolThreshold )
{
	IteratorType inputIt(input,localRegion);
	FloatIteratorType gradientMagnitudeIt(gradientMagnitude,localRegion);
	VoxelIteratorType vmapIt(vmap,localRegion);
	ByteIteratorType colonIt(colon,localRegion);

	for ( inputIt.GoToBegin(), gradientMagnitudeIt.GoToBegin(), vmapIt.GoToBegin(), colonIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++gradientMagnitudeIt, ++vmapIt, ++colonIt)
	{
		if (colonIt.Get() == 255)
		{
			VoxelType voxel;

			short I = inputIt.Get();
			float G = gradientMagnitudeIt.Get();

			if ( ( I >= tissueStoolThreshold && G < 0.8*I ) || I > 1000 )
			{
				voxel = Stool;
			} else if ( I <= -700 ) {
				voxel = Air;
			} else if ( I < tissueStoolThreshold && I > -300 && G <= 300 ) {
				voxel = Tissue;
			} else {
				voxel = Unclassified;
			}

			vmapIt.Set( voxel );
		}
	}
}

ImageType::Pointer Range(ImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// set radius
	ImageType::SizeType rad;
	rad.Fill(radius);

	// allocate output image
	ImageType::Pointer out = ImageType::New();
	out->SetSpacing(input->GetSpacing());
	out->SetRegions(region);
	out->Allocate();
	out->FillBuffer(0);

	// iterate image
	IteratorType outIt(out,region);
	ByteIteratorType maskIt(mask,region);
	
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nIt(rad,input,region);

	for (nIt.GoToBegin(), outIt.GoToBegin(), maskIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt, ++outIt, ++maskIt)
	{
		// inside mask
		if (maskIt.Get() != 0)
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
	}

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

FloatImageType::Pointer RescaledRange(ImageType::Pointer &input, unsigned int radius)
{

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Make map image
	ImageType::Pointer map = ImageType::New();
	
	ImageType::RegionType mapRegion;
	
	ImageType::IndexType mapIndex;
	ImageType::SizeType mapSize;

	for (int i=0; i<3; i++)
	{
		mapIndex[i] = 0;
		mapSize[i] = 2*radius+1;
	}

	mapRegion.SetIndex(mapIndex);
	mapRegion.SetSize(mapSize);

	map->SetRegions(mapRegion);

	map->Allocate();
	map->FillBuffer(0);

	IteratorType mapIt(map,mapRegion);

	// Setup neighborhood and fill map with distance squared
	typedef itk::Neighborhood<PixelType,3> NeighborhoodType;
	NeighborhoodType hood;
	hood.SetRadius(radius);

	typedef ImageType::OffsetType OffsetType;
	
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

	Write(map,"map.nii");

	// Convert map to label map to get number of distance classes
	typedef itk::LabelImageToLabelMapFilter<ImageType> LabelImageToLabelMapFilterType;
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
	
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,ImageType> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToImageFilter->SetInput(relabeler->GetOutput());
	labelMapToImageFilter->Update();

	map = labelMapToImageFilter->GetOutput();
	mapIt = IteratorType(map,mapRegion);

	labelMap.~SmartPointer();

	Write(map,"map2.nii");

	// Allocate vectors to store min and max for each class
	std::vector<PixelType> minVector;
	std::vector<PixelType> maxVector;

	minVector.resize(numClasses);
	maxVector.resize(numClasses);

	for (int i=0; i<numClasses; i++)
	{
		minVector[i] = ( itk::NumericTraits<PixelType>::max() );
		maxVector[i] = ( itk::NumericTraits<PixelType>::NonpositiveMin() );
	}
	
	/*struct srange {
		PixelType min;
		PixelType max;
	};
	
	std::vector<srange> rangeVector;
	rangeVector.resize(numClasses);

	for (int i=0; i<numClasses; i++)
	{
		struct srange r;
		r.min = itk::NumericTraits<PixelType>::max();
		r.max = itk::NumericTraits<PixelType>::NonpositiveMin();
		rangeVector.push_back(r);
	}*/

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
	ImageType::SizeType rad;
	rad.Fill(radius);
	
	typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(rad,input,region);

	for (nit.GoToBegin(), oit.GoToBegin(); !nit.IsAtEnd(); ++nit, ++oit)
	{

		// Reset min and max vectors
		for (int i=0; i<numClasses; i++)
		{
			minVector[i] = itk::NumericTraits<PixelType>::max();
			maxVector[i] = itk::NumericTraits<PixelType>::NonpositiveMin();

			rv[i] = 0;
		}

		// Iterate neighborhood and map image together
		int j=0;

		for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
		{
			if ( mapIt.Get() > 0 )
			{
				PixelType val = nit.GetPixel(j);
				PixelType m = mapIt.Get() - 1;

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