#include "RemoveStool4.h"
#include "io2.cxx"
//#include "io.cxx"
#include "utility.cxx"
#include "scattercorrection.cxx"
#include "QR.cxx"
#include "EM.cxx"
#include "HessianFunctions2.cxx"


int main(int argc, char * argv[])
{
	// Start clock
	init = clock();

	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " InputDicomDirectory OutputDicomDirectory";
		system("pause");
		return EXIT_FAILURE;
	}
	
	std::string inDir = argv[1];
	std::string outDir = argv[2];

	std::cout << "Input Directory: " << inDir << std::endl;
	
	////////////////////////////////////////////
	// Load input
	////////////////////////////////////////////
	typedef itk::ImageSeriesReader< ImageType >     ReaderType;

	typedef itk::GDCMImageIO                        ImageIOType;
	typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

	namesGenerator->SetInputDirectory( inDir );

	const ReaderType::FileNamesContainer & filenames = 
						namesGenerator->GetInputFileNames();

	unsigned int numberOfFilenames =  filenames.size();
	std::cout << numberOfFilenames << std::endl; 

	// copy names
	std::vector<std::string> names;
	names.resize( filenames.size() );

	for (int i=0; i<names.size(); i++)
	{
	names[i] = filenames[i];
	}

	// reverse names
	reverse(names.begin(),names.end());

	/*for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
	{
	std::cout << "filename # " << fni << " = ";
	std::cout << names[fni] << std::endl;
	}*/

	ReaderType::Pointer reader = ReaderType::New();

	reader->SetImageIO( gdcmIO );
	reader->SetFileNames( names );

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << excp << std::endl;
		system("pause");
	}

	ImageType::Pointer output = reader->GetOutput();
	Write(output,"inputDicom.nii");

	////////////////////////////////////////////
	// Run the algorithm
	////////////////////////////////////////////
	
	output = RunAlgorithm( output );

	Write2(output,"output.nii");
	
	////////////////////////////////////////////
	// Write to DICOM
	////////////////////////////////////////////

	// Setup output folder
	// Make output directory
	std::vector<std::string> datasetArray = explode( "\\", inDir );
	std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;
	
	outDir += "\\" + dsname + "\\dcm";
	
	std::cout << "Output Directory: " << outDir << std::endl;

	// Change input path names to output path names
	for (int i=0; i<names.size(); i++)
	{
		std::vector<std::string> namesArray = explode("/",names[i]);
		names[i] = outDir + "\\" + namesArray[ namesArray.size() - 1 ];
	}


	typedef itk::ImageSeriesWriter< 
                             ImageType, ImageType2D >  SeriesWriterType;

	const char * outputDirectory = outDir.c_str();

	itksys::SystemTools::MakeDirectory( outputDirectory );

	 SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

	 seriesWriter->SetInput( output );
	 
	 seriesWriter->SetImageIO( gdcmIO );
	
	 namesGenerator->SetOutputDirectory( outputDirectory );

	 seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
	 
	 seriesWriter->SetMetaDataDictionaryArray( 
	                       reader->GetMetaDataDictionaryArray() );

	 try
	   {
	   seriesWriter->Update();
	   }
	 catch( itk::ExceptionObject & excp )
	   {
		   std::cerr << "Exception thrown while writing the series " << std::endl;
		   std::cerr << excp << std::endl;
		   system("pause");
			//return EXIT_FAILURE;
	   }

	return 0;

}


ImageType::Pointer RunAlgorithm( ImageType::Pointer &inputOriginal )
{

	//ImageType::Pointer				inputOriginal			= ImageType::New();
	
	ImageType::Pointer				input					= ImageType::New();
	ByteImageType::Pointer			colon					= ByteImageType::New();
	FloatImageType::Pointer			gradientMagnitude		= FloatImageType::New();
	VoxelImageType::Pointer			vmap					= VoxelImageType::New();

	// Load images and segment colon
	Setup(inputOriginal,input,colon);

	// Smooth input
	/*typedef itk::FastBilateralImageFilter<ImageType,ImageType> BilateralFilterType;
	BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();
	bilateralFilter->SetInput(input);
	bilateralFilter->SetRangeSigma(200);
	bilateralFilter->SetDomainSigma(5);
	bilateralFilter->Update();
	input = bilateralFilter->GetOutput();
	Write(input,"inputSmoothed.nii");
	*/

	// Get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	//// Scatter correction
	//ImageType::Pointer inputScatter = ScatterCorrection(inputOriginal,colon);
	//Write(inputScatter,"inputScatter.nii");

	//// Make changes only in tagged regions (HU > 200)
	//IteratorType is(inputScatter,region);
	//IteratorType it(input,region);
	//ByteIteratorType cit(colon,region);

	//for (it.GoToBegin(), is.GoToBegin(), cit.GoToBegin(); !it.IsAtEnd(); ++it, ++is, ++cit)
	//{
	//	if (cit.Get() != 0)
	//	{
	//		if (it.Get() > 200)
	//		{
	//			it.Set( is.Get() );
	//		}
	//	}
	//}

	//inputScatter.~SmartPointer();

	//Write(input,"inputAfterScatter.nii");

	LocalThreshold(input,colon,gradientMagnitude,vmap);

	//// Single material classification, i.e. conservative intensity/gradient thresholding
	//typedef itk::RescaleIntensityImageFilter<ImageType,ByteImageType> RescalerType;
	//RescalerType::Pointer rescaler = RescalerType::New();
	//rescaler->SetInput(input);
	//rescaler->SetOutputMaximum(255);
	//rescaler->SetOutputMinimum(0);
	//rescaler->Update();
	//ByteImageType::Pointer inputByte = rescaler->GetOutput();
	//Write(inputByte,"inputByte.nii");
	//
	//
	//std::vector<FloatImageType::Pointer> rrVector = RescaledRange( inputByte ,colon,3);
	//Write(rrVector[0],"rrSlope.nii");

	//ImageType::Pointer mip = LocalMIP(input,colon,1);
	//Write(mip,"mip.nii");
	//
	//LocalThreshold(mip,colon,gradientMagnitude,vmap);*/

	//ByteImageType::Pointer unclassifiedMask = BinaryThreshold(vmap,Unclassified);
	//unclassifiedMask = Mask(unclassifiedMask,colon);
	//unclassifiedMask = Mask(unclassifiedMask, BinaryThreshold(gradientMagnitude,0,300));

	//FloatImageType::Pointer inputFloat = Cast <ByteImageType,FloatImageType> (inputByte);

	//typedef otb::ScalarImageToHaralicksCorrelationMaskFilter<FloatImageType,FloatImageType,ByteImageType> TexturesFilterType;
	//TexturesFilterType::Pointer textureFilter = TexturesFilterType::New();
	//textureFilter->SetInput( inputFloat );
	//textureFilter->SetMaskImage( unclassifiedMask );
	//
	//ImageType::SizeType radiusSize;
	//radiusSize.Fill(0);
	//radiusSize[0] = 2;
	//radiusSize[1] = 2;
	//
	//textureFilter->SetRadius(radiusSize);

	//ImageType::OffsetType offset;
	//offset[0] = 1;
	//offset[1] = 0;
	//offset[2] = 0;

	//textureFilter->SetOffset(offset);
	//textureFilter->SetInputImageMinimum(0);
	//textureFilter->SetInputImageMaximum(255);
	//textureFilter->SetNumberOfBinsPerAxis( 32 );
	//textureFilter->Update();

	//std::stringstream ss;
	//Write(textureFilter->GetOutput(),"haralick.nii");


	//PixelType tst = SingleMaterialClassification(input, gradientMagnitude, vmap, colon);

	//// Get black top hat
	//typedef itk::BlackTopHatImageFilter<ImageType,ImageType,StructuringElementType> BlackTopHatType;
	//BlackTopHatType::Pointer bthF = BlackTopHatType::New();
	//bthF->SetInput( input );

	//float radArray[3] = {1,2,3};

	//for (int i=0; i<3; i++)
	//{
	//	StructuringElementType se;
	//	
	//	ByteImageType::SizeType rad;
	//	rad.Fill(radArray[i]);

	//	se.SetRadius( rad );
	//	se.CreateStructuringElement();

	//	bthF->SetKernel(se);
	//	bthF->Update();
	//	ImageType::Pointer bth = bthF->GetOutput();

	//	bth = Mask(bth,BinaryThreshold(vmap,Stool));

	//	std::stringstream ss;
	//	ss << "bth" << radArray[i] << ".nii";
	//	Write(bth,ss.str());
	//}

	// Determine boundary types
	ArrayImageType::Pointer partial = QuadraticRegression(input,colon,vmap,gradientMagnitude,200);

	gradientMagnitude.~SmartPointer();
	
	//FixATT(input,partial,vmap,colon,smax);

	//smax.~SmartPointer();

	FloatImageType::Pointer hessian = SatoResponse4(input,colon,vmap,partial,0.25,1);
	Write(hessian,"hessian.nii");	

	// Expectation Maximization
	EM(partial,colon,input);

	// Perform subtraction
	ImageType::Pointer carstonOutput = Subtraction(input,inputOriginal,colon,partial,vmap);

	return carstonOutput;

	/*Write(carstonOutput,"carstonOutput.nii");
	Write(inputOriginal,"inputOriginal.nii");*/

	// Write dcm output
	
	//// Get dir names
	//std::string inDir = argv[1];
	//std::string outDir = argv[2];

	//// Make output directory
	//std::vector<std::string> datasetArray = explode( "\\", inDir );
	//std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	//std::cout << "Dataset: " << dsname << std::endl;
	//outDir += "\\" + dsname + "\\dcm";
	//std::cout << "Output Directory: " << outDir << std::endl;

	//const char * outputDirectory = outDir.c_str();
	//
	//itksys::SystemTools::MakeDirectory( outputDirectory );

	//// Create series writer
	//typedef itk::ImageSeriesWriter<ImageType,ImageType2D> SeriesWriterType;
	//SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	//
	//itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
	//seriesWriter->SetImageIO( dicomIO );

	//// Change names to output folder names
	//for (int i=0; i<names.size(); i++)
	//{
	//	std::vector<std::string> namesArray = explode("/",names[i]);
	//	names[i] = outDir + "\\" + namesArray[ namesArray.size() - 1 ];
	//}

	//seriesWriter->SetFileNames( names );

	//std::cout << "Number of filenames: " << names.size() << std::endl;
	//std::cout << "Number of dictionary entries: " << metaDataDictionaryArray.size() << std::endl;

	//for(int i=0; i<metaDataDictionaryArray.size(); i++)
	//{
	//	metaDataDictionaryArray[i] = &(dicomIO->GetMetaDataDictionary());
	//}*/

	//seriesWriter->SetMetaDataDictionaryArray( &metaDataDictionaryArray );//&metaDataDictionaryArray );
	//seriesWriter->SetInput( inputOriginal );

	//try
	//{
	//	seriesWriter->Update();
	//}
	//catch( itk::ExceptionObject & excp )
	//{
	//	std::cerr << "Exception thrown while writing the series " << std::endl;
	//	std::cerr << excp << std::endl;
	//	system("pause");
	//	return EXIT_FAILURE;
	//}

	//// End clock
	//std::cout << (double) (clock() - init) / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;

	//system("pause");
	//return 0;
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

void LocalThreshold(ImageType::Pointer input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap)
{	
	// Compute gradient magnitude
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientMagnitudeFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientMagnitudeFilter->SetInput( input );
	gradientMagnitudeFilter->SetSigma( input->GetSpacing()[0] );
	gradientMagnitudeFilter->Update();
	gradientMagnitude = gradientMagnitudeFilter->GetOutput();
	Write(gradientMagnitude,"gradientMagnitude.nii");
	gradientMagnitudeFilter.~SmartPointer();

	// Smooth input
	input = Median(input,1);
	Write(input,"inputMedian.nii");

	// Mask input with colon
	input = Mask(input,colon);

	// Get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	// Initialize vmap
	vmap = VoxelImageType::New();
	vmap->SetRegions(region);
	vmap->Allocate();
	vmap->FillBuffer(Unclassified);

	// Set initial threshold
	PixelType initialThreshold = 200;

	// Find tagged
	typedef itk::BinaryThresholdImageFilter<ImageType,ByteImageType> BinarizerType;
	BinarizerType::Pointer binarizer = BinarizerType::New();
	binarizer->SetInput(input);
	binarizer->SetOutsideValue(0);
	binarizer->SetInsideValue(255);
	binarizer->SetLowerThreshold(initialThreshold);
	binarizer->Update();
	ByteImageType::Pointer tag = binarizer->GetOutput();
	Write(tag,"allTagged.nii");

	// Open to separate thin stool connections
	tag = BinaryOpen(tag,3);
	Write(tag,"tagOpened.nii");
	
	// Remove small components
	typedef itk::BinaryShapeOpeningImageFilter<ByteImageType> BinaryOpeningType;
	BinaryOpeningType::Pointer opener = BinaryOpeningType::New();
	opener->SetInput(tag);
	opener->SetAttribute("Size");
	opener->SetForegroundValue(255);
	opener->SetBackgroundValue(0);
	opener->SetLambda(5);
	opener->Update();
	Write(opener->GetOutput(),"allTaggedRemovedSmall.nii");

	tag = opener->GetOutput();

	// Setup otsu
	typedef itk::OtsuThresholdImageCalculatorModified< ImageType > OtsuThresholdImageCalculatorModifiedType;
	OtsuThresholdImageCalculatorModifiedType::Pointer otsuCalculator = OtsuThresholdImageCalculatorModifiedType::New();
	otsuCalculator->SetImage(input);
	otsuCalculator->SetMinMax(true);
	otsuCalculator->SetHistogramMin(-300);
	otsuCalculator->SetHistogramMax(1600); //1200
	otsuCalculator->SetNumberOfHistogramBins(128);
	otsuCalculator->Compute();
	PixelType otsuGlobal = otsuCalculator->GetThreshold();

	std::cout << "global otsu: " << otsuGlobal << std::endl;

	// Setup threshold image
	ImageType::Pointer thresholdImage = ImageType::New();
	thresholdImage->SetRegions(region);
	thresholdImage->CopyInformation(input);
	thresholdImage->Allocate();

	// Set initial threshold
	thresholdImage->FillBuffer(initialThreshold);

	// Get label attributes
	typedef itk::BinaryImageToShapeLabelMapFilter<ByteImageType> BinaryImageToShapeLabelMapFilterType;
	typedef BinaryImageToShapeLabelMapFilterType::OutputImageType LabelMapType;
	typedef LabelMapType::LabelObjectType LabelObjectType;

	BinaryImageToShapeLabelMapFilterType::Pointer labeler = BinaryImageToShapeLabelMapFilterType::New();
	labeler->SetInput(tag);
	labeler->SetInputForegroundValue(255);
	labeler->SetOutputBackgroundValue(0);
	labeler->Update();
	LabelMapType::Pointer labelMap = labeler->GetOutput();

	// Convert label map to label image for viewing
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,LabelImageType> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer converter = LabelMapToLabelImageFilterType::New();
	converter->SetInput(labelMap);
	converter->Update();
	Write(converter->GetOutput(),"labelImage.nii");

	unsigned int padding = 5;

	for( unsigned int label=1; label<=labelMap->GetNumberOfLabelObjects(); label++ )
	{
		const LabelObjectType * labelObject = labelMap->GetLabelObject( label );
		
		// Get and label region
		ByteImageType::RegionType labelRegion = labelObject->GetRegion();
		ByteImageType::IndexType labelIndex = labelRegion.GetIndex();
		ByteImageType::SizeType labelSize = labelRegion.GetSize();

		// Expand region in XY plane
		
		for ( unsigned int j=0; j<2; j++ )
		{
			labelIndex[j] = labelIndex[j] - padding > 0 ? labelIndex[j] - padding : 0;

			labelSize[j] += 2*padding;

			int over = labelIndex[j] + labelSize[j] > region.GetSize()[j];

			if (over > 0)
				labelSize[j] -= over;
		}

		labelRegion.SetIndex(labelIndex);
		labelRegion.SetSize(labelSize);

		// Get otsu for sub region
		otsuCalculator->SetRegion(labelRegion);

		//std::stringstream ss;
		//ss << label << "_intensity.csv";
		//otsuCalculator->SetPrintHistogram(ss.str());

		otsuCalculator->Compute();

		PixelType ot = otsuCalculator->GetThreshold();

		std::cout << "Label: " << label << "\tOtsu: " << ot << std::endl;

		// Store largest threshold in image
		IteratorType it(thresholdImage,labelRegion);

		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			/*if (ot > it.Get())
			{
				it.Set(ot);
			}*/

			if (it.Get() == initialThreshold && ot > initialThreshold)
			{
				it.Set( ot );
			} else if (it.Get() > initialThreshold && ot < it.Get())
			{
				it.Set( ot );
			}
		}

	}

	Write(thresholdImage,"thresholdImage.nii");

	// Get high intensity mask
	//ByteImageType::Pointer highIMask = 

	//// Get texture 
	//typedef itk::RescaleIntensityImageFilter<ImageType,FloatImageType> RescalerType;
	//RescalerType::Pointer rescaler = RescalerType::New();
	//rescaler->SetInput(input);
	//rescaler->SetOutputMaximum(255);
	//rescaler->SetOutputMinimum(0);
	//
	//typedef otb::ScalarImageToHaralicksCorrelationMaskFilter<FloatImageType,FloatImageType,ByteImageType> TexturesFilterType;
	//TexturesFilterType::Pointer textureFilter = TexturesFilterType::New();
	//textureFilter->SetInput( rescaler->GetOutput() );
	//textureFilter->SetMaskImage(  );
	//
	//ImageType::SizeType radiusSize;
	//radiusSize.Fill(0);
	//radiusSize[0] = 1;
	//radiusSize[1] = 1;
	//
	//textureFilter->SetRadius(radiusSize);

	//ImageType::OffsetType offset;
	//offset[0] = 1;
	//offset[1] = 0;
	//offset[2] = 0;

	//textureFilter->SetOffset(offset);
	//textureFilter->SetInputImageMinimum(0);
	//textureFilter->SetInputImageMaximum(255);
	//textureFilter->SetNumberOfBinsPerAxis( 8 );
	//textureFilter->Update();

	//std::stringstream ss;
	//Write(textureFilter->GetOutput(),"haralick.nii");

	// Threshold top 80% of haralick as stool
	
	//// Show classification with hard threshold at 200
	//PixelType tissueStoolThreshold = 200;
	//ApplyThresholdRules( input, input->GetLargestPossibleRegion(), gradientMagnitude, vmap, colon, tissueStoolThreshold );
	//Write(vmap,"vmap200.nii");

	//vmap->FillBuffer(Unclassified);

	//// Show classification with global otsu threshold
	//ApplyThresholdRules( input, input->GetLargestPossibleRegion(), gradientMagnitude, vmap, colon, otsuGlobal );
	//Write(vmap,"vmapGlobalOtsu.nii");

	//vmap->FillBuffer(Unclassified);

	FloatIteratorType gmIt(gradientMagnitude,region);
	IteratorType inputIt(input,region);
	IteratorType thIt(thresholdImage,region);
	ByteIteratorType colonIt(colon,region);
	VoxelIteratorType vmapIt(vmap,region);

	// Show classification with adaptive threshold
	for (gmIt.GoToBegin(), inputIt.GoToBegin(), colonIt.GoToBegin(), vmapIt.GoToBegin(), thIt.GoToBegin(); !gmIt.IsAtEnd(); ++gmIt, ++inputIt, ++colonIt, ++vmapIt, ++thIt)
	{
		if (colonIt.Get() != 0)
		{
			PixelType I = inputIt.Get();
			PixelType T = thIt.Get();
			float G = gmIt.Get();

			if ( ( I >= T && G < 0.8*I ) || I > 1000 )
			{
				vmapIt.Set(Stool);
			} else if ( I <= -700 ) {
				vmapIt.Set(Air);
			} else if ( I < T && I > -300 && G <= 300 ) {
				vmapIt.Set(Tissue);
			}
		}
	}

	Write(vmap,"vmapLocal.nii");	
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
		if (partialIt.Get()[1] > 0.05)
		{
			ptbIt.Set(255);
		}
	}

	Write(ptb,"partialTissueBinaryMask.nii");

	// keep only largest components
	ptb = BinaryShapeOpen(ptb,"Size",200);
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
	air = BinaryDilate(air,3);
	Write(air,"airEdgesDilated.nii");

	// smooth only near edge tissue interface
	typedef itk::DiscreteGaussianImageFilter<ImageType,ImageType> SmoothType;
	SmoothType::Pointer smoother = SmoothType::New();
	smoother->SetInput(input);
	
	/*ArrayType sigmaArray;
	sigmaArray[0] = 0.7;
	sigmaArray[1] = 0.7;
	sigmaArray[2] = 0;*/

	smoother->SetVariance(.6*.6);
	smoother->Update();
	ImageType::Pointer inputSmooth = smoother->GetOutput();

	Write(inputSmooth,"inputSmooth.nii");

	IteratorType inputSmoothIt(inputSmooth,region);
	VoxelIteratorType vmapIt(vmap,region);

	airIt = ByteIteratorType(air,region);
	inputIt = IteratorType(input,region);

	// Smooth only the surrounding air next to the tissue
	for (inputIt.GoToBegin(), airIt.GoToBegin(), inputSmoothIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputSmoothIt, ++airIt)
	{
		if (airIt.Get() != 0 /*&& vmapIt.Get() != TissueAir*/) // do not smooth tissue air
		{
			if (inputIt.Get() < -400)
			{
				inputIt.Set( inputSmoothIt.Get() );
			}
			
		}
	}

	Write(input,"inputPostSmoothNearTissueOnly.nii");

	// Smooth air and tissue together at a lower sigma
	smoother = SmoothType::New();
	smoother->SetInput(input);
	smoother->SetVariance(0.8*0.8);
	smoother->Update();
	inputSmooth = smoother->GetOutput();
	
	inputSmoothIt = IteratorType(inputSmooth,region);

	for (inputIt.GoToBegin(), airIt.GoToBegin(), inputSmoothIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputSmoothIt, ++airIt)
	{
		if (airIt.Get() != 0 /*&& vmapIt.Get() != TissueAir*/) // do not smooth tissue air
		{
			inputIt.Set( inputSmoothIt.Get() );
		}
	}

	// Median filter in air/tissue region
	//inputSmooth = Median(input,1);

	//for (inputIt.GoToBegin(), airIt.GoToBegin(), inputSmoothIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputSmoothIt, ++airIt)
	//{
	//	if (airIt.Get() != 0 /*&& vmapIt.Get() != TissueAir*/) // do not smooth tissue air
	//	{
	//		inputIt.Set( inputSmoothIt.Get() );
	//	}
	//}

	Write(input,"inputPostSmoothAll.nii");


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

	//Write2(inputOriginal,"output.nii");

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

	BinaryDilate(maskDilated,3);
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

ByteImageType::Pointer SegmentColon(ImageType::Pointer &input)
{
	// Get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	
	// Find body
	ByteImageType::Pointer body = AllocateByteImage(input);
	ByteIteratorType bit(body,region);
	IteratorType it(input,region);

	for (it.GoToBegin(),bit.GoToBegin(); !it.IsAtEnd(); ++it,++bit)
	{
		if (it.Get() > -250 && it.Get() < 200)
		{
			bit.Set(255);
		}
	}

	Write(body,"bodyThreshold.nii");

	// Median to remove table
	body = Median(body,4);
	Write(body,"bodyMedian.nii");

	// Fill holes
	BinaryFillHoles2D(body);
	Write(body,"bodyFilledHoles.nii");

	// Shrink body so air on border is removed
	body = BinaryErode(body,4);
	bit = ByteIteratorType(body,region);

	Write(body,"bodyEroded.nii");
	
	// Find air
	ByteImageType::Pointer air = AllocateByteImage(input);
	ByteIteratorType ait(air,region);

	for (it.GoToBegin(), ait.GoToBegin(), bit.GoToBegin(); !it.IsAtEnd(); ++it, ++ait, ++bit)
	{
		if (bit.Get() != 0)
		{
			if (it.Get() < -600)
			{
				ait.Set(255);
			}
		}
	}

	// Remove lungs 
	// Detect lungs using region growing with seeds from first slice
	// Assumption: only low intensity object at z=0 is lung
	typedef itk::ConnectedThresholdImageFilter<ByteImageType, ByteImageType> ConnectedThresholdImageFilterByteType;
	ConnectedThresholdImageFilterByteType::Pointer lungGrower = ConnectedThresholdImageFilterByteType::New();
	lungGrower->SetInput( air );
	lungGrower->SetLower(255);
	lungGrower->SetUpper(255);
	lungGrower->SetReplaceValue(255);

	// Set air regions within dilated area as seeds of region growing
	for (ait.GoToBegin(); !ait.IsAtEnd(); ++ait)
	{
		ByteImageType::IndexType idx = ait.GetIndex();

		if (idx[2] == 0)
		{
			if (ait.Get() == 255)
			{
				lungGrower->AddSeed(idx);
			}
		}
	}

	lungGrower->Update();
	ByteImageType::Pointer lung = lungGrower->GetOutput();
	ByteIteratorType lit(lung,region);

	// Remove lung from air mask
	for (ait.GoToBegin(), lit.GoToBegin(); !ait.IsAtEnd(); ++ait, ++lit)
	{
		if (lit.Get() == 255)
			ait.Set( 0 );
	}

	lung.~SmartPointer();
	lungGrower.~SmartPointer();

	//Write(air,"air.nii");

	// Find bone
	PixelType bone_threshold = 250;
	PixelType stool_threshold = 180;

	typedef itk::ConnectedThresholdImageFilter< ImageType, ByteImageType> ConnectedThresholdImageFilterType;
	ConnectedThresholdImageFilterType::Pointer grower = ConnectedThresholdImageFilterType::New();
	grower->SetInput( input );
	grower->SetLower( bone_threshold );
	grower->SetUpper(1500); //1200
	grower->SetReplaceValue(255);

	// Set tagged regions within dilated area as seeds of region growing
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		ImageType::IndexType idx = it.GetIndex();
	
		if (idx[2] == 0)
		{
			if (it.Get() > bone_threshold )
			{
				grower->AddSeed(idx);
			}
		}
	}

	grower->Update();
	ByteImageType::Pointer bone = grower->GetOutput();
	
	grower.~SmartPointer();

	// Dilate bone
	bone = BinaryDilate(bone,4);
	ByteIteratorType boneit(bone,region);

	//Write(bone,"bone.nii");

	// Find all tagged with high intensity, adjacent to air, but are not bone
	// Get dilated air region
	ByteImageType::Pointer airNear = BinarySubtract( BinaryDilate(air,3) , air );
	ByteIteratorType anit(airNear,region);


	// Get all high intensity regions
	ByteImageType::Pointer tagged = BinaryThreshold(input,stool_threshold);

	// Remove bone
	tagged = BinarySubtract(tagged,bone);

	// Smooth
	tagged = Median(tagged,3);
	ByteIteratorType tit(tagged,region);

	//Write(tagged,"tagged.nii");

	// Remove all tagged components which are not close to air
	typedef itk::ConnectedComponentImageFilter<ByteImageType, LabelImageType> CCFilterType;
	CCFilterType::Pointer ccFilter = CCFilterType::New();
	ccFilter->SetInput(tagged);
	ccFilter->SetMaskImage(body);
	ccFilter->Update();
	LabelImageType::Pointer tcc = ccFilter->GetOutput();
	ccFilter.~SmartPointer();
	
	Relabel(tcc);

	LabelIteratorType ccit(tcc,region);

	typedef itk::MinimumMaximumImageCalculator<LabelImageType> CalcType;
	CalcType::Pointer calc = CalcType::New();
	calc->SetImage(tcc);
	calc->ComputeMaximum();

	unsigned int numComponents = calc->GetMaximum();

	std::vector<bool> compNearAirVec;
	compNearAirVec.resize(numComponents);

	for (int i=0; i<numComponents; i++)
		compNearAirVec[i] = false;

	for (ccit.GoToBegin(),anit.GoToBegin(); !ccit.IsAtEnd(); ++ccit,++anit)
	{
		if ( ccit.Get() > 0 )
		{
			if (anit.Get() != 0)
			{
				compNearAirVec[ ccit.Get() - 1] = true;
			}
		}
	}

	for (ccit.GoToBegin(), tit.GoToBegin(); !ccit.IsAtEnd(); ++ccit, ++tit)
	{
		if ( ccit.Get() > 0)
		{
			if ( compNearAirVec[ ccit.Get() - 1] == false )
			{
				tit.Set(0);
			}
		}
	}

	//Write(tagged,"taggedNearAir.nii");

	tcc.~SmartPointer();

	/*

	grower = ConnectedThresholdImageFilterType::New();
	grower->SetInput( input );
	grower->SetLower( stool_threshold );
	grower->SetUpper(1500); //1200
	grower->SetReplaceValue(255);

	// Set tagged regions within dilated area as seeds of region growing
	for (it.GoToBegin(), anit.GoToBegin(); !it.IsAtEnd(); ++it, ++anit)
	{
		ImageType::IndexType idx = it.GetIndex();

		if (it.Get() > stool_threshold)
		{
			if (anit.Get() != 0)
			{
				grower->AddSeed(it.GetIndex());
			}
		}
	}

	grower->Update();
	ByteImageType::Pointer tagged = grower->GetOutput();
	Write(tagged,"tagged.nii");
	
	grower.~SmartPointer();

	// Mask tagged to remove bone
	tagged = Mask(tagged, BinaryInvert(bone) );
	Write(tagged,"taggedNoBone.nii");

	// Find all tagged with high intensity but are not bone
	ByteImageType::Pointer tagged = AllocateByteImage(input);
	ByteIteratorType tit(tagged,region);

	for (it.GoToBegin(), tit.GoToBegin(), bit.GoToBegin(), boneit.GoToBegin(); !it.IsAtEnd(); ++it, ++tit, ++bit, ++boneit)
	{
		if (bit.Get() != 0)
		{
			if ( boneit.Get() == 0 )
			{
				if (it.Get() > bone_threshold)
				{
					tit.Set(255);
				}
			}
		}
	}


	*/

	// Add air and tagged
	tagged = BinaryOr(tagged,air);

	// Median smooth
	tagged = Median(tagged, 3);

	//Write(tagged,"taggedSmooth.nii");

	// Dilate
	tagged = BinaryDilate(tagged,12);

	return tagged;
}




/*********************************************************
- Load image
- Shift input so that air is 0 HU
- Segment colon
- Mask input with colon
- Crop input and colon in XY plane
*********************************************************/
void Setup(ImageType::Pointer &inputOriginal, ImageType::Pointer &input, ByteImageType::Pointer &colon)
{
	//----------------------------------------------
	// Load image
	//----------------------------------------------
	/*std::vector<std::string> datasetArray = explode( "\\", dataset );
	std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;*/

	//if (truncateOn)
	//{
	//	std::cout << "Truncating data" << std::endl;
	//	std::cout << "Slices: " << truncateArray[0] << " to " << truncateArray[1] << std::endl;
	//	
	//	//std::stringstream ss;
	//	//ss << dsname << "_" << truncateArray[0] << "_" << truncateArray[1];
	//	//dsname = ss.str();
	//}

	// Set writer prefix
	//note = dsname;

	// Make output directory
	/*outputDirectory += "\\" + dsname;

	itksys::SystemTools::MakeDirectory( outputDirectory.c_str() );

	outputDirectory += "\\dcm";

	std::cout << "Output Directory: " << outputDirectory << std::endl;

	itksys::SystemTools::MakeDirectory( outputDirectory.c_str() );*/

	
	//// Load dicom files
	//if (!truncateOn)
	//{
	//	inputOriginal = ReadDicom < ImageType > ( dataset );
	//} else {
	//	inputOriginal = ReadDicom < ImageType > ( dataset, truncateArray[0], truncateArray[1] );
	//}

	//// Setup output file names
	//for (int i=0; i<names.size(); i++)
	//{
	//	names[i] = outputDirectory + "/" + names[i];
	//}

	//seriesWriter->SetFileNames(names);

	//inputOriginal = ReadDicom <ImageType> (dataset, outputDirectory);

	//Write(inputOriginal,"inputOriginal.nii");


	//// Make output directory
	//std::vector<std::string> datasetArray = explode( "\\", inDir );
	//std::string dsname = datasetArray[ datasetArray.size() - 2 ];
	//std::cout << "Dataset: " << dsname << std::endl;
	//outDir += "\\" + dsname + "\\dcm";
	//std::cout << "Output Directory: " << outDir << std::endl;

	//const char * outputDirectory = outDir.c_str();
	//
	//itksys::SystemTools::MakeDirectory( outputDirectory );

	//// Create series writer
	//typedef itk::ImageSeriesWriter<ImageType,ImageType2D> SeriesWriterType;
	//SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	//
	//itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
	//seriesWriter->SetImageIO( dicomIO );

	//// Change names to output folder names
	//for (int i=0; i<names.size(); i++)
	//{
	//	std::vector<std::string> namesArray = explode("/",names[i]);
	//	names[i] = outDir + "\\" + namesArray[ namesArray.size() - 1 ];
	//}

	//seriesWriter->SetFileNames( names );

	//std::cout << "Number of filenames: " << names.size() << std::endl;
	//std::cout << "Number of dictionary entries: " << metaDataDictionaryArray.size() << std::endl;

	//for(int i=0; i<metaDataDictionaryArray.size(); i++)
	//{
	//	metaDataDictionaryArray[i] = &(dicomIO->GetMetaDataDictionary());
	//}*/

	//seriesWriter->SetMetaDataDictionaryArray( &metaDataDictionaryArray );//&metaDataDictionaryArray );
	//seriesWriter->SetInput( inputOriginal );

	//try
	//{
	//	seriesWriter->Update();
	//}
	//catch( itk::ExceptionObject & excp )
	//{
	//	std::cerr << "Exception thrown while writing the series " << std::endl;
	//	std::cerr << excp << std::endl;
	//	system("pause");
	//	return EXIT_FAILURE;
	//}





















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
	/*typedef itk::ColonSegmentationFilter< ImageType, ByteImageType > ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer colonSegmenter = ColonSegmentationFilterType::New();
	colonSegmenter->SetInput( inputOriginal );
	colonSegmenter->SetOutputForegroundValue( 255 );
	colonSegmenter->SetOutputBackgroundValue( 0 );
	colonSegmenter->SetPrintImages(true);

	if ( truncateOn )
		colonSegmenter->SetRemoveBoneLung( false );

	colonSegmenter->Update();
	colon = colonSegmenter->GetOutput();
	Write(colon,"colonFull.nii");
	BinaryFillHoles2D(colon);
	Write(colon,"colonFullFilled.nii");*/

	colon = SegmentColon(inputOriginal);

	std::cout << "Colon segmented: " << (double) (clock() - init) / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;

	//Write(colon,"colonFull.nii");
	
	BinaryFillHoles2D(colon);
	Write(colon,"colonFilled.nii");

	//----------------------------------------------
	// Mask input with colon
	//----------------------------------------------
	/*typedef itk::MaskImageFilter< ImageType, ByteImageType, ImageType > MaskImageFilterType;
	MaskImageFilterType::Pointer masker = MaskImageFilterType::New();
	masker->SetInput1( inputOriginal );
	masker->SetInput2( colon );
	masker->SetOutsideValue( BACKGROUND );
	masker->Update();
	input = masker->GetOutput();*/

	//----------------------------------------------
	// Calculate gradient magnitude
	//----------------------------------------------
	
	// 3d
	/*typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
	
	GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientMagnitudeFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientMagnitudeFilter->SetInput( inputOriginal );
	gradientMagnitudeFilter->SetSigma( inputOriginal->GetSpacing()[0] );
	gradientMagnitudeFilter->Update();
	gradientMagnitude = gradientMagnitudeFilter->GetOutput();*/

	// 2d
	/*
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType2D,FloatImageType2D > GradientMagnitudeRecursiveGaussianImageFilterType2D;
	typedef itk::SliceBySliceImageFilter< ImageType, FloatImageType, GradientMagnitudeRecursiveGaussianImageFilterType2D > SliceGradientFilterType;

	GradientMagnitudeRecursiveGaussianImageFilterType2D::Pointer gm2DF = GradientMagnitudeRecursiveGaussianImageFilterType2D::New();
	gm2DF->SetSigma( inputOriginal->GetSpacing()[0] );
	
	SliceGradientFilterType::Pointer sliceGradientFilter = SliceGradientFilterType::New();
	sliceGradientFilter->SetInput(input);
	sliceGradientFilter->SetFilter(gm2DF);
	sliceGradientFilter->Update();
	gradientMagnitude = sliceGradientFilter->GetOutput();
	*/

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
	
	ByteIteratorType colonIt(colon,REGION);

	for (colonIt.GoToBegin(); !colonIt.IsAtEnd(); ++colonIt)
	{
		if ( colonIt.Get() != 0 )
		{
			ByteImageType::IndexType idx = colonIt.GetIndex();

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

	colon = Crop(colon,extractRegion);

	input = Crop(inputOriginal,extractRegion);

	/*
	typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> RegionOfInterestImageFilterType;
	RegionOfInterestImageFilterType::Pointer cropper = RegionOfInterestImageFilterType::New();
	cropper->SetInput( input );
	cropper->SetRegionOfInterest( extractRegion );
	cropper->Update();
	input = cropper->GetOutput();
	*/

	/*
	typedef itk::RegionOfInterestImageFilter<FloatImageType,FloatImageType> RegionOfInterestImageFilterFloatType;
	RegionOfInterestImageFilterFloatType::Pointer cropperFloat = RegionOfInterestImageFilterFloatType::New();
	cropperFloat->SetInput( gradientMagnitude );
	cropperFloat->SetRegionOfInterest( extractRegion );
	cropperFloat->Update();
	gradientMagnitude = cropperFloat->GetOutput();
	*/

	/*cropperFloat = RegionOfInterestImageFilterFloatType::New();
	cropperFloat->SetInput( gradientMagnitude2D );
	cropperFloat->SetRegionOfInterest( extractRegion );
	cropperFloat->Update();
	gradientMagnitude2D = cropperFloat->GetOutput();*/

	/*
	typedef itk::RegionOfInterestImageFilter<ByteImageType,ByteImageType> RegionOfInterestImageFilterByteType;
	RegionOfInterestImageFilterByteType::Pointer cropperByte = RegionOfInterestImageFilterByteType::New();
	cropperByte->SetInput( colon );
	cropperByte->SetRegionOfInterest( extractRegion );
	cropperByte->Update();
	colon = cropperByte->GetOutput();
	*/

	// Set cropped region globally
	REGION = colon->GetLargestPossibleRegion();

	Write(input,"input.nii");
	Write(colon,"colon.nii");

	//Write(gradientMagnitude,"gradientMagnitudeSmoothed.nii");
	//Write(gradientMagnitude2D,"gradientMagnitudeSmoothed2D.nii");
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

	PixelType tissueStoolThreshold = 200;

	//----------------------------------------------
	// Apply initial threshold rules
	//----------------------------------------------

	// Use smoothed input to eliminate noise
	ApplyThresholdRules( Median(input,1), input->GetLargestPossibleRegion(), gradientMagnitude, vmap, colon, tissueStoolThreshold );

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

			if ( ( I >= tissueStoolThreshold && G < 0.8*I ) || I > 800 )
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

std::vector<float> ComputeLogSlope(std::vector<float> x, std::vector<float> y)
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

	// output slope and intercept
	std::vector<float> out;
	out.resize(2);
	
	if ( ( s0*s2 - s1*s1 ) == 0 )
	{
		out[0] = 0;
		out[1] = 0;
	} else {
		out[0] = ( s0*t1 - s1*t0 ) / (s0*s2 - s1*s1);
		out[1] = ( s2*t0 - s1*t1 ) / (s0*s2 - s1*s1);
	}

	return out;
}

std::vector<FloatImageType::Pointer> RescaledRange(ByteImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius)
{
	// Get region
	ByteImageType::RegionType region = input->GetLargestPossibleRegion();

	// Get mask
	ByteIteratorType maskIt(mask,region);

	// Make map image
	ByteImageType2D::Pointer map = ByteImageType2D::New();
	
	ByteImageType2D::RegionType mapRegion;
	
	ByteImageType2D::IndexType mapIndex;
	ByteImageType2D::SizeType mapSize;


	for (int i=0; i<2; i++)
	{
		mapIndex[i] = 0;
	}

	mapSize[0] = 2*radius+1;
	mapSize[1] = 2*radius+1;

	mapRegion.SetIndex(mapIndex);
	mapRegion.SetSize(mapSize);

	map->SetRegions(mapRegion);

	map->Allocate();
	map->FillBuffer(0);

	Write(map,"map.nii");

	ByteIteratorType2D mapIt(map,mapRegion);

	// Setup neighborhood and fill map with distance squared
	typedef itk::Neighborhood<BytePixelType,Dimension> NeighborhoodType;
	NeighborhoodType hood;
	
	ImageType::SizeType radiusSize;
	radiusSize.Fill(0);
	radiusSize[0] = radius;
	radiusSize[1] = radius;

	hood.SetRadius(radiusSize);

	std::cout << "hood.Size(): " << hood.Size() << std::endl;

	typedef ByteImageType::OffsetType OffsetType;
	
	int count=0;

	for (mapIt.GoToBegin(); !mapIt.IsAtEnd(); ++mapIt)
	{
		OffsetType off = hood.GetOffset(count);
		
		PixelType d = 0;

		for (int i=0; i<2; i++)
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

	Write(map,"map2.nii");

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
	
	typedef itk::LabelMapToLabelImageFilter<LabelMapType,ByteImageType2D> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToImageFilter->SetInput(relabeler->GetOutput());
	labelMapToImageFilter->Update();

	map = labelMapToImageFilter->GetOutput();
	mapIt = ByteIteratorType2D(map,mapRegion);

	labelMap.~SmartPointer();

	Write(map,"map3.nii");

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

	// Allocate output images
	std::vector<FloatImageType::Pointer> outVector;
	outVector.resize(2);

	for (int i=0; i<2; i++)
	{
		outVector[i] = FloatImageType::New();
		outVector[i]->SetRegions(region);
		outVector[i]->CopyInformation(input);
		outVector[i]->Allocate();
		outVector[i]->FillBuffer(0);
		
	}

	FloatIteratorType oit1(outVector[0],region);
	FloatIteratorType oit2(outVector[1],region);

	// Iterate through image
	ByteImageType::SizeType rad;
	rad.Fill(0);
	rad[0] = radius;
	rad[1] = radius;
	
	typedef itk::NeighborhoodIterator<ByteImageType> NeighborhoodIteratorType;
	NeighborhoodIteratorType nit(rad,input,region);

	for (nit.GoToBegin(), oit1.GoToBegin(), oit2.GoToBegin(), maskIt.GoToBegin(); !nit.IsAtEnd(); ++nit, ++oit1, ++oit2, ++maskIt)
	{
		if (maskIt.Get() != 0)
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

			// Get slope and y intercept
			std::vector<float> line = ComputeLogSlope(dv,rv);

			oit1.Set(line[0]);
			oit2.Set(line[1]);
		}
	}

	/*std::stringstream ss;
	ss << "rr" << radius << ".nii";
	WriteITK <FloatImageType2D> (out,ss.str());*/
							
	return outVector; 									
} 		

ImageType::Pointer LocalMIP(ImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int Radius)
{
	// get region
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	
	// get mask
	ByteIteratorType maskIt(mask,region);

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
	ImageType::Pointer maxOut = AllocateImage(input);
	
	IteratorType maxIt(maxOut,region);

	for (it.GoToBegin(), maxIt.GoToBegin(), maskIt.GoToBegin(); !it.IsAtEnd(); ++it, ++maxIt, ++maskIt)
	{
		if (maskIt.Get() != 0)
		{
			ShapedNeighborhoodIteratorType::ConstIterator ci;
			
			PixelType max = itk::NumericTraits<PixelType>::NonpositiveMin();
			
			for (ci = it.Begin(); ci != it.End(); ci++)
			{
				if (ci.Get() > max)
					max = ci.Get();
			}

			maxIt.Set( max );
		}
	}

	std::stringstream ss;
	ss << "localMax_radius_" << Radius << ".nii";
	Write(maxOut,ss.str());
					
	return maxOut; 									
}