#include "RemoveStool3.h"

const double ALPHA = 0.5;
const double BETA = 0.3;
const double GAMMA = 0.3;
const double ETA = 0.2;

int writeCount = 1;

std::vector<std::string> explode( const std::string &delimiter, const std::string &str);

int main(int argc, char * argv[])
{

	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory\n";
		system("PAUSE");
		return EXIT_FAILURE;
	}

	std::string dataset = argv[1];

	std::vector<std::string> datasetArr = explode( "/", dataset );
	std::string ds = datasetArr[ datasetArr.size() - 2 ];
	std::cout << "Dataset: " << ds << std::endl;

	// Read and write dicom input
	ImageType::Pointer input = ReadDicom( dataset );

	WriteITK(input, ds+".hdr");

	// Get input info
	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SpacingType spacing = input->GetSpacing();

	// Threshold to compute only in tagged regions
	IteratorTypeFloat4WithIndex input_iter(input,region);

	ByteImageType::Pointer tagged = AllocateNewByteImage(region);
	IteratorTypeByteWithIndex tagged_iter(tagged,region);

	for (input_iter.GoToBegin(), tagged_iter.GoToBegin();
		 !input_iter.IsAtEnd() && !tagged_iter.IsAtEnd();
		 ++input_iter, ++tagged_iter)
	{
		if (input_iter.Get() > 200) 
		{
			tagged_iter.Set(1);
		} else {
			tagged_iter.Set(0);
		}
	}

	// Create response image
	ImageType::Pointer Himage = AllocateNewImage(region);
	Himage->FillBuffer(0.0);
	IteratorTypeFloat4WithIndex Himage_iter(Himage,region);

	// Compute hessian across sigma scales
	double sigma[2] = {0.8*spacing[0], 2*spacing[0]};

	for (int k=0; k<2; k++)
	{
		// Compute smoothed Hessian
		HessianGaussianFilterType::Pointer hessianFilter = HessianGaussianFilterType::New();
		hessianFilter->SetInput(input);
		hessianFilter->SetNormalizeAcrossScale(true);
		hessianFilter->SetSigma(sigma[k]);
		hessianFilter->Update();
		itk::ImageRegionConstIterator<HessianGaussianFilterType::OutputImageType> hessian_iter(hessianFilter->GetOutput(),region);

		hessian_iter.GoToBegin();
		tagged_iter.GoToBegin();
		Himage_iter.GoToBegin();

		while (!hessian_iter.IsAtEnd()) 
		{
			EigenValueArrayType lambda;

			if (tagged_iter.Get() == 1) { // compute only in tagged region
				
				hessian_iter.Get().ComputeEigenValues(lambda);
				std::sort(lambda.Begin(),lambda.End(),OrderByMagnitude);

				if ( lambda[2] < 0 )
				{
					float val = vnl_math_max( fRut(lambda, ALPHA, BETA, GAMMA), fCup(lambda, ETA) );
					
					if ( val > Himage_iter.Get() )
						Himage_iter.Set(val);
				}
			}

			++hessian_iter;
			++tagged_iter;
			++Himage_iter;
		}

		std::stringstream ss;
		ss << ds << "_H.hdr";
		WriteITK(Himage, ss.str());

	}

	return 0;
}

void WriteITK(ImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
    typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

bool OrderByMagnitude (double a, double b) {
	return ( abs(a) < abs(b) );
}

bool OrderByValue (double a, double b) {
	return ( a < b );
}

bool OrderByValueDesc (double a, double b) {
	return ( a > b );
}

double fA(EigenValueArrayType lambda, double alpha) {
	double Ra;
	Ra = abs(lambda[0])/sqrt(abs(lambda[1]*lambda[2]));
	return exp(-vnl_math_sqr(Ra)/(2*vnl_math_sqr(alpha)));
}

double fB(EigenValueArrayType lambda, double beta, double gamma) {
	double Rb;
	Rb = abs(lambda[1]/lambda[2]);
	return exp(-vnl_math_sqr(Rb - gamma)/(2*vnl_math_sqr(beta)));
}

double fC(double ev1, double ev2, double eta) {
	double Rc;
	Rc = abs(ev1)/abs(ev2);
	return 1.0 - exp(-vnl_math_sqr(Rc)/(2*vnl_math_sqr(eta)));
}

double fRut(EigenValueArrayType lambda, double alpha, double beta, double gamma) {
	return fA(lambda, alpha)*fB(lambda, beta, gamma);
}

double fCup(EigenValueArrayType lambda, double eta) {
	return fC(lambda[0], lambda[1], eta)*fC(lambda[1], lambda[2], eta);
}

void ReadITK(ImageType::Pointer &image, char * fileName) {
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	image = reader->GetOutput();
}

std::vector<std::string> explode( const std::string &delimiter, const std::string &str)
{
	std::vector<std::string> arr;

    int strleng = str.length();
    int delleng = delimiter.length();
    if (delleng==0)
        return arr;//no change

    int i=0; 
    int k=0;
    while( i<strleng )
    {
        int j=0;
        while (i+j<strleng && j<delleng && str[i+j]==delimiter[j])
            j++;
        if (j==delleng)//found delimiter
        {
            arr.push_back(  str.substr(k, i-k) );
            i+=delleng;
            k=i;
        }
        else
        {
            i++;
        }
    }
    arr.push_back(  str.substr(k, i-k) );
    return arr;
}

ImageType::Pointer ReadDicom( std::string path )
{	
	// Create reader
	itk::ImageSeriesReader<ImageType>::Pointer reader = itk::ImageSeriesReader<ImageType>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	// Truncate data
	//names.erase( names.begin()+10, names.end() );

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return 0;
	}

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI);
    orienter->SetInput(reader->GetOutput());
    orienter->Update();

    return orienter->GetOutput();

}

ByteImageType::Pointer AllocateNewByteImage(ImageType::RegionType fullRegion) {
    ByteImageType::Pointer newImage = ByteImageType::New();
    newImage->SetLargestPossibleRegion( fullRegion );
    newImage->SetBufferedRegion( fullRegion );
    newImage->SetRequestedRegion( fullRegion );
    newImage->Allocate();
    return newImage;
}

ImageType::Pointer AllocateNewImage(ImageType::RegionType fullRegion) {
    ImageType::Pointer newImage = ImageType::New();
    newImage->SetLargestPossibleRegion( fullRegion );
    newImage->SetBufferedRegion( fullRegion );
    newImage->SetRequestedRegion( fullRegion );
	newImage->Allocate();
    return newImage;
}