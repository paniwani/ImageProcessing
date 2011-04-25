/****************************************************************************************   
   1.  Convolve the image with a suitable statistical operator, i.e. the mean or median.
   2. Subtract the original from the convolved image.
   3. Threshold the difference image with C.
   4. Invert the thresholded image. 
****************************************************************************************/

#include <itkImage.h>
#include <utils.h>
#include <itkMeanImageFilter.h>
#include <itkPNGImageIO.h>

int main()
{
    unsigned int rad = 3;
    unsigned int C = 7;    

	typedef itk::ImageFileReader<ByteImageType2D> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO( itk::PNGImageIO::New() );
	reader->SetFileName( "son1.png" );

	try
	{
		reader->Update();
	} catch( itk::ExceptionObject & err ) {
		std::cerr << "Error reading image: " << err << std::endl;
		system("pause");
		return EXIT_FAILURE;
	}

	ByteImageType2D::Pointer input = reader->GetOutput();

	ByteImageType2D::RegionType region = input->GetLargestPossibleRegion();

    // 1.  Convolve the image with a suitable statistical operator, i.e. the mean or median.
    
    typedef itk::MeanImageFilter<ByteImageType2D,ByteImageType2D> MeanImageFilterType;
    MeanImageFilterType::Pointer meanFilter = MeanImageFilterType::New();
    meanFilter->SetInput(input);
    
    ByteImageType2D::SizeType radius;
	radius.Fill(0);
    radius[0] = rad;
    radius[1] = rad;

    meanFilter->SetRadius(radius);
    meanFilter->Update();

    ByteImageType2D::Pointer input_diff = meanFilter->GetOutput();

	WriteITK <ByteImageType2D> (input_diff,"mean.nii"); 

    ByteIteratorType2D input_iter(input,region);
    ByteIteratorType2D input_diff_iter(input_diff,region);

	// 3. Threshold the difference image with C.

    for (input_iter.GoToBegin(), input_diff_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter,++input_diff_iter)
    {
		unsigned char diff = 0;

		if ( input_diff_iter.Get() > input_iter.Get() )
			diff = input_diff_iter.Get() - input_iter.Get();

		if ( diff > C )
		{
			input_diff_iter.Set(255);
		} else {
			input_diff_iter.Set(0);
		}
    }

    WriteITK <ByteImageType2D> (input_diff,"local_thresold.nii");    
}
