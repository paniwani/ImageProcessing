#include "itkExceptionObject.h"
#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbVectorImage.h"
#include "otbImageToVectorImageCastFilter.h"
#include "otbVectorRescaleIntensityImageFilter.h"
#include "otbScalarImageToTexturesFilter.h"

int main(int argc, char * argv[])
{
  // Parse command line parameters
  if (argc != 7)
    {
    std::cerr << "Usage: " << argv[0] << " <inputImage> ";
    std::cerr << " <outputImage> <outputRescaled> ";
    std::cerr << " <radius> <xOffset> <yOffset> ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  const char* infname   = argv[1];
  const char* outfname  = argv[2];
  const char* outprettyfname  = argv[3];

  const unsigned int radius  =  static_cast<unsigned int>(atoi(argv[4]));
  const unsigned int xOffset =  static_cast<unsigned int>(atoi(argv[5]));
  const unsigned int yOffset =  static_cast<unsigned int>(atoi(argv[6]));

  typedef double PixelType;
  const int Dimension = 2;
  typedef otb::Image<PixelType, Dimension> ImageType;

  typedef otb::ScalarImageToTexturesFilter<ImageType, ImageType> TexturesFilterType;
 
  typedef otb::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader  = ReaderType::New();

  reader->SetFileName(infname);

  TexturesFilterType::Pointer texturesFilter = TexturesFilterType::New();

  typedef ImageType::SizeType SizeType;
  SizeType sradius;
  sradius.Fill(radius);

  texturesFilter->SetRadius(sradius);

  typedef ImageType::OffsetType OffsetType;
  OffsetType offset;
  offset[0] =  xOffset;
  offset[1] =  yOffset;

  texturesFilter->SetOffset(offset);
 
  texturesFilter->SetInputImageMinimum(0);
  texturesFilter->SetInputImageMaximum(255);
 
  texturesFilter->SetInput(reader->GetOutput());
  texturesFilter->Update();

  for (int i=0; i < 8; i++)
  {	  
	  typedef otb::VectorImage<double, 2>        VectorImageType;
	  typedef otb::VectorImage<unsigned char, 2> PrettyVectorImageType;
	  typedef otb::ImageFileWriter<PrettyVectorImageType> WriterPrettyOutputType;

	  typedef otb::ImageToVectorImageCastFilter<ImageType, VectorImageType> VectorCastFilterType;
	  typedef otb::VectorRescaleIntensityImageFilter<VectorImageType, PrettyVectorImageType> RescalerOutputType;

	  RescalerOutputType::Pointer     outputRescaler     = RescalerOutputType::New();
	  WriterPrettyOutputType::Pointer prettyOutputWriter = WriterPrettyOutputType::New();
	  VectorCastFilterType::Pointer vectorCastFilter = VectorCastFilterType::New();
	  
	  vectorCastFilter->SetInput( texturesFilter->GetOutput( i ) );
	  
	  outputRescaler->SetInput(vectorCastFilter->GetOutput());

	  PrettyVectorImageType::PixelType min(1), max(1);
	  min.Fill(0);
	  max.Fill(255);

	  outputRescaler->SetOutputMinimum(min);
	  outputRescaler->SetOutputMaximum(max);

	  std::stringstream ss;
	  ss << "out" << i << ".png";

	  prettyOutputWriter->SetFileName(ss.str().c_str());

	  prettyOutputWriter->SetInput(outputRescaler->GetOutput());

	  prettyOutputWriter->Update();
  }

  return EXIT_SUCCESS;
}
