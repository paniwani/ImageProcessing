// IO 
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include <stdio.h>

// Math
#include <vnl/vnl_math.h>
#include <vnl/vnl_real_polynomial.h>
#include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_erf.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <complex>
#include <limits>

#include <itkColonSegmentationFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkImageDuplicator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkOtsuThresholdImageCalculatorModified.h>
#include <itkImageDuplicator.h>
#include <itkSubtractImageFilter.h>
#include <itkCannyEdgeDetectionImageFilterModified.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>


#define CDF_SIGMA 0.27

enum VoxelType {
	Stool=1,
	Air=2,
	Tissue=3,
	Unclassified=4,
	StoolAir=5,
	TissueAir=6,
	TissueStool=7,
    ThinStool=8
};

typedef unsigned short	                                                            PixelType;
typedef unsigned char																BytePixelType;
typedef itk::BinaryBallStructuringElement<PixelType, 3>								StructuringElementType;
typedef itk::CovariantVector<float,3>												VectorType;


typedef itk::Image<PixelType, 3>													ImageType;
typedef itk::Image<BytePixelType, 3>												ByteImageType;
typedef itk::Image< short int, 3>													ShortImageType;
typedef itk::Image< float, 3>														FloatImageType;
typedef itk::Image< VoxelType, 3>													VoxelImageType;
typedef itk::Image< VectorType, 3>													VectorImageType;
typedef itk::BSplineInterpolateImageFunction<FloatImageType>::ContinuousIndexType   ContinuousIndexType;
	
typedef itk::ImageRegionIteratorWithIndex<ImageType>							    IteratorType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType>							ByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex<ShortImageType>							ShortIteratorType;
typedef itk::ImageRegionIteratorWithIndex<FloatImageType>							FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex<VoxelImageType>							VoxelIteratorType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType>							VectorIteratorType;

void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon);
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon);
void ApplyThresholdRules( ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissue_stool_threshold );
void QuadraticRegression(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap, FloatImageType::Pointer &gradient_magnitude);
void Dilate(ByteImageType::Pointer &img, unsigned int radius);

// Global vars
ImageType::RegionType REGION;
bool write_num = true;
int write_count = 1;
bool truncate_on = true;
unsigned int truncate_ar[2] = {85,90};
std::string note;
short BACKGROUND = 0;