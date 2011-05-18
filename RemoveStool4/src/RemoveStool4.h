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
#include <itkAddConstantToImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>


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

typedef short			                                                            PixelType;
typedef unsigned char																BytePixelType;
typedef itk::BinaryBallStructuringElement<PixelType, 3>								StructuringElementType;
typedef itk::CovariantVector<float,3>												VectorType;
typedef itk::FixedArray<float,3>													ArrayType;


typedef itk::Image<PixelType, 3>													ImageType;
typedef itk::Image<BytePixelType, 3>												ByteImageType;
//typedef itk::Image< short int, 3>													ShortImageType;
typedef itk::Image< float, 3>														FloatImageType;
typedef itk::Image< VoxelType, 3>													VoxelImageType;
typedef itk::Image< VectorType, 3>													VectorImageType;
typedef itk::Image< ArrayType, 3>													ArrayImageType;
typedef itk::BSplineInterpolateImageFunction<FloatImageType>::ContinuousIndexType   ContinuousIndexType;
	
typedef itk::ImageRegionIteratorWithIndex<ImageType>							    IteratorType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType>							ByteIteratorType;
//typedef itk::ImageRegionIteratorWithIndex<ShortImageType>							ShortIteratorType;
typedef itk::ImageRegionIteratorWithIndex<FloatImageType>							FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex<VoxelImageType>							VoxelIteratorType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType>							VectorIteratorType;
typedef itk::ImageRegionIteratorWithIndex<ArrayImageType>							ArrayIteratorType;

void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradient_magnitude);
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon);
void ApplyThresholdRules( ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissue_stool_threshold );
void Dilate(ByteImageType::Pointer &img, unsigned int radius);

vnl_matrix<float> GetNeighbor(ArrayImageType::Pointer partialVector, ImageType::IndexType index);

// Global vars
ImageType::RegionType OLDREGION;
ImageType::RegionType REGION;
bool write_num = true;
int write_count = 1;
bool truncateOn = true;
unsigned int truncateArray[2] = {85,90};
std::string note;
PixelType BACKGROUND = 0;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
std::ofstream debug;