#include <itkImage.h>

// IO 
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

typedef float			                                                            PixelType;
typedef unsigned char																BytePixelType;

typedef itk::Image<PixelType, 3>													ImageType;
typedef itk::Image<BytePixelType, 3>												ByteImageType;
typedef itk::Image< short int, 3>													ShortImageType;
typedef itk::Image< float, 3>														FloatImageType;
typedef itk::Image< VoxelType, 3>													VoxelImageType;

typedef itk::ImageRegionIteratorWithIndex<ImageType>							    IteratorType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType>							ByteIteratorType;
typedef itk::ImageRegionIteratorWithIndex<ShortImageType>							ShortIteratorType;
typedef itk::ImageRegionIteratorWithIndex<FloatImageType>							FloatIteratorType;
typedef itk::ImageRegionIteratorWithIndex<VoxelImageType>							VoxelIteratorType;

void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon);
void SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon);

// Global vars
ImageType::RegionType region;
bool write_num = true;
int write_count = 1;
bool truncate_on = true;
unsigned int truncate_ar[2] = {85,90};
std::string note;