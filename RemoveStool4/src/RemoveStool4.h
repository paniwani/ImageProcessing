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
//#include <vnl/algo/vnl_rpoly_roots.h>
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
#include <itkOtsuThresholdImageCalculator.h>
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
#include <itkDirectionalGradientImageFilter.h>
#include <itkDirectionalGradientImageFilter2.h>
//#include <otbScalarImageToTexturesFilter2.h>
#include <time.h>
#include <itkMorphologicalDistanceTransformImageFilter.h>
#include <itkMorphologicalSignedDistanceTransformImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkCurvesLevelSetImageFilter.h>
#include <itkZeroCrossingImageFilter.h>
#include <itkGaussianBlurImageFunction.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkImageToHistogramFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionMultidimensionalSplitter.h>
#include <itkAdaptiveOtsuThresholdImageFilter.h>
#include <itkFastBilateralImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkConfidenceConnectedImageFilter.h>
#include <itkNeighborhood.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkRelabelLabelMapFilter.h>

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
typedef  itk::FixedArray< double, 3 >												EigenValueArrayType;

typedef itk::Image<unsigned long,1> ImageType1D;
typedef itk::ImageRegionIteratorWithIndex<ImageType1D> IteratorType1D;
typedef itk::Image<PixelType, 3>													ImageType;
typedef itk::Image<BytePixelType, 3>												ByteImageType;
typedef itk::Image<PixelType, 2>													ImageType2D;
typedef itk::Image<BytePixelType, 2>												ByteImageType2D;
typedef itk::Image<float, 2>														FloatImageType2D;
typedef itk::Image<unsigned int,3>													LabelImageType;
	

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
typedef itk::ImageRegionIteratorWithIndex<LabelImageType>							LabelIteratorType;

void Setup(std::string dataset, ImageType::Pointer  &input_original, ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradient_magnitude);
PixelType SingleMaterialClassification(ImageType::Pointer &input, FloatImageType::Pointer &gradient_magnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon);
void ApplyThresholdRules( ImageType::Pointer &input, ImageType::RegionType localRegion, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, PixelType tissueStoolThreshold );
void Dilate(ByteImageType::Pointer &img, unsigned int radius);
void DirectionalGradient(ImageType::Pointer &input, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap);
ImageType::Pointer Subtraction(ImageType::Pointer &input, ImageType::Pointer &inputOriginal, ByteImageType::Pointer &colon, ArrayImageType::Pointer &partial, VoxelImageType::Pointer &vmap);
//void TextureAnalysis(ImageType::Pointer &input);

void HeteroStoolRemoval(ImageType::Pointer &cOutput, ByteImageType::Pointer &colon, VoxelImageType::Pointer &vmap);
ImageType::Pointer Range(ImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius);
void LevelSet(ImageType::Pointer &input, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradient);
void FixATT(ImageType::Pointer &input, ArrayImageType::Pointer &partial, VoxelImageType::Pointer &vmap, ByteImageType::Pointer &colon, FloatImageType::Pointer &smax, PixelType tissueStoolThreshold);


FloatImageType::Pointer StandardDeviation(ImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius);
FloatImageType::Pointer StandardDeviation(FloatImageType::Pointer &input, ByteImageType::Pointer &mask, unsigned int radius);

void LocalThreshold(ImageType::Pointer &input, ByteImageType::Pointer &colon, FloatImageType::Pointer &gradientMagnitude, VoxelImageType::Pointer &vmap);

void AdaptiveThreshold(ImageType::Pointer &input, ByteImageType::Pointer &colon);

void ConnectedTest(ImageType::Pointer &input, FloatImageType::Pointer &gradientMagnitude, ByteImageType::Pointer &colon);

void TextureTest(ImageType::Pointer &input, ByteImageType::Pointer &colon);

float ComputeLogSlope(std::vector<float> x, std::vector<float> y);

FloatImageType::Pointer RescaledRange(ImageType::Pointer &input, unsigned int radius);

struct point{
	int intensity;
	int size;
	int order;
};

typedef struct point ptype;


// Global vars
ImageType::RegionType OLDREGION;
ImageType::RegionType REGION;
bool write_num = false;
int write_count = 1;
bool truncateOn = true;
//int truncateArray[2] = {85,-1};
unsigned int truncateArray[2] = {130,135};
//unsigned int truncateArray[2] = {85,90};
//unsigned int truncateArray[2] = {0,105};
std::string note;
PixelType BACKGROUND = 0;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
std::ofstream debug;