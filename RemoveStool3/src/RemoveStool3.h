#include <limits>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkGrayscaleMorphologicalClosingImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBSplineInterpolateImageFunction.h>
#include "itkChamferDistanceTransformImageFilter.h"
#include <itkGaussianBlurImageFunction.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_real_polynomial.h>
#include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_erf.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_qr.h>
#include <math.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkNeighborhoodIterator.h>
#include <stdio.h>
#include <itkSimpleContourExtractorImageFilter.h>
#include <itkCurvesLevelSetImageFilter.h>
#include <itkContinuousIndex.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkBinaryMedianImageFilter.h>

// Outdated includes
//#include "Typedef.h"
//#include "utils.h"
//#include "itkLocalRoughnessImageFilter.h"
//#include <NIH_Algo_TextureAnalysis.h>
//#include <DynArray.h>

// Regex
#include <boost/xpressive/xpressive.hpp>

#include <itkCastImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkFastMarchingImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkSimilarityIndexImageFilter.h>

#include <fstream>
#include <iostream>
#include <complex>

#define CDF_SIGMA 0.27
#define DEBUG true
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

//Pixel Type
typedef float                                                               PixelType;
typedef itk::ContinuousIndex< PixelType, 3 >								ContinuousIndexType;

typedef unsigned long														LabelType;
typedef itk::ShapeLabelObject< LabelType, 3 >								LabelObjectType;
typedef itk::LabelMap< LabelObjectType >									LabelMapType;

//CovariantVector Type
typedef itk::CovariantVector<PixelType,3>                                   CovariantVectorType;
typedef itk::Image<PixelType, 3>                                            ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType>                        IteratorTypeFloat4WithIndex;
typedef itk::Image<CovariantVectorType, 3>                                  ImageVectorType;

typedef itk::Image<BYTE, 3>													ByteImageType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType>                    IteratorTypeByteWithIndex;

typedef itk::Image<int, 3>													IntImageType;
typedef itk::ImageRegionIteratorWithIndex<IntImageType>						IteratorTypeIntWithIndex;


typedef itk::ImageRegionIterator<ImageVectorType>                           IteratorImageVectorType;
typedef itk::Image<VoxelType, 3>                                            VoxelTypeImage;
typedef itk::ImageRegionIteratorWithIndex<VoxelTypeImage>                   IteratorTypeVoxelType;
//Axillary Types
typedef itk::BinaryBallStructuringElement<PixelType, 3>                     StructuringElementType;
typedef itk::BSplineInterpolateImageFunction<ImageType, float>				InterpolationType;
typedef itk::GaussianBlurImageFunction<ImageType>                           GaussianFilterType;

//Image Filters
typedef itk::GradientImageFilter< ImageType >								GradientFilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType>             GradientMagnitudeFilterType;
typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, 
        ImageType, StructuringElementType >									ClosingFilterType;

typedef itk::GrayscaleMorphologicalClosingImageFilter<ByteImageType, 
        ByteImageType, StructuringElementType >								ByteClosingFilterType;


typedef itk::GradientAnisotropicDiffusionImageFilter<
		ImageType, ImageType >												AnisotropicFilterType;


typedef itk::ChamferDistanceTransformImageFilter<
	ByteImageType, ByteImageType>											ChamferDistanceFilterType;


typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>				GaussianFilterType2;
typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType>				GaussianFilterType3;
// define the fillhole filter
typedef itk::GrayscaleFillholeImageFilter<ImageType, ImageType>				FillholeFilterType;
typedef itk::NeighborhoodIterator<ImageType>								NeighborhoodIteratorType;
typedef itk::NeighborhoodIterator<ImageVectorType>							NeighborhoodVectorIteratorType;


typedef itk::SimpleContourExtractorImageFilter<ImageType, ImageType>		SimpleContourExtractorImageFilterType;

//typedef itk::LocalRoughnessImageFilter<ImageType, ImageType>				LocalRoughnessImageFilterType;

typedef itk::ConnectedComponentImageFilter<ByteImageType, IntImageType>		ConnectedComponentFilterType;
typedef itk::BinaryShapeOpeningImageFilter<ByteImageType>					BinaryShapeOpeningFilterType;
typedef itk::LabelImageToShapeLabelMapFilter< IntImageType, LabelMapType > 	LabelImageToShapeLabelMapFilterType;
typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType >		BinaryMedianFilterType;

//Main Operation Function
ImageType::Pointer RemoveStool(ImageType::Pointer input);
VoxelType SingleMaterialClassification(ImageType::PixelType input_pixel, 
                                       ImageType::PixelType input_gradient_pixel);
ImageType::Pointer AllocateNewImage(ImageType::RegionType fullRegion);
ImageVectorType::Pointer AllocateNewVectorImage(ImageType::RegionType fullRegion);

void VoxelEdgeClassification(float * threshold, VoxelType * previous, double d2, double d1, 
                                      InterpolationType::Pointer input_interpolator, 
                                      IteratorTypeFloat4WithIndex input_smax,
                                      ImageType::IndexType index);
float PolyDist(float X, vnl_real_polynomial poly, float x, float y);
float PolyMinDist(vnl_real_polynomial poly, float x, float y);
float ComputeSmaxFit(float intensity[], float gradient_magnitude[], float Smax);
float ComputeSmax(float intensity[], float gradient_magnitude[], int size);
float Stool_Air_ComputeSmax(float intensity[], float gradient_magnitude[], int size);
float AverageTissueAirDist(float intensity[], float gradient_magnitude[]);
float AverageTissueStoolDist(float Smax, float intensity[], float gradient_magnitude[]);
float AverageStoolAirDist(float Smax, float intensity[], float gradient_magnitude[]);
//ImageType::Pointer RemoveStool(ImageType::Pointer input, int x_offset, int y_offset, int z_offset, int x_size, int y_size, int z_size/*, SSlice*& DAB*/);
vnl_matrix<float> GetNeighbor(ImageVectorType::Pointer partialVector, ImageType::IndexType index);
ByteImageType::Pointer AllocateNewByteImage(ImageType::RegionType fullRegion) ;

vnl_vector<float> expectation(double Y, double mean[], double variance[], float weight[],  vnl_matrix<float> neighbor, float current_partial[], float gradient);
vnl_vector<float> expectation(double Y, float mean[], float variance[], float weight[], float current_partial[]);
void EMClassification(IteratorTypeFloat4WithIndex input_iter, IteratorTypeVoxelType voxel_type_iter, IteratorImageVectorType partialVector_iter,
					  IteratorTypeByteWithIndex chamfer_colon_iter, IteratorTypeFloat4WithIndex temp_iter);
	

//New functions
double round(float d);
void FindVoxelsByGradient(VoxelTypeImage::Pointer voxelEdge, ImageType::IndexType &index, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex, CovariantVectorType &grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector);
void FindVoxelsByGradient2(VoxelTypeImage::Pointer voxelEdge, ImageType::IndexType &index, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex, CovariantVectorType &grad, int numOfVoxels, std::vector<ImageType::IndexType> &indexVector);
void OptimizeVoxelEdge(ImageType::Pointer input, IteratorTypeFloat4WithIndex inputIt, VoxelTypeImage::Pointer &voxelEdge, IteratorTypeVoxelType &voxelEdgeIt );
void WriteITK(ImageType::Pointer image, std::string name);
void WriteITK(ByteImageType::Pointer image, std::string name);
void WriteITK(VoxelTypeImage::Pointer vimage, std::string name);
void WriteITK(IntImageType::Pointer image, std::string name);
void ReadITK(ImageType::Pointer &image, char * fileName);
void ReadITK(VoxelTypeImage::Pointer &vimage, char * fileName);
bool compareSizeOnBorder(LabelObjectType::Pointer a, LabelObjectType::Pointer b);	
bool compareSize(LabelObjectType::Pointer a, LabelObjectType::Pointer b);
int VoxelTypeToNum(VoxelType type);

template <class T> int solveNormalizedCubic (T r, T s, T t, T x[3]);
template <class T> int solveCubic (T a, T b, T c, T d, T x[3]);

typedef unsigned char SSlice[512][512]; 

