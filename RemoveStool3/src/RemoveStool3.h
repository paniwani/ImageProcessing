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
//#include <itkBinaryShapeOpeningImageFilter.h>
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkOrientImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <itkHessian3DToVesselnessMeasureImageFilter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>


//#include <algorithm>
//#include <vector>

// Outdated includes
//#include "Typedef.h"
//#include "utils.h"
//#include "itkLocalRoughnessImageFilter.h"
//#include <NIH_Algo_TextureAnalysis.h>
//#include <DynArray.h>

// Regex
//#include <boost/xpressive/xpressive.hpp>
//#include "pcre.h"

#include <itkCastImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkFastMarchingImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkSimilarityIndexImageFilter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkDerivativeImageFilter.h>
#include "itkSimpleFuzzyConnectednessScalarImageFilter.h"
#include <itkRegionalMinimaImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkLabelShapeOpeningImageFilter.h>

#include <itkAnisotropicEdgeEnhancementDiffusionImageFilter.h>
#include <itkDiscreteHessianGaussianImageFunction.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkShapeRelabelImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkSliceBySliceImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkAnisotropicHybridDiffusionImageFilter.h>
#include <itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkMucosalReconstructionFilter.h>


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

typedef  itk::FixedArray< double, 3 >										EigenValueArrayType;

typedef unsigned long														LabelType;
typedef itk::ShapeLabelObject< LabelType, 3 >								LabelObjectType;
typedef itk::LabelMap< LabelObjectType >									LabelMapType;

//CovariantVector Type
typedef itk::CovariantVector<PixelType,3>                                   CovariantVectorType;
typedef itk::Image<PixelType, 3>                                            ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType>                        IteratorTypeFloat4WithIndex;
typedef itk::Image<CovariantVectorType, 3>                                  ImageVectorType;
typedef itk::Image< EigenValueArrayType, 3 >								EigenValueImageType;


typedef itk::Image<BYTE, 3>													ByteImageType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType>                    IteratorTypeByteWithIndex;

typedef itk::Image<int, 3>													IntImageType;
typedef itk::ImageRegionIteratorWithIndex<IntImageType>						IteratorTypeIntWithIndex;

typedef itk::Image<unsigned int, 3>											UintImageType;
typedef itk::ImageRegionIteratorWithIndex<UintImageType>					IteratorTypeUintWithIndex;

//typedef itk::Image<unsigned short,3>										UshortImageType;
//typedef itk::ImageRegionIteratorWithIndex<UshortImageType>					UshortImageTypeWithIndex;

typedef itk::ImageRegionIteratorWithIndex<EigenValueImageType>				EigenValueIterType;


typedef itk::ImageRegionIterator<ImageVectorType>                           IteratorImageVectorType;
typedef itk::Image<VoxelType, 3>                                            VoxelTypeImage;
typedef itk::ImageRegionIteratorWithIndex<VoxelTypeImage>                   IteratorTypeVoxelType;
//Axillary Types
typedef itk::BinaryBallStructuringElement<PixelType, 3>                     StructuringElementType;
typedef itk::BSplineInterpolateImageFunction<ImageType>						InterpolationType;
typedef itk::BSplineInterpolateImageFunction<ImageVectorType, float>		InterpolationVectorType;


typedef itk::GaussianBlurImageFunction<ImageVectorType>                     GaussianBlurVectorType;

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


typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>				DiscreteGaussianFilterType;
//typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType>				GaussianFilterType3;
// define the fillhole filter
typedef itk::GrayscaleFillholeImageFilter<ImageType, ImageType>				FillholeFilterType;
typedef itk::NeighborhoodIterator<ImageType>								NeighborhoodIteratorType;
typedef itk::NeighborhoodIterator<VoxelTypeImage>							NeighborhoodIteratorVoxelType;
typedef itk::NeighborhoodIterator<ImageVectorType>							NeighborhoodVectorIteratorType;


typedef itk::SimpleContourExtractorImageFilter<ImageType, ImageType>		SimpleContourExtractorImageFilterType;

//typedef itk::LocalRoughnessImageFilter<ImageType, ImageType>				LocalRoughnessImageFilterType;

typedef itk::ConnectedComponentImageFilter<ByteImageType, IntImageType>		ConnectedComponentFilterType;
//typedef itk::BinaryShapeOpeningImageFilter<ByteImageType>					BinaryShapeOpeningFilterType;
typedef itk::LabelImageToShapeLabelMapFilter< IntImageType, LabelMapType >  LabelImageToShapeLabelMapFilterType;
//typedef itk::LabelImageToLabelMapFilter<IntImageType>						LabelImageToLabelMapFilterType;
typedef itk::BinaryMedianImageFilter< ByteImageType, ByteImageType >		BinaryMedianFilterType;
typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType >	    HoleFillingFilterType;
typedef itk::HessianRecursiveGaussianImageFilter<ImageType>					HessianGaussianFilterType;
typedef itk::SymmetricEigenAnalysisImageFilter< HessianGaussianFilterType::OutputImageType, EigenValueImageType> EigenAnalysisFilterType;
typedef itk::RescaleIntensityImageFilter< ImageType >						RescaleIntensityFilterType;
typedef itk::ResampleImageFilter< ImageType, ImageType >					ResampleImageFilterType;
typedef itk::DerivativeImageFilter<ImageType,ImageType>						DerivativeImageFilterType;
typedef itk::Hessian3DToVesselnessMeasureImageFilter<float>					HessianVesselnessFilterType;
typedef itk::HessianToObjectnessMeasureImageFilter< HessianGaussianFilterType::OutputImageType, ImageType> ObjectnessFilterType;

typedef itk::SimpleFuzzyConnectednessScalarImageFilter<ImageType,ImageType> FuzzyFilterType;
typedef FuzzyFilterType::FuzzySceneType										FuzzySceneType;

typedef itk::SymmetricSecondRankTensor<float>								SymmetricSecondRankTensorType;
typedef itk::Image<SymmetricSecondRankTensorType, 3>						HessianImageType;			

//Main Operation Function
ImageType::Pointer RemoveStool(ImageType::Pointer input);
VoxelType SingleMaterialClassification(ImageType::PixelType input_pixel, 
                                       ImageType::PixelType input_gradient_pixel);
ImageType::Pointer AllocateNewImage(ImageType::RegionType fullRegion);
ImageVectorType::Pointer AllocateNewVectorImage(ImageType::RegionType fullRegion);



void VoxelEdgeClassification(float * threshold, VoxelType * previous, double d2, double d1,
                                      InterpolationType::Pointer &input_interpolator, 
									  InterpolationType::Pointer &gradient_magnitude_interpolator,
                                      IteratorTypeFloat4WithIndex &input_smax,
                                      ImageType::IndexType &index,
									  CovariantVectorType &gradient);

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

vnl_vector<float> expectation(double Y, double mean[], double variance[], float weight[],  vnl_matrix<float> neighbor, float current_partial[]);
vnl_vector<float> expectation(double Y, float mean[], float variance[], float weight[], float current_partial[]);
void EMClassification(IteratorTypeFloat4WithIndex input_iter, IteratorTypeVoxelType voxel_type_iter, IteratorImageVectorType partialVector_iter,
					  IteratorTypeByteWithIndex chamfer_colon_iter, IteratorTypeFloat4WithIndex temp_iter);
	

//New functions
double round(float d);
void WriteITK(ImageType::Pointer image, std::string name);
void WriteITK(FuzzyFilterType::FuzzySceneType, std::string name);
void WriteITK(ByteImageType::Pointer image, std::string name);
void WriteITK(VoxelTypeImage::Pointer vimage, std::string name);
void WriteITK(IntImageType::Pointer image, std::string name);
void WriteITK(UintImageType::Pointer image, std::string name);
void ReadITK(ImageType::Pointer &image, char * fileName);
void ReadITK(VoxelTypeImage::Pointer &vimage, char * fileName);
void ReadITK(ByteImageType::Pointer &image, char * fileName);
bool compareSizeOnBorder(LabelObjectType::Pointer a, LabelObjectType::Pointer b);	
bool compareSize(LabelObjectType::Pointer a, LabelObjectType::Pointer b);
int VoxelTypeToNum(VoxelType type);
//VoxelType NumToVoxelType(int num);
ImageType::Pointer ReadDicom( std::string path );
void SmoothPartialVector(ImageVectorType::Pointer pv, ByteImageType::Pointer chamfer_colon, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex);

void RunFuzzy(ImageType::Pointer input, ByteImageType::Pointer chamfer_colon, VoxelTypeImage::Pointer voxel_type);
void RunConnectedThresholdGrowing(ImageType::Pointer input, ByteImageType::Pointer chamfer_colon, VoxelTypeImage::Pointer voxel_type);

float vesselness(const float lambda1, const float lambda2, const float lambda3, const float alpha, const float gamma12, const float gamma23);
float thinness( float l[3] );

bool OrderByMagnitude (double a, double b);
bool OrderByValue (double a, double b);
bool OrderByValueDesc (double a, double b);
ByteImageType::Pointer FindBinaryEdge(ByteImageType::Pointer im, ImageType::IndexType &startIndex, ImageType::IndexType &endIndex);

double fA(EigenValueArrayType lambda, double alpha);
double fB(EigenValueArrayType lambda, double beta, double gamma);
double fC(double ev1, double ev2, double eta);
double fRut(EigenValueArrayType lambda, double alpha, double beta, double gamma);
double fCup(EigenValueArrayType lambda, double eta);
ImageType::Pointer ResampleImage(ImageType::Pointer input, ImageType::SpacingType output_spacing);
bool MatchVoxels(std::string s);
void VoteVoxels(VoxelTypeImage::Pointer v, ByteImageType::Pointer mask);
ImageType::Pointer ComputeNeighborhoodSmax(ImageType::Pointer input, VoxelTypeImage::Pointer v, IteratorTypeByteWithIndex &mask_iter, ImageType::IndexType &sidx, ImageType::IndexType &eidx);
void WritePartialImages(ImageVectorType::Pointer partialVector, ByteImageType::Pointer chamfer_colon, std::string name);
vnl_matrix_fixed< float, 3, 3> EvaluateAtNeighborhood (const vnl_matrix_fixed<float,3,3> Next, const vnl_matrix_fixed<float,3,3> Previous, float weights[3]);

template <class T> int solveNormalizedCubic (T r, T s, T t, T x[3]);
template <class T> int solveCubic (T a, T b, T c, T d, T x[3]);

typedef unsigned char SSlice[512][512]; 

void MeasureObjectness(ImageType::Pointer input);
void ComputeSatoHessian(ImageType::Pointer input);
void ComputeFrangiHessian(ImageType::Pointer input);

float omega(float ls, float lt, float gamma);
float weight(float ls, float lt, float alpha, float gamma);
float omega_dark(float ls, float lt, float gamma);
float weight_dark(float ls, float lt, float alpha, float gamma);

float S_sheet(float lambda[3], float alpha, float gamma);
float S_sheet_dark(float lambda[3], float alpha, float gamma);
float S_blob(float lambda[3], float gamma);
float S_line(float lambda[3], float alpha, float gamma);
float S_line_dark(float lambda[3], float alpha, float gamma);
float S_fold(float lambda[3], float alpha, float gamma);

void HeuristicClosing(VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon);
std::vector<std::string> explode( const std::string &delimiter, const std::string &str);

std::string VoxelTypeToString(VoxelType type);
ImageType::Pointer HessianResponse(ImageType::Pointer input_aniso, double alpha, double beta, double gamma, double eta);
void EnhanceVoxelType(ImageType::Pointer input, ImageType::Pointer hessian, VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon);
ImageType::Pointer SatoResponse(ImageType::Pointer input_aniso, ByteImageType::Pointer chamfer_colon, double alpha, double gamma);
//void ProcessHessian(ImageType::Pointer hessian, VoxelTypeImage::Pointer voxel_type);

//ImageType::Pointer SatoHessianEdgeEnhancingDiffusion(ImageType::Pointer input);
void GetHessian(ImageType::Pointer input_aniso);
void UnderstandHessian(ImageType::Pointer input_aniso);

ByteImageType::Pointer SegmentColon(ImageType::Pointer input);
void CleanIsolatedStool(VoxelTypeImage::Pointer voxel_type);
ImageType::Pointer SatoModifiedResponse(ImageType::Pointer input_aniso, ByteImageType::Pointer chamfer_colon, double sigma, double gamma);
ImageType::Pointer Sharpen(ImageType::Pointer input);
void DiffusionTest( ImageType::Pointer input, double ContrastParameter, double sigma);
void SubtractStool( ImageType::Pointer input, VoxelTypeImage::Pointer voxel_type, ByteImageType::Pointer chamfer_colon);

/// PARAMS
float Modified=false;
double PI=3.1415926;
double neighbor_weight[3]={1,1,.5};
double beta=.7;
double weight_sum=2.5;
int chamfer_weights[3]={3,4,5};

bool writeNum = true;
int writeCount=1;

std::string note = "";
const float N_SVAL = 400;

const double ALPHA = 0.5;
const double BETA = 0.3;
const double GAMMA = 0.3;
const double ETA = 0.2;

bool truncateOn = true;
unsigned int truncate_ar[2] = {85,100};