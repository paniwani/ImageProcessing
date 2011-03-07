#ifndef itkLocalRoughnessImageFilter_h
#define itkLocalRoughnessImageFilter_h
#include <itkImageToImageFilter.h>
#include <itkImage.h>

#include <itkNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkDerivativeOperator.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkConstNeighborhoodIterator.h>

#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkDerivativeImageFilter.h>

#include <itkImageFileWriter.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix.h>
#include <math.h>
#include "utils.h"
namespace itk
{
	template<typename InputImage, typename OutputImage>
	class ITK_EXPORT LocalRoughnessImageFilter : 
		public ImageToImageFilter<InputImage, OutputImage>
	{
		public :
		    typedef LocalRoughnessImageFilter Self;
			typedef ImageToImageFilter<InputImage, OutputImage > Superclass;
			typedef SmartPointer<Self> Pointer;
			typedef SmartPointer<Self const> ConstPointer;
			itkNewMacro(Self);	
			itkTypeMacro(LocalRoughnessImageFilter, ImageToImageFilter);

			/** Typedef support of method types. */
			typedef InputImage InputImageType;
			typedef typename InputImageType::Pointer InputImagePointer;
			typedef typename InputImageType::ConstPointer InputImageConstPointer;
			typedef typename InputImageType::RegionType RegionType; 
			typedef typename InputImageType::PixelType PixelType; 
			typedef typename InputImageType::IndexType IndexType;
			typedef typename InputImageType::SizeType  SizeType;

			typedef OutputImage OutputImageType;

			LocalRoughnessImageFilter();
		protected :
			
			


			
			typedef itk::ConstNeighborhoodIterator<InputImage> ConstNeighborhoodIteratorType;
			typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;  









			

			typedef itk::NeighborhoodIterator<OutputImageType> OutputImageNeighborhoodIteratorType;

			typedef itk::NeighborhoodInnerProduct<InputImageType> IP;
			typedef itk::DerivativeOperator<PixelType, itk::GetImageDimension<InputImage>::ImageDimension> DerivativeOperatorType;	

			typedef itk::SymmetricSecondRankTensor< PixelType, itk::GetImageDimension<InputImage>::ImageDimension > HessianValueType;

			typedef itk::Image<HessianValueType,InputImageType::ImageDimension> HessianImageType;

			typedef itk::HessianRecursiveGaussianImageFilter<InputImageType,HessianImageType > HessianImageFilterType;

			typedef itk::ImageRegionConstIteratorWithIndex<InputImageType> InputImageIteratorType;

			typedef itk::ImageRegionIteratorWithIndex<InputImageType> TempImageIteratorType;
			typedef itk::ImageRegionIteratorWithIndex<OutputImageType> OutputImageIteratorType;
			typedef itk::ImageRegionIteratorWithIndex<HessianImageType> HessianImageIteratorType;

			typedef itk::DerivativeImageFilter<InputImageType,InputImageType> DerivativeImageFilterType;

		    void PrintSelf(std::ostream& os, Indent indent) const;
		    void GenerateData();    

			

			Curvature GetCurvatures(Derivative  &D);
			
			typename InputImageType::Pointer AllocateNewImage(typename InputImageType::RegionType region) {
				typename InputImageType::Pointer newImage = InputImageType::New();
				newImage->SetLargestPossibleRegion( region );
				newImage->SetBufferedRegion( region );
				newImage->SetRequestedRegion( region );
				newImage->Allocate();
				return newImage;
			}
			float WeightingFunction(float input);

			//Functions for the Deirche Filter

			float ComputeF0(float index, float alpha, vnl_vector<float> constant);
			float ComputeF1(float index, float alpha, vnl_vector<float> constant);
			float ComputeF2(float index, float alpha, vnl_vector<float> constant);

			typename InputImage::Pointer GetSubImage(int filter_width, int filter_height, int filter_depth, typename InputImage::IndexType index);

			float DoKernel(typename InputImage::Pointer filter, typename InputImage::Pointer sub_image);
	
			vnl_vector<float> ComputeConstants(float alpha, float step, float range);
			void ComputeFilter(float alpha1, float alpha2, int filter_width, int filter_height, int filter_depth, 
				typename InputImage::Pointer FilterX, typename InputImage::Pointer FilterY, typename InputImage::Pointer FilterZ,
				typename InputImage::Pointer FilterXX, typename InputImage::Pointer FilterXY, typename InputImage::Pointer FilterXZ,
				typename InputImage::Pointer FilterYY, typename InputImage::Pointer FilterYZ, typename InputImage::Pointer FilterZZ);

			/** Directly Set/Get the array of weights used in the gradient calculations.
			  Note that calling UseImageSpacingOn will clobber these values.*/
			void SetDerivativeWeights(float data[]);
			itkGetVectorMacro(DerivativeWeights, const float, itk::GetImageDimension<TInputImage>::ImageDimension);
			
			vnl_matrix_fixed< float, 3, 3> EvaluateAtNeighborhood (const vnl_matrix_fixed<float ,3,3> Next, const vnl_matrix_fixed<float ,3,3> Previous) const{
			
				unsigned i, j;
				vnl_matrix_fixed<float ,3,3> J;

				for (i = 0; i < 3; ++i) {
					for (j = 0; j < 3; ++j) {
						J[i][j] = m_DerivativeWeights[i]
							* 0.5 * (Next[i][j] - Previous[i][j]);
					}
				}

				return J;
			}

			/** The weights used to scale partial derivatives during processing */
			float m_DerivativeWeights[3];
		private :
			void WriteITK(typename InputImageType::Pointer image, char * name);
		    LocalRoughnessImageFilter(Self const &); // not implemented
		    Self & operator=(Self const &); // not implemented



	};
}


#ifndef ITK_MANUAL_INSTANTIATION

#include "itkLocalRoughnessImageFilter.txx"

#endif



#endif // itkLocalRoughnessImageFilter