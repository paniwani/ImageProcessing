#ifndef itkCurvatureImageFilter_h
#define itkCurvatureImageFilter_h
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

#include "utils.h"
namespace itk
{
	template<typename InputImage>
	class ITK_EXPORT CurvatureImageFilter : 
		public ImageToImageFilter<InputImage, Image<Curvature,itk::GetImageDimension<InputImage>::ImageDimension> >
	{
		public :
		    typedef CurvatureImageFilter Self;
			typedef ImageToImageFilter<InputImage, Image<Curvature,itk::GetImageDimension<InputImage>::ImageDimension> > Superclass;
			typedef SmartPointer<Self> Pointer;
			typedef SmartPointer<Self const> ConstPointer;
			itkNewMacro(Self);	
			itkTypeMacro(CurvatureImageFilter, ImageToImageFilter);

			/** Typedef support of method types. */
			typedef InputImage InputImageType;
			typedef typename InputImageType::Pointer InputImagePointer;
			typedef typename InputImageType::ConstPointer InputImageConstPointer;
			typedef typename InputImageType::RegionType RegionType; 
			typedef typename InputImageType::PixelType PixelType; 
			typedef typename InputImageType::IndexType IndexType;
			typedef typename InputImageType::SizeType  SizeType;

			typedef Image<Curvature,itk::GetImageDimension<InputImage>::ImageDimension> OutputImageType;

			CurvatureImageFilter();
		protected :
			typedef itk::ConstNeighborhoodIterator<InputImageType> InputImageNeighborhoodIteratorType;
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

			Curvature GetCurvatures(Derivative D);
		private :
			
			void WriteITK(typename InputImageType::Pointer image, char * name);

			
		    CurvatureImageFilter(Self const &); // not implemented
		    Self & operator=(Self const &); // not implemented
	};
}



#ifndef ITK_MANUAL_INSTANTIATION

#include "itkCurvatureImageFilter.txx"

#endif



#endif // itkCurvatureImageFilter_h
