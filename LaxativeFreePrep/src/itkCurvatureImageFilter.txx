#ifndef itkCurvatureImageFilter_txx
#define itkCurvatureImageFilter_txx
#include "itkCurvatureImageFilter.h"
namespace itk {
	template<typename InputImage>
	Curvature
	CurvatureImageFilter<InputImage>
	::GetCurvatures(Derivative D) {
		Curvature C;
		float denominator=0.0;    // denominator for the Mean and Gaussian calculations.
	    // calculate the Gaussian curvature
		denominator = (D.x*D.x) + (D.y*D.y) + (D.z*D.z);
		denominator *= denominator;
		if ( denominator == 0.0 ) { C.gaussian = 0.0; }
		else
		{
			C.gaussian =	((D.x*D.x)*((D.yy*D.zz) - (D.yz*D.yz)))
							+ (2*D.y*D.z*((D.xz*D.xy) - (D.xx*D.yz)))
							+ ((D.y*D.y)*((D.xx*D.zz) - (D.xz*D.xz)))
							+ (2*D.x*D.z*((D.yz*D.xy) - (D.yy*D.xz)))
							+ ((D.z*D.z)*((D.xx*D.yy) - (D.xy*D.xy)))
							+ (2*D.x*D.y*((D.xz*D.yz) - (D.zz*D.xy)));
			C.gaussian /= denominator;
		}
	    // calculate the Mean curvature
		denominator = 2*powf(((D.x*D.x) + (D.y*D.y) + (D.z*D.z)), 1.5);
		if ( denominator == 0.0 ) { C.mean = 0.0; }
		else
		{
			denominator  =  1.0/denominator;
			C.mean = ((D.x*D.x)*(D.yy + D.zz)) - (2*D.y*D.z*D.yz)
				     + ((D.y*D.y)*(D.xx + D.zz)) - (2*D.x*D.z*D.xz)
					 + ((D.z*D.z)*(D.xx + D.yy)) - (2*D.x*D.y*D.xy);
			C.mean *= denominator;
		}
	    // calculate the min & max principle curvatures
		float square = (C.mean*C.mean) - C.gaussian;
		if ( square < 0.0 ) square *= -1.0;           // use absolute value
		square = sqrtf(square);
	    C.min = C.mean - square;
		C.max = C.mean + square;
		return C;
	}  // end of get_curvatures


	template<typename InputImage>
	CurvatureImageFilter<InputImage>
	::CurvatureImageFilter()
	{
	}

	template<typename InputImage>
	void 
	CurvatureImageFilter<InputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	
	
	template<typename InputImage>
	void
	CurvatureImageFilter<InputImage>
	::WriteITK(typename InputImageType::Pointer image, char * name) {
		typedef itk::ImageFileWriter< InputImageType >  WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(name);
		writer->SetInput(image);
		std::cerr<<"before update"<<std::endl;
		writer->Update();
		writer.~SmartPointer();
	}
	
	template<typename InputImage>
	void
	CurvatureImageFilter<InputImage>
	::GenerateData()
	{
		typename OutputImageType::Pointer output = GetOutput();
		output->SetRegions(GetInput()->GetRequestedRegion());
		output->Allocate();
		
		InputImageIteratorType input_iter(GetInput(),GetInput()->GetRequestedRegion());
		OutputImageIteratorType output_iter(GetOutput(),GetInput()->GetRequestedRegion());
		
		typename DerivativeImageFilterType::Pointer derivative_filter = DerivativeImageFilterType::New();
		derivative_filter->SetInput(GetInput());
		derivative_filter->SetOrder(1);

		derivative_filter->SetDirection(0);
		//derivative_filter->SetUseImageSpacing(true);
		derivative_filter->Update();

		typename InputImageType::Pointer x_partial;
		x_partial = derivative_filter->GetOutput();
		TempImageIteratorType x_iter(x_partial,x_partial->GetRequestedRegion());


		derivative_filter.~SmartPointer();
		derivative_filter = DerivativeImageFilterType::New();
		derivative_filter->SetInput(GetInput());
		derivative_filter->SetOrder(1);
		derivative_filter->SetDirection(1);
		//derivative_filter->SetUseImageSpacing(true);
		derivative_filter->Update();

		typename InputImageType::Pointer y_partial;
		y_partial= derivative_filter->GetOutput();
		TempImageIteratorType y_iter(y_partial,y_partial->GetRequestedRegion());
		

		derivative_filter.~SmartPointer();
		derivative_filter = DerivativeImageFilterType::New();
		derivative_filter->SetInput(GetInput());
		derivative_filter->SetOrder(1);
		derivative_filter->SetDirection(2);
		//derivative_filter->SetUseImageSpacing(true);
		derivative_filter->Update();

		typename InputImageType::Pointer z_partial;
		z_partial= derivative_filter->GetOutput();
		TempImageIteratorType z_iter(z_partial,z_partial->GetRequestedRegion());
		//WriteITK(x_partial, "x.mhd");
		//WriteITK(y_partial, "y.mhd");
		//WriteITK(z_partial, "z.mhd");

		derivative_filter.~SmartPointer();

		typename HessianImageFilterType::Pointer hessian_image_filter = HessianImageFilterType::New();
		hessian_image_filter->SetInput(GetInput());
		hessian_image_filter->SetSigma(0);
		hessian_image_filter->Update();

		typename HessianImageType::Pointer hessian_image = hessian_image_filter->GetOutput();
		HessianImageIteratorType hessian_image_iter(hessian_image, hessian_image->GetRequestedRegion());
		
		for (input_iter.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(), hessian_image_iter.GoToBegin(), output_iter.GoToBegin(); 
			!input_iter.IsAtEnd(); 
			++input_iter, ++x_iter, ++y_iter, ++z_iter, ++hessian_image_iter, ++output_iter)
		{
			Derivative temp_data;
			temp_data.x=x_iter.Get();
			temp_data.y=y_iter.Get();
			temp_data.z=z_iter.Get();

			HessianValueType hessian_value = hessian_image_iter.Get();
			temp_data.xx=hessian_value(0,0);
			temp_data.xy=hessian_value(0,1);
			temp_data.xz=hessian_value(0,2);

			temp_data.yy=hessian_value(1,1);
			temp_data.yz=hessian_value(1,2);

			temp_data.zz=hessian_value(2,2);

			output_iter.Set(GetCurvatures(temp_data));
		
		}
		

		//InputImageNeighborhoodIteratorType::SizeType neighborhood_size={10,10,10};
		//InputImageNeighborhoodIteratorType input_neighborhood_iter(neighborhood_size,GetInput(),GetInput()->GetRequestedRegion());
		//OutputImageNeighborhoodIteratorType output_neighborhood_iter(neighborhood_size,output,output->GetRequestedRegion());

		//for (input_neighborhood_iter.GoToBegin(), output_neighborhood_iter.GoToBegin();
		//	!input_neighborhood_iter.IsAtEnd() && !output_neighborhood_iter.IsAtEnd(); 
		//	++input_neighborhood_iter, ++output_neighborhood_iter)
		//{
		//	if (input_neighborhood_iter.IsInBound()) {
		//		InputImageType::
		//		input_neighborhood_iter.GetBoundingBoxAsImageRegion()
		//	}
		//}
		//HessianImageFilterType::Pointer hessian_image_filter = HessianImageFilterType::New();
		//hessian_image_filter->SetSigma(.5);
		//hessian_image_filter->SetInput();
		//hessian_image_filter->Update();

	}
}
#endif // itkCurvatureImageFilter_txx
