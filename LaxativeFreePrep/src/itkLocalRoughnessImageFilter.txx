#ifndef itkLocalRoughnessImageFilter_txx
#define itkLocalRoughnessImageFilter_txx
#include "itkLocalRoughnessImageFilter.h"

namespace itk {

	template<typename InputImage, typename OutputImage>
	float
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::WeightingFunction(float value)
	{
		float sigma=2;
		return 100*exp(-vnl_math_sqr(value)/(2*sigma))/(sqrt(2*PI*sigma));
	}

	template<typename InputImage, typename OutputImage>
	float
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::ComputeF1(float index, float alpha, vnl_vector<float> constant)
	{
		return constant[1]*index*vnl_math_sqr(alpha)*exp(-alpha*abs(index));
	}
	
	template<typename InputImage, typename OutputImage>
	float
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::ComputeF2(float index, float alpha, vnl_vector<float> constant)
	{
		return constant[2]*(1 - constant[3]*alpha*abs(index))*exp(-alpha*abs(index));
	}

	template<typename InputImage, typename OutputImage>
	float
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::ComputeF0(float index, float alpha, vnl_vector<float> constant)
	{
		return constant[0]*(1+alpha*abs(index))*exp(-alpha*abs(index));
	}

	template<typename InputImage, typename OutputImage>
	typename InputImage::Pointer
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::GetSubImage(int filter_width, int filter_height, int filter_depth, typename InputImage::IndexType index){
		typename InputImage::RegionType imageRegion = GetInput()->GetRequestedRegion();

		int min_x=imageRegion.GetIndex()[0];
		int max_x=min_x+imageRegion.GetSize()[0]-1;

		int min_y=imageRegion.GetIndex()[1];
		int max_y=min_y+imageRegion.GetSize()[1]-1;

		int min_z=imageRegion.GetIndex()[2];
		int max_z=min_z+imageRegion.GetSize()[2]-1;

		typename InputImage::IndexType filter_Index={-filter_width/2,-filter_height/2,-filter_depth/2};
		typename InputImage::SizeType filter_Size={filter_width,filter_height,filter_depth};
		typename InputImage::RegionType filter_region(filter_Index,filter_Size);
		typename InputImage::Pointer temp = AllocateNewImage(filter_region);


		TempImageIteratorType temp_iter(temp, filter_region);

		for(temp_iter.GoToBegin(); !temp_iter.IsAtEnd(); ++temp_iter) {
			typename InputImage::IndexType filter_index=temp_iter.GetIndex();
			int x= std::min<int>(std::max<int>(min_x, index[0]+filter_index[0]),max_x);
			int y= std::min<int>(std::max<int>(min_y, index[1]+filter_index[1]),max_y);
			int z= std::min<int>(std::max<int>(min_z, index[2]+filter_index[2]),max_z);
			//std::cerr<<x<<" "<<y<<" "<<z<<std::endl;
			typename InputImage::IndexType temp_index={x,y,z};
			temp_iter.Set(GetInput()->GetPixel(temp_index));
		}
		return temp;
	}

	template<typename InputImage, typename OutputImage>
	float
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::DoKernel(typename InputImage::Pointer filter, typename InputImage::Pointer sub_image)
	{
		typename InputImage::RegionType imageRegion = filter->GetRequestedRegion();
		TempImageIteratorType filter_iter(filter, imageRegion);
		TempImageIteratorType sub_image_iter(sub_image, imageRegion);

		float totalSum=0;

		for(filter_iter.GoToBegin(), sub_image_iter.GoToBegin(); !filter_iter.IsAtEnd() && !sub_image_iter.IsAtEnd(); ++filter_iter, ++sub_image_iter) {
			totalSum+=filter_iter.Get()*sub_image_iter.Get();
		}

		return totalSum;
	}

	/*---------------------------------------------------
	 *  ComputeConstants -- This function computes the constatns for the Deriche Filter
	 *
	 */
	template<typename InputImage, typename OutputImage>
	vnl_vector<float>
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::ComputeConstants(float alpha, float step, float range)
	{
		float integral0  = 0.0;
		float integral1  = 0.0;
		float integral2a = 0.0;
		float integral2b = 0.0;
		float integral3a = 0.0;
		float integral3b = 0.0;

		float start=(-1*step*(range-1)/2);
		float end=step*((range-1)/2);
		for (float x = start; x <= end; x++) {
			
			float e_val = exp(-alpha*abs(x));
			
			float f0 = (1 + alpha*abs(x))*e_val;
			float f1 = x*vnl_math_sqr(alpha)*e_val;
			float f2a = e_val;
			float f2b = alpha*abs(x)*e_val;

			integral0 = integral0 + f0;
			integral1  = integral1 + x*f1;
			integral2a = integral2a + f2a;
			integral2b = integral2b + f2b;
			integral3a = integral3a + vnl_math_sqr(x)*0.5*f2a;
			integral3b = integral3b + vnl_math_sqr(x)*0.5*f2b;

		}
		float return_value[] = {1/integral0, 1/integral1, 1/((integral3a - (integral2a/integral2b*integral3b))), integral2a/integral2b};
		return vnl_vector<float>(4,4, return_value);
	}

	/*----------------------------------------------------
	 *  ComputeFilter -- Computes the filter for the partial derivative
	 *                   calculation by convolution.
	 *
	 *  Inputs:
	 *      alpha1, alpha2 -- parameters for smoothing in the filter
	 *      filter_width, filter_height
	 */

	template<typename InputImage, typename OutputImage>
	void
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::ComputeFilter(float alpha1, float alpha2, int filter_width, int filter_height, int filter_depth, 
		typename InputImage::Pointer FilterX, typename InputImage::Pointer FilterY, typename InputImage::Pointer FilterZ,
		typename InputImage::Pointer FilterXX, typename InputImage::Pointer FilterXY, typename InputImage::Pointer FilterXZ,
		typename InputImage::Pointer FilterYY, typename InputImage::Pointer FilterYZ, typename InputImage::Pointer FilterZZ)
	{

		float del_x=1;
		float del_y=1;
		float del_z=1;

		//% get constants for first order derivatives

		vnl_vector<float> constants_first_width = ComputeConstants(alpha1, del_x, filter_width);
		vnl_vector<float> constants_first_height = ComputeConstants(alpha1, del_y, filter_height);
		vnl_vector<float> constants_first_depth = ComputeConstants(alpha1, del_z, filter_depth);

		//% get constants for second order derivatives

		vnl_vector<float> constants_second_width = ComputeConstants(alpha2, del_x, filter_width);
		vnl_vector<float> constants_second_height = ComputeConstants(alpha2, del_y, filter_height);
		vnl_vector<float> constants_second_depth = ComputeConstants(alpha2, del_z, filter_depth);

		//% calculate filter values through loops

		//counter = 1;
		//z = -1*del_z*(n_depth-1)/2;

		typename InputImage::IndexType index = {-filter_width/2, -filter_height/2, -filter_depth/2};
		typename InputImage::SizeType size = {filter_width, filter_height, filter_depth};
		
		typename InputImage::RegionType iteration_region(index, size);
	
		TempImageIteratorType filter_x = TempImageIteratorType(FilterX, iteration_region);
		TempImageIteratorType filter_y = TempImageIteratorType(FilterY, iteration_region);
		TempImageIteratorType filter_z = TempImageIteratorType(FilterZ, iteration_region);

		TempImageIteratorType filter_xx = TempImageIteratorType(FilterXX, iteration_region);
		TempImageIteratorType filter_xy = TempImageIteratorType(FilterXY, iteration_region);
		TempImageIteratorType filter_xz = TempImageIteratorType(FilterXZ, iteration_region);

		TempImageIteratorType filter_yy = TempImageIteratorType(FilterYY, iteration_region);
		TempImageIteratorType filter_yz = TempImageIteratorType(FilterYZ, iteration_region);
		TempImageIteratorType filter_zz = TempImageIteratorType(FilterZZ, iteration_region);

		for (filter_x.GoToBegin(), filter_y.GoToBegin(), filter_z.GoToBegin(), filter_xx.GoToBegin(), filter_xy.GoToBegin(),
			filter_xz.GoToBegin(), filter_yy.GoToBegin(), filter_yz.GoToBegin(), filter_zz.GoToBegin();
			!filter_x.IsAtEnd() && !filter_y.IsAtEnd()&& ! filter_z.IsAtEnd()&& ! filter_xx.IsAtEnd()&& ! filter_xy.IsAtEnd()&& 
			!filter_xz.IsAtEnd()&& ! filter_yy.IsAtEnd()&& ! filter_yz.IsAtEnd()&& ! filter_zz.IsAtEnd();
			++filter_x, ++filter_y, ++filter_z, ++filter_xx, ++filter_xy, ++filter_xz, ++filter_yy, ++filter_yz, ++filter_zz)
		{
			//Can save computation time if stored differently, reduce cubic to linear.
			typename InputImage::IndexType temp_index = filter_x.GetIndex();
			int x = temp_index[0];
			int y = temp_index[1];
			int z = temp_index[2];
			
			filter_x.Set(ComputeF1(x,alpha1, constants_first_width)*ComputeF0(y, alpha1, constants_first_height)*ComputeF0(z,alpha1, constants_first_depth));
			filter_y.Set(ComputeF0(x,alpha1, constants_first_width)*ComputeF1(y, alpha1, constants_first_height)*ComputeF0(z,alpha1, constants_first_depth));
			filter_z.Set(ComputeF0(x,alpha1, constants_first_width)*ComputeF0(y, alpha1, constants_first_height)*ComputeF1(z,alpha1, constants_first_depth));

			filter_xx.Set(ComputeF2(x,alpha2, constants_second_width)*ComputeF0(y, alpha1, constants_first_height)*ComputeF0(z,alpha1, constants_first_depth));
			filter_yy.Set(ComputeF0(x,alpha1, constants_first_width)*ComputeF2(y, alpha2, constants_second_height)*ComputeF0(z,alpha1, constants_first_depth));
			filter_zz.Set(ComputeF0(x,alpha1, constants_first_width)*ComputeF0(y, alpha1, constants_first_height)*ComputeF2(z,alpha2, constants_second_depth));

			filter_xy.Set(ComputeF1(x,alpha1, constants_first_width)*ComputeF1(y, alpha1, constants_first_height)*ComputeF0(z,alpha1, constants_first_depth));
			filter_yz.Set(ComputeF0(x,alpha1, constants_first_width)*ComputeF1(y, alpha1, constants_first_height)*ComputeF1(z,alpha1, constants_first_depth));
			filter_xz.Set(ComputeF1(x,alpha1, constants_first_width)*ComputeF0(y, alpha1, constants_first_height)*ComputeF1(z,alpha1, constants_first_depth));
		}
	}
	
	
	/*-----------------------------------------------------
	 *  get_curvatures --  function which calculates the curvatures
	 *             of the surface at the vertex, using the
	 *             derivatives calculated at the vertex.
	 *
	 *  Inputs:
	 *      Derivative   &D          -  the partial derivatives
	 *
	 *  Outputs:
	 *      Curvature    &C          -  the curvatures calculated from
	 *                                  the derivatives
	 *
	 *  Description:
	 *
	 *  The Gaussian curvature, and the Mean curvature of the isosurface
	 *  are calculated using the calculated derivatives for a point.
	 *  The equations used to calculate the curvatures are from the 
	 *  paper "Computing the Differential Characteristics of Isointensity 
	 *  Surfaces" by Jean-Philippe Thirion and Alexis Gourdon, published
	 *  by Computer Vision and Image Understanding, Vol. 61, No. 2,
	 *  March 1995.
	 *
	 *  These are the formulas for the Gaussian and Mean curvatures:
	 *  
	 *                       _
	 *      Gaussian =   1  |
	 *                  --- | (D.x)^2(D.yy*D.zz - (D.yz)^2) + 2D.y*D.z(D.xz*D.xy - D.xx*D.yz)
	 *                  h^2 |_ 
	 *
	 *                      + (D.y)^2(D.xx*D.zz - (D.xz)^2) + 2D.x*D.z(D.yz*D.xy - D.yy*D.xz)
	 *                                                                                       _
	 *                                                                                        |
	 *                      + (D.z)^2(D.xx*D.yy - (D.xy)^2) + 2D.x*D.y(D.xz*D.yz - D.zz*D.xy) |
	 *                                                                                       _|
	 *                       _
	 *      Mean  =   1     | 
	 *              ------- | (D.x)^2(D.yy + D.zz) -2D.y*D.z*D.yz + (D.y)^2(D.xx + D.zz)
	 *              2(h^1.5)|_                                                        _
	 *                                                                                 |
	 *                         - 2D.x*D.z*D.xz + (D.z)^2(D.xx + D.yy) - 2D.x*D.y*D.xy  |
	 *                                                                                _|
	 *
	 *  where h = (D.x)^2 + (D.y)^2 + (D.z)^2
	 *
	 *  The principle curvatures (the min and the max) are found using the calculated
	 *  Gaussian and Mean curvatures:
	 *
	 *
	 *       Min  = Mean - sqrt(Mean^2 - Gaussian)
	 *
	 *       Max  = Mean + sqrt(Mean^2 - Gaussian)
	 *
	 *-----------------------------------------------------*/
	template<typename InputImage, typename OutputImage>
	Curvature
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::GetCurvatures(Derivative &D) {
		//std::cerr<<D.x<<" "<<D.y<<" "<<D.z<<std::endl;
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
		//std::cerr<<C.min<<" "<<C.max<<std::endl;
		return C;
	}  // end of get_curvatures



	template<typename InputImage, typename OutputImage>
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::LocalRoughnessImageFilter()
	{
	}

	template<typename InputImage, typename OutputImage>
	void
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	template<typename InputImage, typename OutputImage>
	void
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::WriteITK(typename InputImageType::Pointer image, char * name) {
		typedef itk::ImageFileWriter< InputImageType >  WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(name);
		writer->SetInput(image);
		std::cerr<<"before update"<<std::endl;
		writer->Update();
		writer.~SmartPointer();
	}
	
	template<typename InputImage, typename OutputImage>
	void
	LocalRoughnessImageFilter<InputImage, OutputImage>
	::GenerateData()
	{
		//InputImage::IndexType index = {-4, -4, -2};
		//InputImage::SizeType size = {9, 9, 5};
		//InputImage::RegionType iteration_region(index, size);
		//
		//InputImage::Pointer FilterX= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterY= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterZ= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterXX= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterXY= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterXZ= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterYY= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterYZ= AllocateNewImage(iteration_region);;
		//InputImage::Pointer FilterZZ= AllocateNewImage(iteration_region);;

		//std::cerr<<"Before Enter Compute"<<std::endl;
		//ComputeFilter(0.7,0.1,size[0],size[1],size[2],FilterX,FilterY,FilterZ,FilterXX,FilterXY,FilterXZ,FilterYY,FilterYZ,FilterZZ);
		//std::cerr<<"After Compute"<<std::endl;

		typename OutputImageType::Pointer output = GetOutput();
		output->SetRegions(GetInput()->GetRequestedRegion());
		output->Allocate();

		InputImageIteratorType input_iter(GetInput(),GetInput()->GetRequestedRegion());
		OutputImageIteratorType output_iter(GetOutput(),GetInput()->GetRequestedRegion());

		typename InputImageType::Pointer curvature_image = AllocateNewImage(GetInput()->GetRequestedRegion());
		TempImageIteratorType curvature_image_iter(curvature_image, GetInput()->GetRequestedRegion());
		
		typename HessianImageFilterType::Pointer hessian_image_filter = HessianImageFilterType::New();
		hessian_image_filter->SetInput(GetInput());
		hessian_image_filter->SetSigma(1);
		hessian_image_filter->Update();

		typename HessianImageType::Pointer hessian_image = hessian_image_filter->GetOutput();
		HessianImageIteratorType hessian_image_iter(hessian_image, hessian_image->GetRequestedRegion());
		
		int count=0;
		float max=0;
		for (input_iter.GoToBegin(), curvature_image_iter.GoToBegin(), output_iter.GoToBegin(), hessian_image_iter.GoToBegin(); 
			!input_iter.IsAtEnd() && !curvature_image_iter.IsAtEnd() && !output_iter.IsAtEnd() && !hessian_image_iter.IsAtEnd(); 
			++input_iter, ++curvature_image_iter, ++output_iter, ++hessian_image_iter)
		{
			//Derivative temp_data;

			//InputImage::IndexType index = input_iter.GetIndex();

			//InputImage::Pointer subImage=GetSubImage(size[0],size[1], size[2], index);

			////temp_data.x=DoKernel(FilterX, subImage);
			////temp_data.y=DoKernel(FilterY, subImage);
			////temp_data.z=DoKernel(FilterZ, subImage);

			//temp_data.xx=DoKernel(FilterXX, subImage);
			//temp_data.xy=DoKernel(FilterXY, subImage);
			//temp_data.xz=DoKernel(FilterXZ, subImage);

			//temp_data.yy=DoKernel(FilterYY, subImage);
			//temp_data.yz=DoKernel(FilterYZ, subImage);
			//temp_data.zz=DoKernel(FilterZZ, subImage);
			HessianValueType hessian_value = hessian_image_iter.Get();
			vnl_matrix<float> temp_return(3,3);
			//std::cerr<<"Got Here"<<std::endl;
			for (int i=0;i<3;i++) {
				for (int j=0;j<3;j++) {
					temp_return.put(i,j,(hessian_value(i,j)+hessian_value(j,i))/2);
				}
			}

			//std::cerr<<"Got Here"<<std::endl;
			vnl_symmetric_eigensystem<float> eigensystem(temp_return);
			float lamda[3]={eigensystem.get_eigenvalue(2), eigensystem.get_eigenvalue(1), eigensystem.get_eigenvalue(0)};

			for (int i=0;i<3;i++) {
				for (int j=i;j<3;j++) {
					if (abs(lamda[j])<abs(lamda[i])) {
						float temp = lamda[j];
						lamda[j] = lamda[i];
						lamda[i] = temp;
					}
				}
			//	std::cerr<<lamda[i]<<" ";
			}

			float mean =(lamda[0]+lamda[1]+lamda[2])/3;
			float variance = (vnl_math_sqr(lamda[0]-mean)+vnl_math_sqr(lamda[2]-mean)+vnl_math_sqr(lamda[1]-mean))/3;


			//std::cerr<<std::endl;
			//if (abs(lamda[1]/lamda[2])<.01 && abs(lamda[2])>.1) {
			//	//output_iter.Set(1);
			//	output_iter.Set(abs(lamda[1]/lamda[2]));
			//} else 
			
			//if (abs(lamda[0]/lamda[1])<.5 && abs(lamda[0])>.05 && lamda[2]<.1) {
			//	output_iter.Set(abs(lamda[0]/lamda[1]));
			//}
			
			if (abs(lamda[0])>20 && lamda[2]>0) {
				output_iter.Set(abs(lamda[1]/lamda[2]));
			}
			
			if (max<output_iter.Get()) {
				max = output_iter.Get();
			}
			
				//if ((0.5<lamda_1/lamda_2 && lamda_1/lamda_2<2) && (lamda_1>0 && lamda_3/lamda_1<3)) {
					
				//}
			/*	else if (lamda_1>-0.5 && lamda_1<0.5 && (lamda_3-lamda_2)>100) {
					curvature_image_iter.Set(0);
				}*/
			//}

//			Curvature value = GetCurvatures(temp_data);
			count++;

			if (count % 100000==0) {
				std::cerr<<count<<std::endl;
			}
			//curvature_image_iter.Set(sqrtf(max));

			//curvature_image_iter.Set(sqrtf((vnl_math_sqr(value.max)+vnl_math_sqr(value.min))/2));
			//curvature_image_iter.Set(value.mean);
		}
		WriteITK(output, "eigenvalue scale 1.mhd");

		

		for (output_iter.GoToBegin(); 
			!output_iter.IsAtEnd(); 
			++output_iter)
		{
			output_iter.Set(1-output_iter.Get()/max);
		}


		//InputImageIteratorType input_iter(GetInput(),GetInput()->GetRequestedRegion());
		//OutputImageIteratorType output_iter(GetOutput(),GetInput()->GetRequestedRegion());
		//
		//typename DerivativeImageFilterType::Pointer derivative_filter = DerivativeImageFilterType::New();
		//derivative_filter->SetInput(GetInput());
		//derivative_filter->SetOrder(1);

		//derivative_filter->SetDirection(0);
		////derivative_filter->SetUseImageSpacing(true);
		//derivative_filter->Update();

		//typename InputImage::Pointer x_partial = derivative_filter->GetOutput();
		//TempImageIteratorType x_iter(x_partial,x_partial->GetRequestedRegion());


		//derivative_filter.~SmartPointer();
		//derivative_filter = DerivativeImageFilterType::New();
		//derivative_filter->SetInput(GetInput());
		//derivative_filter->SetOrder(1);
		//derivative_filter->SetDirection(1);
		////derivative_filter->SetUseImageSpacing(true);
		//derivative_filter->Update();

		//typename InputImage::Pointer y_partial = derivative_filter->GetOutput();
		//TempImageIteratorType y_iter(y_partial,y_partial->GetRequestedRegion());
		//

		//derivative_filter.~SmartPointer();
		//derivative_filter = DerivativeImageFilterType::New();
		//derivative_filter->SetInput(GetInput());
		//derivative_filter->SetOrder(1);
		//derivative_filter->SetDirection(2);
		////derivative_filter->SetUseImageSpacing(true);
		//derivative_filter->Update();

		//typename InputImage::Pointer z_partial = derivative_filter->GetOutput();
		//TempImageIteratorType z_iter(z_partial,z_partial->GetRequestedRegion());
		//WriteITK(x_partial, "x.mhd");
		//WriteITK(y_partial, "y.mhd");
		//WriteITK(z_partial, "z.mhd");

		//derivative_filter.~SmartPointer();

		//typename InputImageType::Pointer curvature_image = AllocateNewImage(GetInput()->GetRequestedRegion());
		//TempImageIteratorType curvature_image_iter(curvature_image, GetInput()->GetRequestedRegion());

		//typename HessianImageFilterType::Pointer hessian_image_filter = HessianImageFilterType::New();
		//hessian_image_filter->SetInput(GetInput());
		//float low=1;
		//hessian_image_filter->SetSigma(low);
		////hessian_image_filter->SetUseImageSpacing(true);
		//hessian_image_filter->Update();

		//typename HessianImageType::Pointer hessian_image = hessian_image_filter->GetOutput();
		//HessianImageIteratorType hessian_image_iter(hessian_image, hessian_image->GetRequestedRegion());

		//for (input_iter.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(), hessian_image_iter.GoToBegin(), curvature_image_iter.GoToBegin(), output_iter.GoToBegin(); 
		//	!input_iter.IsAtEnd() && !x_iter.IsAtEnd() && !y_iter.IsAtEnd() && !z_iter.IsAtEnd() && !hessian_image_iter.IsAtEnd() && !curvature_image_iter.IsAtEnd() && !output_iter.IsAtEnd(); 
		//	++input_iter, ++x_iter, ++y_iter, ++z_iter, ++hessian_image_iter, ++curvature_image_iter, ++output_iter)
		//{
		//	Derivative temp_data;
		//	

		//	temp_data.x=x_iter.Get();
		//	temp_data.y=y_iter.Get();
		//	temp_data.z=z_iter.Get();

		//	HessianValueType hessian_value = hessian_image_iter.Get();
		//	temp_data.xx=hessian_value(0,0);
		//	temp_data.xy=hessian_value(0,1);
		//	temp_data.xz=hessian_value(0,2);

		//	temp_data.yy=hessian_value(1,1);
		//	temp_data.yz=hessian_value(1,2);

		//	temp_data.zz=hessian_value(2,2);
		//	

		//	Curvature value = GetCurvatures(temp_data);
		//	curvature_image_iter.Set(sqrtf((vnl_math_sqr(value.max)+vnl_math_sqr(value.min))/2));
		//	//curvature_image_iter.Set(value.mean);

		//	output_iter.Set(0);
		//}
		//WriteITK(curvature_image, "min curve.mhd");
		//WriteITK(x_partial, "mean curve.mhd");
		//WriteITK(y_partial, "gaussian curve.mhd");
		//WriteITK(z_partial, "max curve.mhd");

		//for (float i=low+0.5;i<=6;i=i+0.5) {
		//	std::cerr<<"Sigma = "<<i<<std::endl;
		//	hessian_image_filter->SetSigma(i);
		//	hessian_image_filter->Update();

		//	hessian_image = hessian_image_filter->GetOutput();
		//	hessian_image_iter = HessianImageIteratorType(hessian_image, hessian_image->GetRequestedRegion());

		//	for (input_iter.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(), hessian_image_iter.GoToBegin(), curvature_image_iter.GoToBegin(), output_iter.GoToBegin(); 
		//		!input_iter.IsAtEnd() && !x_iter.IsAtEnd() && !y_iter.IsAtEnd() && !z_iter.IsAtEnd() && !hessian_image_iter.IsAtEnd() && !curvature_image_iter.IsAtEnd() && !output_iter.IsAtEnd(); 
		//		++input_iter, ++x_iter, ++y_iter, ++z_iter, ++hessian_image_iter, ++curvature_image_iter, ++output_iter)
		//	{
		//		Derivative temp_data;
		//		temp_data.x=x_iter.Get();
		//		temp_data.y=y_iter.Get();
		//		temp_data.z=z_iter.Get();

		//		HessianValueType hessian_value = hessian_image_iter.Get();
		//		temp_data.xx=hessian_value(0,0);
		//		temp_data.xy=hessian_value(0,1);
		//		temp_data.xz=hessian_value(0,2);

		//		temp_data.yy=hessian_value(1,1);
		//		temp_data.yz=hessian_value(1,2);

		//		temp_data.zz=hessian_value(2,2);

		//		Curvature value = GetCurvatures(temp_data);

		//		float curvature=sqrtf((vnl_math_sqr(value.max)+vnl_math_sqr(value.min))/2);
		//		//float curvature=value.mean;
		//		//output_iter.Set(output_iter.Get()+WeightingFunction(i)*vnl_math_sqr(curvature-curvature_image_iter.Get()));
		//		output_iter.Set(output_iter.Get()+vnl_math_abs(curvature-curvature_image_iter.Get()));

		//		curvature_image_iter.Set(curvature);
		//	}
		//}
		//for (input_iter.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(), hessian_image_iter.GoToBegin(), curvature_image_iter.GoToBegin(), output_iter.GoToBegin(); 
		//		!input_iter.IsAtEnd() && !x_iter.IsAtEnd() && !y_iter.IsAtEnd() && !z_iter.IsAtEnd() && !hessian_image_iter.IsAtEnd() && !curvature_image_iter.IsAtEnd() && !output_iter.IsAtEnd(); 
		//		++input_iter, ++x_iter, ++y_iter, ++z_iter, ++hessian_image_iter, ++curvature_image_iter, ++output_iter)
		//{
		//	output_iter.Set(1-atan(output_iter.Get())/(0.5*PI));
		//	
		//}
	}
}
#endif // itkCurvatureImageFilter_txx