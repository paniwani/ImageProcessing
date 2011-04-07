#ifndef FastDilateFilter_txx
#define FastDilateFilter_txx

#include "itkFastDilateFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage >
FastDilateFilter<TInputImage,TOutputImage >
::FastDilateFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_NumOfIterations = 1;
  m_FullyConnected = false;
}


/**
 * GenerateData Performs the reconstruction
 */
template <class TInputImage, class TOutputImage >
void
FastDilateFilter<TInputImage,TOutputImage>
::GenerateData( void )
{
	// Setup
	InputImagePointer input = const_cast<InputImageType *>(this->GetInput());
	typename Superclass::OutputImagePointer output = this->GetOutput(0);

	InputImageType::RegionType region = input->GetLargestPossibleRegion();
	InputImageType::SizeType size = region.GetSize();
	InputImageType::SpacingType spacing = input->GetSpacing();

	output->SetRegions( region );
	output->SetOrigin( input->GetOrigin() );
	output->SetSpacing( input->GetSpacing() );
	output->CopyInformation( input );
	output->Allocate();	
	
	typedef ImageRegionIterator<InputImageType> InputIteratorType;
	typedef ImageRegionIterator<OutputImageType> OutputIteratorType;
	
	int x, y, z, sizex, sizey, sizez, sizexy, i, k;
    sizex = size[0];
    sizey = size[1];
    sizez = size[2];
    sizexy = sizex*sizey;
    
	int march[8], n, connect;
    march[0]=-1; march[1]=1; march[2]=-sizex; march[3]=sizex;
    march[4]=-sizex-1; march[5]=-sizex+1; march[6]=sizex-1; march[7]=sizex+1;
	
	if(m_FullyConnected) connect=8;
    else connect=4;

    // now only for binary image
    for(i=0; i<m_NumOfIterations; i++)
    {
		// Copy input to output
        InputIteratorType input_iter(input,region);
		OutputIteratorType output_iter(output,region);
		
		for (input_iter.GoToBegin(), output_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++output_iter)
		{
			output_iter.Set( input_iter.Get() );
		}

		InputImagePixelType *imArray = input->GetBufferPointer();
		InputImagePixelType *omArray = output->GetBufferPointer();
		
		// March to neighbors
        omArray = output->GetBufferPointer();
        for(z=1; z<sizez-1; z++)
        {
            for(y=1; y<sizey-1; y++)
            {
                k = z*sizexy+y*sizex;
                for(x=1; x<sizex-1; x++)
                {
                    k++;
                    if(omArray[k]!=0)
                    {
                        for(n=0; n<connect; n++) 
                            if(omArray[k+march[n]]==0) imArray[k+march[n]] = omArray[k];
                    }   // if oim
                }   // for x
            }   // for y
        }   // for z
    }   // for i

    output = input;
	
}

template <class TInputImage, class TOutputImage >
void
FastDilateFilter<TInputImage,TOutputImage>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "NumOfIterations: " << m_NumOfIterations << std::endl;
  os << indent << "FullyConnected: " << m_FullyConnected << std::endl;
}

} // end namespace itk

#endif
