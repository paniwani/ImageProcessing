#include <itkImage.h> 						
#include <iostream>
#include <utils.h>
#include <itkBinaryShapeKeepNObjectsImageFilter.h>
#include <time.h>

ByteImageType::Pointer AllocateByteImage(ImageType::Pointer &in)
{
	ByteImageType::Pointer out = ByteImageType::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(0);
	return out;
}

ByteImageType::Pointer BinaryThreshold(ImageType::Pointer &in, PixelType t1=itk::NumericTraits<PixelType>::NonpositiveMin(), 
									   PixelType t2=itk::NumericTraits<PixelType>::max())
{
	ByteImageType::Pointer out = AllocateByteImage(in);
	
	IteratorType it(in,in->GetLargestPossibleRegion());
	ByteIteratorType oit(out,out->GetLargestPossibleRegion());

	for (it.GoToBegin(), oit.GoToBegin(); !it.IsAtEnd(); ++it, ++oit)
	{
		if ( it.Get() > t1 && it.Get() < t2)
		{
			oit.Set(255);
		}
	}
	
	return out;
}

void BinaryKeeper(ByteImageType::Pointer &in, std::string attribute, unsigned int numOfObjects, bool reverse=false)
{
	typedef itk::BinaryShapeKeepNObjectsImageFilter<ByteImageType> KeeperType;
	KeeperType::Pointer keeper = KeeperType::New();
	keeper->SetInput(in);
	keeper->SetAttribute(attribute);
	keeper->SetNumberOfObjects(numOfObjects);
	keeper->SetReverseOrdering(reverse);
	keeper->SetForegroundValue(255);
	keeper->SetBackgroundValue(0);
	keeper->Update();
	in = keeper->GetOutput();
}
 												
int main(int argc, char * argv[])				
{ 	
	clock_t init;
	init=clock();

	ImageType::Pointer input = ReadDicom <ImageType> ("C:/imagedata/mr10-uncleansed/mr10_092_13p.i0344/dcm");
	
	std::cout << (double) (clock()-init) / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;
	init=clock();

	WriteITK <ImageType> (input,"input.nii");

	std::cout << (double) (clock()-init) / ((double) CLOCKS_PER_SEC) <<  " seconds" << std::endl;

	/*ByteImageType::Pointer bone = BinaryThreshold(input,400);
	BinaryKeeper(bone,"Size",1);

	WriteITK <ByteImageType> (bone,"bone.nii");*/
	return 0; 									
} 												
