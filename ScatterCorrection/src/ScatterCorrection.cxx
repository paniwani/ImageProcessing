#include "itkImage.h"
#include "scattercorrection.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "scattercorrection.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkRegularExpressionSeriesFileNames.h"
#include "itkOrientImageFilter.h"
#include "itkColonSegmentationFilter.h"

typedef int Int4;
typedef unsigned char SSlice[512][512];

typedef float PixelType;
typedef unsigned char BytePixelType;
typedef itk::Image<PixelType,3> ImageType;
typedef itk::Image<BytePixelType,3> ByteImageType;

typedef itk::ImageRegionIterator<ImageType> IteratorType;
typedef itk::ImageRegionIteratorWithIndex<ByteImageType> IteratorTypeByteWithIndex;

typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleImageFilterType;

int  mean_density_value;
int  background = 0;
float *scale_map;

int pcol,prow,pslice;
int *sphere_no_points;
short sphere_points[100][100][3];
double mask[2*FILTER+1][2*FILTER+1], mask_total;
unsigned short *data_scale16; 

unsigned char *scale_image;


bool truncateOn = false;
unsigned int truncate_ar[2] = {280, 320};

PixelType getValueFromImageArray(Int4 Counter, PixelType *lattice);
void putValueInImageArray(Int4 Counter, float val, PixelType *lattice);

void compute_scale(PixelType *lat16BitCTVolume, float *scale_map, int NI, int NJ, int NK, int xMin, int xMax, int yMin, int yMax);
void scatterCorrection(PixelType *lat16BitCTVolume, SSlice* DAB, int NI, int NJ, int NK); //temparily added by Jiamin
float  do_conv (const int width, const int height, const int depth, float *F, float *image); //temparily added by Jiamin
int compute_scale_image(PixelType *lat16BitCTVolume, float *VoxelSize, SSlice* DAB, int NI, int NJ, int NK, double a,int b);

// NEW FUNCTIONS
ImageType::Pointer ReadDicom( std::string path );

void ReadITK(ImageType::Pointer &img, char * fileName);
void ReadITK(ByteImageType::Pointer &img, char * fileName);

void WriteITK(ImageType::Pointer image, std::string name);
void WriteITK(ByteImageType::Pointer image, std::string name);

void ConvertByteImageTo3DArray(SSlice *& DAB, ByteImageType::Pointer img);
PixelType getValueFromImageArray(int counter, PixelType *img);
void putValueInImageArray(Int4 counter, float val, PixelType *img);

ImageType::Pointer ResampleImage(ImageType::Pointer input, ImageType::SpacingType output_spacing);
std::vector<std::string> explode( const std::string &delimiter, const std::string &str);

// Parameters
float filterP = 0.01;
int scaleP = 2;
int SCALE = 5;
std::string note = "";

int main(int argc, char * argv[])
{

	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " DicomDirectory SCALE filterP";
		system("pause");
		return EXIT_FAILURE;
	}
	std::cerr<<"Started"<<std::endl;

	std::string dataset = argv[1];

	std::vector<std::string> datasetArr = explode( "/", dataset );
	std::string dsname = datasetArr[ datasetArr.size() - 2 ];
	std::cout << "Dataset: " << dsname << std::endl;

	// Set writer prefix
	note = dsname;

	ImageType::Pointer input = ReadDicom( argv[1] );

	// Set parameters
	SCALE = atoi ( argv[2] );
	filterP = atof( argv[3] );
	




	ImageType::SpacingType spacing = input->GetSpacing();
	
	/*
	// Make isotropic
	spacing[1] = spacing[0];
	spacing[2] = spacing[0];
	input = ResampleImage(input, spacing);
	*/

	// Get size and spacing
	float VoxelSize[3];
	VoxelSize[0] = spacing[0];
	VoxelSize[1] = spacing[1];
	VoxelSize[2] = spacing[2];
	ImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();

	ImageType::RegionType region = input->GetLargestPossibleRegion();

	IteratorType input_iter(input,region);

	WriteITK(input,"input.nii");

	// Copy input image
	ImageType::Pointer input2 = ImageType::New();
	input2->SetRegions(region);
	input2->SetSpacing(input->GetSpacing());
	input2->CopyInformation(input);
	input2->Allocate();
	IteratorType input2_iter(input2,region);

	for(input_iter.GoToBegin(), input2_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++input2_iter)
	{
		input2_iter.Set( input_iter.Get() );
	}


	// Segment colon
	typedef itk::ColonSegmentationFilter<ImageType, ByteImageType> ColonSegmentationFilterType;
	ColonSegmentationFilterType::Pointer segmenter = ColonSegmentationFilterType::New();
	segmenter->SetInput( input );
	segmenter->SetBackgroundValue(0);
	segmenter->SetForegroundValue(255);
	segmenter->SetRemoveBoneLung(true);
	segmenter->Update();

	ByteImageType::Pointer colon = segmenter->GetOutput();
	IteratorTypeByteWithIndex colon_iter(colon,region);

	WriteITK(colon,"colon.nii");


	// Shift image so that air is 0
	typedef itk::MinimumMaximumImageCalculator< ImageType > MinimumMaximumImageCalculatorType;
	MinimumMaximumImageCalculatorType::Pointer minMaxCalc = MinimumMaximumImageCalculatorType::New();
	minMaxCalc->SetImage( input );
	minMaxCalc->ComputeMinimum();
	PixelType min = minMaxCalc->GetMinimum();

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		input_iter.Set( input_iter.Get() + abs(min) );
	}

	PixelType *img = input->GetBufferPointer();

	// Create mask
	ByteImageType::Pointer tag_mask = ByteImageType::New();
	tag_mask->SetRegions(region);
	tag_mask->SetDirection(input->GetDirection());
	tag_mask->Allocate();
	tag_mask->FillBuffer(0);
	IteratorTypeByteWithIndex tag_mask_iter(tag_mask,region);

	for(input_iter.GoToBegin(), tag_mask_iter.GoToBegin(), colon_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++tag_mask_iter, ++colon_iter)
	{
		/*ByteImageType::IndexType idx = tag_mask_iter.GetIndex();

		if ( idx[0] > 50 && idx[0] < 512-50 && idx[1] > 50 && idx[1] < 512-50 )
		{
			tag_mask_iter.Set(255);
		}*/

		if (colon_iter.Get() == 255 && input_iter.Get() >= 200+abs(min))
		{
			tag_mask_iter.Set( 255 );
		}
	}

	WriteITK(tag_mask,"mask.nii");

	tag_mask_iter.GoToBegin();
	input_iter.GoToBegin();

	SSlice *DAB;
	ConvertByteImageTo3DArray(DAB, tag_mask);

	// Allocate scale_image
	ByteImageType::Pointer scale_image_itk = ByteImageType::New();
	scale_image_itk->SetDirection(input->GetDirection());
	scale_image_itk->SetRegions(input->GetLargestPossibleRegion());
	scale_image_itk->Allocate();
	scale_image_itk->FillBuffer(0);
	scale_image = scale_image_itk->GetBufferPointer();

	// Compute object scale
	compute_scale_image(img, VoxelSize, DAB, size[0], size[1], size[2], HIST_THRESHOLD,0);

	// Copy and rescale scale image
	typedef itk::ImageDuplicator<ByteImageType> ImageDuplicatorByteType;
	ImageDuplicatorByteType::Pointer duplicator = ImageDuplicatorByteType::New();
	duplicator->SetInputImage( scale_image_itk );
	duplicator->Update();

	typedef itk::RescaleIntensityImageFilter<ByteImageType> RescaleIntensityImageFilterByteType;
	RescaleIntensityImageFilterByteType::Pointer rescaler = RescaleIntensityImageFilterByteType::New();
	rescaler->SetInput( duplicator->GetOutput() );
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->Update();

	std::stringstream ss;
	ss << "scale_image" << "_SCALE_" << SCALE << "_filterP_" << filterP << ".nii";

	WriteITK(rescaler->GetOutput(), ss.str());

	// Scatter correction
	scatterCorrection(img,DAB,size[0],size[1],size[2]);

	// Shift output back to air as -1024
	input_iter = IteratorType(input,region);

	for (input_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter)
	{
		input_iter.Set( input_iter.Get() - abs(min) );
	}

	// Save image
	ss.str("");
	ss << "input_corrected" << "_SCALE_" << SCALE << "_filterP_" << filterP << ".nii";
	WriteITK(input,ss.str());

	// Write change only image
	input_iter = IteratorType(input,region);

	for(input_iter.GoToBegin(), input2_iter.GoToBegin(); !input_iter.IsAtEnd(); ++input_iter, ++input2_iter)
	{
		input2_iter.Set( input2_iter.Get() - input_iter.Get() );	
	}

	ss.str("");
	ss << "scatter_change" << "_SCALE_" << SCALE << "_filterP_" << filterP << ".nii";
	WriteITK(input2,ss.str());

	return 0;
}

void ConvertByteImageTo3DArray(SSlice *& DAB, ByteImageType::Pointer img) 
{
	ByteImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
	BytePixelType *imgPtr = img->GetBufferPointer();

	DAB = new SSlice[ size[2] ];

	for (int z=0; z < size[2] ; z++)
	{
		for (int y=0; y < size[1] ; y++)
		{
			for (int x=0; x < size[0]; x++)
			{
				DAB[z][y][x] = imgPtr[ z*size[0]*size[1] + y*size[0] + x ];
			}
		}
	}
}

ImageType::Pointer ReadDicom( std::string path )
{	
	// Create reader
	itk::ImageSeriesReader<ImageType>::Pointer reader = itk::ImageSeriesReader<ImageType>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	//fit->SetRegularExpression("[^.]*.(.*)");
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	if (truncateOn)
	{
		names.erase( names.begin(), names.begin()+truncate_ar[0] );
		names.erase( names.begin() + (int) (truncate_ar[1]-truncate_ar[0]) , names.end() );
	}

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
		return 0;
	}

	ImageType::Pointer output = reader->GetOutput();

	/*

    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput( output );
    orienter->Update();

	// Set direction cosines to identity matrix
	output = orienter->GetOutput();

	*/
	
	/*
	ImageType::DirectionType direction;
	direction.Fill(0);
	direction[0][0] = 1;
	direction[1][1] = 1;
	direction[2][2] = 1;
	
	output->SetDirection(direction);
	*/

	// Write direction cosines
	ImageType::DirectionType direction = output->GetDirection();

	std::cout << direction[0][0] << " " << direction[1][1] << " " << direction[2][2];

    return output;

}

void ReadITK(ByteImageType::Pointer &img, char * fileName) 
{
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ByteImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	img = reader->GetOutput();
}

void ReadITK(ImageType::Pointer &img, char * fileName) 
{
	std::cout << "Reading " <<  fileName << std::endl;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	img = reader->GetOutput();
}



PixelType getValueFromImageArray(Int4 counter, PixelType *img)
{
	return img[counter];
}

void putValueInImageArray(Int4 counter, float val, PixelType *img)
{
	img[counter] = val;
}

void WriteITK(ImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	name = note + "_" + name;
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

void WriteITK(ByteImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	name = note + "_" + name;
	WriterType::SetGlobalWarningDisplay(false);
	writer->SetFileName(name.c_str());
	writer->SetInput(image);
	std::cout<<"Writing: "<<name<<std::endl;
	writer->Update();
}

ImageType::Pointer ResampleImage(ImageType::Pointer input, ImageType::SpacingType output_spacing)
{
	ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();

	resampleFilter->SetDefaultPixelValue( -1024 );

	ImageType::SpacingType input_spacing = input->GetSpacing();
	ImageType::SizeType input_size = input->GetLargestPossibleRegion().GetSize();

	// Adjusto 
	ImageType::SizeType output_size;

	output_size[0] = input_size[0] * input_spacing[0] / output_spacing[0];
	output_size[1] = input_size[1] * input_spacing[1] / output_spacing[1];
	output_size[2] = input_size[2] * input_spacing[2] / output_spacing[2];

	resampleFilter->SetOutputSpacing( output_spacing );
	resampleFilter->SetOutputOrigin( input->GetOrigin() );
	resampleFilter->SetOutputDirection( input->GetDirection() );

	resampleFilter->SetSize( output_size );
	resampleFilter->SetInput( input );
	resampleFilter->Update();

	return resampleFilter->GetOutput();
}

float  do_conv (const int         width,
                const int         height,
                const int         depth,
                float              *F,
                float       *image);
/********************************************************************************************
 * FUNCTION: compute_scale_image
 * DESCRIPTION:
 * PARAMETERS:
 *     
 * SIDE EFFECTS: None
 * ENTRY CONDITIONS: None
 * RETURN VALUE: 
 * 
 * 
 * EXIT CONDITIONS: None
 * HISTORY:
 *    Created: 03/07/07 by Jiamin Liu
 *    Modified : 
 *
 *******************************************************************************************/
int compute_scale_image(PixelType *lat16BitCTVolume, float* VoxelSize, SSlice* DAB, int NI, int NJ, int NK, double a,int b)
{

  int i,j,k,l,size,error_code, count, **histogram, hist_sum[FEATURES], feature_thr[FEATURES];
  int tti1,tti2,tti3,tti4,tti5,tti6,tti7,tti8,min,max;
  int x,y,z,xx,yy,zz,largest_density_value, ***ppptti1;
  double tt1,tt2;
  double anisotropy_col,anisotropy_row,anisotropy_slice,homogeneity_sigma,inv_scale_sigma,sigma_constant = 1;
  double sigma_mean,sigma_std,sigma_sum;
  FILE *fp_out, *outfile, *s_file;
  int pow_value[FEATURES],diff_value_max[FEATURES];
  Int4 counterGlobal,counterOffset;
  float imgValue, imgValue2;
  int xMin, xMax, yMin, yMax, zMin, zMax;
  
   
  sigma_constant=a;
  background=0;  


  pcol = NI;
  prow = NJ;
  pslice = NK;
  Int4 slice_size = pcol * prow;
  Int4 volume_size = slice_size * pslice;
  
  //*****************************************************************************
  tti1 = 2 * (SCALE + 5);
  ppptti1 = (int ***) malloc(tti1 * sizeof(int **));
    
  ppptti1[0] = (int **) malloc(tti1 * tti1 * sizeof(int *));
    
  for (i = 0; i < tti1; i++)
    ppptti1[i] = ppptti1[0] + i * tti1;
  ppptti1[0][0] = (int *) malloc(tti1 * tti1 * tti1 * sizeof(int));
  
  for (i = 0; i < tti1; i++)
    for (j = 0; j < tti1; j++)
      ppptti1[i][j] = ppptti1[0][0] + (i * tti1 + j) * tti1;
  
  for (i = 0; i < tti1; i++)
    for (j = 0; j < tti1; j++)
      for (k = 0; k < tti1; k++)
	ppptti1[i][j][k] = 0;
  
  sphere_no_points = (int *) malloc((SCALE + 1) * sizeof(int));

  //sphere_points = (void *) malloc((SCALE + 1) * sizeof(void *));
  
  anisotropy_col = VoxelSize[0];
  anisotropy_row = VoxelSize[1];
  anisotropy_slice = VoxelSize[2];
  tt1 = anisotropy_col;
  if (tt1 > anisotropy_row)
    tt1 = anisotropy_row;
  if (tt1 > anisotropy_slice)
    tt1 = anisotropy_slice;
  anisotropy_col = anisotropy_col / tt1;
  anisotropy_row = anisotropy_row / tt1;
  anisotropy_slice = anisotropy_slice / tt1;
  
  tti1 = SCALE + 5;
  if (pslice > 1)
    printf("Anisotropy: slice = %f, row = %f, column = %f\n", anisotropy_slice, anisotropy_row, anisotropy_col);
    
  for (k = 0; k <= SCALE; k++)
  {
	  sphere_no_points[k] = 0;
	  for (i = 0; i < 1; i++) 
		  for (j = -k - 2; j <= k + 2; j++)
			  for (l = -k - 2; l <= k + 2; l++)
				  if (ppptti1[tti1 + i][tti1 + j][tti1 + l] == 0)
				  {
					  tt1 = sqrt(pow(((double) j) * anisotropy_row, 2.0) + pow(((double) l) * anisotropy_col, 2.0));
					  if (tt1 <= ((double) k) + 0.5)
					  {
						  sphere_no_points[k] = sphere_no_points[k] + 1;
						  ppptti1[tti1 + i][tti1 + j][tti1 + l] = 2;
					  }
				  }

//				  sphere_points[k] = (short *) malloc(3 * sphere_no_points[k] * sizeof(short));

				  if (sphere_points[k] == NULL)
				  {
					  printf("Couldn't allocate memory (execution terminated)\n");
					  exit(-1);
				  }

				  tti2 = 0;
				  for (i = 0; i < 1; i++)
					  for (j = -k - 2; j <= k + 2; j++)
						  for (l = -k - 2; l <= k + 2; l++)
							  if (ppptti1[tti1 + i][tti1 + j][tti1 + l] == 2)
							  {
								  ppptti1[tti1 + i][tti1 + j][tti1 + l] = 1;
								  sphere_points[k][tti2][0] = 0;
								  sphere_points[k][tti2][1] = j;
								  sphere_points[k][tti2][2] = l;
								  tti2 = tti2 + 1;
							  }
  }

  
  fflush(stdout);
  free(ppptti1[0][0]);
  free(ppptti1[0]);
  free(ppptti1);
  
  mask_total = 0.0;
  for (yy = -FILTER; yy <= FILTER; yy++)
	  for (xx = -FILTER; xx <= FILTER; xx++)
		  mask[yy + FILTER][xx + FILTER] = 0;

  for (yy = -FILTER; yy <= FILTER; yy++)
	  for (xx = -FILTER; xx <= FILTER; xx++)
	  {
		  tt2 = pow(anisotropy_col * xx, 2.0);
		  tt2 = tt2 + pow(anisotropy_row * yy, 2.0);
		  tt2 = 1 / (1 + tt2);
		  mask[yy + FILTER][xx + FILTER] = tt2;
		  mask_total = mask_total + tt2;
	  }
    
 
  //**************************************************************************************

	  xMin=NI;
	  xMax=0;
	  yMin=NJ;
	  yMax=0;
	  zMin=NK;
	  zMax=0;

	for (k=0,l=0;k<NK;k++ ){
		for (i=0;i<NI;i++ ){
			for (j=0;j<NJ;j++,l++ )
			{
				int val = DAB[k][i][j];
				//if ((val>=InitLumenPhantom && val<InitFluidVol)||(val>=InitFluidVol && val<InitFluidWall))
				//if (val>=InitFluidVol && val<InitFluidWall)
				if (val==255)
				{
					if (j>xMax) xMax=j;
					if (j<xMin) xMin=j;
					if (i>yMax) yMax=i;
					if (i<yMin) yMin=i;
					if (k>zMax) zMax=k;
					if (k<zMin) zMin=k;
				}											
			}
		}
	}

	  tt1 = 0;
	  largest_density_value = 0;
	  for(z = 0;z<pslice;z++) {
		  for(y = 0;y<prow;y++) {
			  for(x = 0;x<pcol;x++)
					{
						counterGlobal=z*(prow*pcol)+y*pcol+x;
						imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
						tt1 = tt1 + imgValue;
						if (imgValue > largest_density_value)
							largest_density_value = imgValue;
					}
		  }
	  }
	  mean_density_value = tt1/volume_size;
	  printf("mean_density_value: %d\n",mean_density_value);

  /**********************************************************************************/
  /*----to compute the histogram for each feature and the threshold for true edge---*/


	  /*****************************************************************

	  SET LARGEST DENSITY VALUE TO MAX-MIN OF IMAGE TO ALLOW ALL POSSIBLE CHANGE VALUES


	  *****************************************************************/
		/*
	  // find min and max of image
	  min = 0;
	  max = 0;


	  for(z = 0;z<pslice;z++) {
		  for(y = 0;y<prow;y++) {
			  for(x = 0;x<pcol;x++)
					{
						counterGlobal=z*(prow*pcol)+y*pcol+x;
						
						int val=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
						
						if (val > max)
							max = val;

						if ( val < min )
							min = val;
					}
		  }
	  }

	  largest_density_value = max - min;
	  */

		



	  histogram = (int **)malloc(FEATURES*sizeof(int *));
	  histogram[0] = (int *)malloc(FEATURES*(largest_density_value+1)*sizeof(int));
	  if(histogram == NULL || histogram[0] == NULL)
		  printf("Memory allocation error \n");
      
	  for(j=0;j<=largest_density_value;j++)
		  histogram[0][j] = 0;

	  for (i=0;i<pslice;i++)
	  {
		  for (j=0;j<prow;j++)
		  {
			  for (k=0;k<pcol-1;k++)
			  {
				  xx = k+1;
				  yy = j;
				  zz = i;
				  Int4 counterGlobal=i*(prow*pcol)+j*prow+k;
				  Int4 counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  float imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  float imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  histogram[0][tti1]++;
			  }
		  }
	  }

	  for (i=0;i<pslice;i++){
		  for (j=0;j<prow-1;j++){
			  for (k=0;k<pcol;k++)
			  {
				  xx = k;
				  yy = j+1;
				  zz = i;
				  counterGlobal=i*(prow*pcol)+j*prow+k;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  histogram[0][tti1]++;
			  }
		  }
	  }

	  for (i=0;i<pslice-1;i++){
		  for(j=0;j<prow;j++){
			  for (k=0;k<pcol;k++)
			  {
				  xx = k;
				  yy = j;
				  zz = i+1;
				  counterGlobal=i*(prow*pcol)+j*prow+k;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  histogram[0][tti1]++;			  }
		  }
	  }

	  hist_sum[0] = 0;
	  for(j=0;j<largest_density_value;j++)
		  hist_sum[0] = hist_sum[0] + histogram[0][j];

	  for(j=0;j<largest_density_value;j++)
	  {
		  tti1 = 0;
		  feature_thr[0] = (double)j;
		  for(k=0;k<=j;k++)
			  tti1 = tti1+histogram[0][k];
		  if (((double)tti1 /(double) hist_sum[0])>=HIST_THRESHOLD)
			  break;
	  }

	  // Write histogram to text file
	  std::ofstream file;

	  std::stringstream ss;
	  ss << "histogram_" << SCALE << ".csv";

	  file.open( ss.str().c_str() );

	  file << "Bin,Value\n";

	  for(j=0;j<=largest_density_value;j++)
		  file << j << "," << histogram[0][j] << "\n";

	  file.close();

      printf("Histogram threshold computation is done \n");
      printf("Features Threshold %d : %f \n", i,(double)feature_thr[0]); 
      
      //--------------------compute the homogeneity sigma-------------------------------
      sigma_sum = 0;
      count = 0;
      sigma_mean = 0;
	 
	  for(z = 0;z<pslice;z++){
		  for(y = 0;y<prow;y++){
			  for(x = 0;x<pcol-1;x++)
			  {
				  zz = z;
				  yy = y;
				  xx = x + 1;

				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  			  
				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + tti1;
					  count ++;
				  }
			  }
		  }
	  }

	  for(z = 0;z<pslice;z++){
		  for(y = 0;y<prow-1;y++){
			  for(x = 0;x<pcol;x++)
			  {
				  zz = z;
				  yy = y+1;
				  xx = x;
				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  
				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + tti1;
					  count ++;
				  }
			  }
		  }
	  }
	  for(z = 0;z<pslice-1;z++){
		  for(y = 0;y<prow;y++){
			  for(x = 0;x<pcol;x++)
			  {
				  zz = z+1;
				  yy = y;
				  xx = x;
				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  
				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + tti1;
					  count ++;
				  }
			  }
		  }
	  }

	  sigma_mean = sigma_sum / count;
	  printf("homogeneity_mean value is: %f \n", sigma_mean);
      //--------------------compute the homogeneity sigma-------------------------------
      sigma_sum = 0;
      count = 0;
	
	  for(z = 0;z<pslice;z++){
		  for(y = 0;y<prow;y++){
			  for(x = 0;x<pcol-1;x++)
			  {
				  zz = z;
				  yy = y;
				  xx = x + 1;
				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);

				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + pow((tti1-sigma_mean),2);
					  count ++;
				  }
			  }
		  }
	  }

	  for(z = 0;z<pslice;z++){
		  for(y = 0;y<prow-1;y++){
			  for(x = 0;x<pcol;x++)
			  {
				  zz = z;
				  yy = y+1;
				  xx = x;
				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  
				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + pow((tti1 - sigma_mean),2);
					  count ++;
				  }
			  }
		  }
	  }
		
	  for(z = 0;z<pslice-1;z++){
		  for(y = 0;y<prow;y++){
			  for(x = 0;x<pcol;x++)
			  {
				  zz = z+1;
				  yy = y;
				  xx = x;
				  counterGlobal=z*(prow*pcol)+y*prow+x;
				  counterOffset=zz*(prow*pcol)+yy*prow+xx;
				  imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				  imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
				  tti1 = abs(imgValue - imgValue2);
				  
				  if(tti1 <= feature_thr[0])
				  {
					  sigma_sum = sigma_sum + pow((tti1 - sigma_mean),2);
					  count ++;
				  }
			  }
		  }
	  }

	sigma_std = sqrt((double)sigma_sum/(double)count);
	homogeneity_sigma = sigma_constant * (sigma_mean + 3*sigma_std);
	printf("homogeneity_sigma value: %f \n",homogeneity_sigma); 

	scale_map = (float *) malloc( (largest_density_value + 1) * sizeof(double));
	
	//------------------ GAUSSIAN ----------------------------------------------------
	inv_scale_sigma = -0.5 / pow(homogeneity_sigma, 2.0);
	for (i = 0; i <= largest_density_value; i++){
	  scale_map[i] = exp(inv_scale_sigma * pow((double) i, 2.0));
    }
  
  	

	
 compute_scale(lat16BitCTVolume, scale_map, NI, NJ, NK, xMin, xMax, yMin, yMax);
    
  /*fp_out = fopen("scaleimage.raw","w+");
  fwrite(scale_image, sizeof(unsigned char), volume_size, fp_out);
  fclose(fp_out);*/

  free(histogram[0]);
  free(histogram);
  // free(scale_image);

  free(sphere_no_points);

  /*for(i=0;i<=SCALE;i++)
    free(sphere_points[i]);
  free(sphere_points);*/

  printf("scale end normally \n");
  
}
/*****************************************************************************
 * FUNCTION: compute_scale
 * DESCRIPTION: Computes the scale values for the entire volume anf store in the 
 *        scale-image array.
 * PARAMETERS: None
 * SIDE EFFECTS: 
 * FUNCTIONS CALEED: None
 * ENTRY CONDITIONS: 1) scale_map array is alloted and 
 *           proper values are assigned
 * RETURN VALUE: None
 * EXIT CONDITIONS: Compute scale values
 * HISTORY:
 *  Created: 02/24/00
 *  Modified:07/25/00 extend to 24 bits color image by Jiamin  
 *
 *****************************************************************************/
void compute_scale(PixelType *lat16BitCTVolume, float *scale_map, int NI, int NJ, int NK, int xMin, int xMax, int yMin, int yMax)
{
	int i, j, k, x, y, z, xx, yy, zz, mean_g, tti5, slice, row, col;
	int flag, tti1,tti2,edge_flag;
	double inv_scale_sigma, count_obj, count_nonobj, tt1, tt2, tt3, mask_f[FEATURES];
	int mean[FEATURES],temp[FEATURES];
	Int4 counterGlobal, counterOffset;
	float imgValue, imgValue2;

	
	Int4 slice_size = NI * NJ;
	Int4 volume_size = slice_size * NK;

	//scale_image = (unsigned char *) malloc(volume_size * sizeof(unsigned char));

	for (slice = 0; slice < pslice; slice++)
		for (row = yMin; row < yMax; row++)
			for (col = xMin; col < xMax; col++)
		//for (row = 0; row < prow; row++)
			//for (col = 0; col < pcol; col++)
			{
				
				counterGlobal=slice*(prow*pcol)+row*prow+col;
				imgValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
				tti1 = imgValue;
				if((background == 0)&&(tti1<mean_density_value))
					scale_image[slice * slice_size + row * pcol + col] = 0;
				else
				{
					flag = 0;
					tt1 = 0.0;
					tt3 = 0.0;


					for (yy = -FILTER; yy <= FILTER; yy++)
						for (xx = -FILTER; xx <= FILTER; xx++)
						{
							x = xx + col;
							y = yy + row;
							z = slice;

							counterOffset=z*(prow*pcol)+y*prow+x;
							imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
							imgValue=getValueFromImageArray(counterOffset,lat16BitCTVolume);

							if (x >= 0 && y >= 0  && x < pcol && y < prow)
								tt3 = tt3 + mask[yy + FILTER][xx + FILTER]	* (double) imgValue2;
							else
								tt3 = tt3 + mask[yy + FILTER][xx + FILTER] 
							* (double) imgValue;
						}
						mean_g = (int) (tt3 / mask_total + 0.5);

						for (k = 1; k < SCALE && !flag; k++)  
						{
							count_obj = 0;
							count_nonobj = 0;
							for (i = 0; i < sphere_no_points[k]; i++) 
							{
								x = col + sphere_points[k][i][2];
								y = row + sphere_points[k][i][1];
								z = slice + sphere_points[k][i][0];
								if (x < 0 || x >= pcol)
									x = col;
								if (y < 0 || y >= prow)
									y = row;
								if (z < 0 || z >= pslice)
									z = slice;

								counterOffset=z*(prow*pcol)+y*prow+x;
								imgValue2=getValueFromImageArray(counterOffset,lat16BitCTVolume);
								tti5 = (int) imgValue2;
								tti5 = tti5 - mean_g;
								if (tti5 < 0)
									tti5 = -tti5;
								count_obj = count_obj + scale_map[tti5];
								count_nonobj = count_nonobj + 1.0 - scale_map[tti5];
							}
							if (100.0 * count_nonobj >= tolerance
								* (count_nonobj + count_obj)) {
									scale_image[slice * slice_size + row * pcol + col] = k;
									flag = 1;
							}
						}
						if (!flag)
							scale_image[slice * slice_size + row * pcol + col] = k;
				}
			}

			if(scale_map)
				free(scale_map);

			printf("\rScale computation is done.     \n"); fflush(stdout);	
}


void RemoveFluidInCT(PixelType *lat16BitCTVolume,SSlice* DAB,int NI, int NJ, int NK)
{
	int k,l,i,j;


	for (k=0,l=0;k<NK;k++ ){
		for (i=0;i<NI;i++ ){
			for (j=0;j<NJ;j++,l++ )
			{

				int val = DAB[k][i][j];
				//if ((val>=InitLumenPhantom && val<InitFluidVol)||(val>=InitFluidVol && val<InitFluidWall))
				if (val==255)
				{
					Int4 counterGlobal=k*(NI*NJ)+i*NJ+j;
					float PhantomValue= (((float) rand() / (float) RAND_MAX) * 50 + 10);
					putValueInImageArray(counterGlobal,PhantomValue,lat16BitCTVolume);

				}											
			}
		}
	}

}

//added by jiamin
void gaussFilter(int size_sf, float std_sf, float *scatterFilter)
{
	
	float sum=0;
	for (int i=0;i<size_sf;i++){
		for (int j=0;j<size_sf; j++)
		{
			float x=(i-(size_sf+1)/2)^2;
			float y=(j-(size_sf+1)/2)^2;
			scatterFilter[i*size_sf+j]=exp(-(x+y)/(2*std_sf*std_sf));
			sum=sum+scatterFilter[i*size_sf+j];
		}
	}

	for (int i=0;i<size_sf;i++){
		for (int j=0;j<size_sf; j++)
		{
			scatterFilter[i*size_sf+j]=filterP*scatterFilter[i*size_sf+j]/sum;
			
		}
	}
	
}
void scatterCorrection(PixelType *lat16BitCTVolume, SSlice* DAB, int NI, int NJ, int NK)
{

	short i,j,k,l,index=0;
	Int4 counterSub,counterGlobal;
	int xMin=1000,xMax=0,yMin=1000,yMax=0,zMin=1000,zMax=0;

	float* scatterFilter;
	float *subImage;

	/*for (k=0,l=0;k<600;k++ )
		for (i=0;i<NI;i++ )
			for (j=0;j<NJ;j++,l++ )
				segmentedFluid[k][i][j]=0;*/

		
	for (k=0,l=0;k<NK;k++ ){
		for (i=0;i<NI;i++ ){
			for (j=0;j<NJ;j++,l++ )
			{
				
				int val = DAB[k][i][j];
				//if ((val>=InitLumenPhantom && val<InitFluidVol)||(val>=InitFluidVol && val<InitFluidWall))
				//if (val>=InitFluidVol && val<InitFluidWall)
				if (val==255)
				{
					if (j>xMax) xMax=j;
					if (j<xMin) xMin=j;
					if (i>yMax) yMax=i;
					if (i<yMin) yMin=i;
					if (k>zMax) zMax=k;
					if (k<zMin) zMin=k;
					
				}											
			}
		}
	}
	float inputValue;
	float scatterValue;

	for (k=zMin,l=0;k<zMax;k++ )
		for (i=yMin;i<yMax;i++ )
			for (j=xMin;j<xMax;j++,l++ )
			{

				int size_sf=scale_image[k * NI*NJ + i * NJ + j];

				if (size_sf >0) 
				{
					size_sf=scaleP*(SCALE-size_sf+1)+1;
					if (size_sf<=0) size_sf=1;

					int PartFullorE=0;

					for (int indexI=i-(size_sf-1)/2;indexI<=i+(size_sf-1)/2;indexI++){
						for (int indexJ=j-(size_sf-1)/2;indexJ<=j+(size_sf-1)/2;indexJ++)
						{
							int val = DAB[k][indexI][indexJ];
							//if (val>=InitFluidVol && val<InitFluidWall)
							//if ((val>=InitLumenPhantom && val<InitFluidVol)||(val>=InitFluidVol && val<InitFluidWall))
							if (val==255)	
								PartFullorE++;

						}
					}

					if ((PartFullorE<size_sf*size_sf) && (PartFullorE>0.2*size_sf*size_sf))
					{
						scatterFilter = new float[size_sf*size_sf];
						int std_sf=size_sf;
						gaussFilter(size_sf,std_sf,scatterFilter);

						index=0;
						subImage = new float[size_sf*size_sf];
						for (int indexI=i-(size_sf-1)/2;indexI<=i+(size_sf-1)/2;indexI++){
							for (int indexJ=j-(size_sf-1)/2;indexJ<=j+(size_sf-1)/2;indexJ++)
							{
								int val = DAB[k][indexI][indexJ];
								counterGlobal=k*(NI*NJ)+i*NJ+j;
								//if ((val>=InitLumenPhantom && val<InitFluidVol)||(val>=InitFluidVol && val<InitFluidWall))
								/*if (val>=InitFluidVol && val<InitFluidWall)
								{*/
								counterSub=k*(NI*NJ)+indexI*NJ+indexJ;			
								subImage[index]=getValueFromImageArray(counterSub,lat16BitCTVolume);
								/*}
								else
								subImage[index]=0;*/
								index++;
							}
						}

						counterGlobal=k*(NI*NJ)+i*NJ+j;
						inputValue=getValueFromImageArray(counterGlobal,lat16BitCTVolume);
						scatterValue=do_conv(size_sf,size_sf,1,scatterFilter,subImage);
						delete scatterFilter;
						delete subImage;

						//int val_t = DAB[k][i][j];
						//////scatter is not corrected in the fluid region
						//if ((val_t>=InitFluidVol) && (val_t<InitFluidWall))
						//{
						//}
						//else
						putValueInImageArray(counterGlobal,inputValue-scatterValue,lat16BitCTVolume); //scatter correction will affect local region of fluid.
					}
				}
				
			}
}

float  do_conv (const int         width,
                const int         height,
                const int         depth,
                float              *F,
                float       *image)
{
    int numVoxels = width*height*depth;
    float * itImage = image;
    float * itFilter = F;
    float * itImageEnd = image + numVoxels;
    float weight_pixel_sum = 0;
    for ( ; itImage != itImageEnd;  ++itImage, ++itFilter)
	{
        weight_pixel_sum += (*itImage) * (*itFilter);
    }
    return weight_pixel_sum;
} // do_conv()

std::vector<std::string> explode( const std::string &delimiter, const std::string &str)
{
	std::vector<std::string> arr;

    int strleng = str.length();
    int delleng = delimiter.length();
    if (delleng==0)
        return arr;//no change

    int i=0; 
    int k=0;
    while( i<strleng )
    {
        int j=0;
        while (i+j<strleng && j<delleng && str[i+j]==delimiter[j])
            j++;
        if (j==delleng)//found delimiter
        {
            arr.push_back(  str.substr(k, i-k) );
            i+=delleng;
            k=i;
        }
        else
        {
            i++;
        }
    }
    arr.push_back(  str.substr(k, i-k) );
    return arr;
}