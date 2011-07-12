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

void Write2(ImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void Write(ImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void Write(ImageType1D::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ImageType1D >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void Write(ByteImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ByteImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void Write(ByteImageType2D::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ByteImageType2D >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

void Write(LabelImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< LabelImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}

/*
void Write(ShortImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< ShortImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	writer->Update();
}
*/

void Write(FloatImageType::Pointer image, std::string name) {
	typedef itk::ImageFileWriter< FloatImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
	WriterType::SetGlobalWarningDisplay(false);
	
	std::stringstream ss;
	if ( !note.empty() )
		ss << note << "_";

	if ( write_num )
		ss << write_count++ << "_";

	ss << name;

	writer->SetFileName(ss.str().c_str());

	writer->SetInput(image);
	std::cout<<"Writing: "<<ss.str()<<std::endl;
	
	try {
		writer->Update();
	} catch (itk::ExceptionObject &excp) {
		std::cerr << "ExceptionObject caught: " << excp << std::endl;
		system("pause");
	}
} 

void Write(VoxelImageType::Pointer vmap, std::string name) 
{
	VoxelIteratorType vit( vmap, vmap->GetLargestPossibleRegion() );
	
	FloatImageType::Pointer temp = FloatImageType::New();
	temp->SetRegions( vmap->GetLargestPossibleRegion() );
	temp->CopyInformation( vmap );
	temp->Allocate();
	FloatIteratorType tit( temp, temp->GetLargestPossibleRegion() );

	for (vit.GoToBegin(), tit.GoToBegin(); !vit.IsAtEnd() && !tit.IsAtEnd(); ++vit, ++tit)
	{
		 switch(vit.Get()) {
            case Stool:
                tit.Set(1);
                break;
			case Air:
				tit.Set(2);
				break;
			case Tissue:
				tit.Set(3);
				break;
			case Unclassified:
				tit.Set(4);
				break;
			case StoolAir:
				tit.Set(5);
				break;
			case TissueAir:
				tit.Set(6);
				break;
			case TissueStool:
				tit.Set(7);
				break;
			case ThinStool:
				tit.Set(8);
				break;
			default:
                tit.Set(0);
                break;
        }
		
		//tit.Set( floor( (float) tit.Get() * 255/8 ) );
	}

	Write(temp, name);

}

void Write(ArrayImageType::Pointer &v, std::string name)
{
	std::vector<std::string> name_ar = explode( ".", name );

	ArrayImageType::RegionType region = v->GetLargestPossibleRegion();

	ArrayIteratorType v_iter(v,region);

	//int i=1; // only write tissue partial

	for (int i=0; i < ArrayImageType::ImageDimension ; i++)
	{
		FloatImageType::Pointer image = FloatImageType::New();
		image->SetRegions( region );
		image->SetSpacing( v->GetSpacing() );
		image->SetDirection( v->GetDirection() );
		image->CopyInformation( v );
		image->Allocate();
		
		FloatIteratorType image_iter(image,region);

		for (v_iter.GoToBegin(), image_iter.GoToBegin(); !v_iter.IsAtEnd(); ++v_iter, ++image_iter)
			image_iter.Set( /*255 **/ v_iter.Get()[i] );


		std::stringstream ss;
		ss << name_ar[0] << "_" << i << "." << name_ar[1];
		
		Write( image, ss.str() );
	}
}

//template <typename T>
//typename T::Pointer ReadDicom( std::string inDir, std::string outDir, int slice1=-1, int slice2=-1)
//{	
//	// Create reader
//	itk::ImageSeriesReader<T>::Pointer reader = itk::ImageSeriesReader<T>::New();
//    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
//    reader->SetImageIO( dicomIO );
//
//	/*
//	
//	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
//	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
//	
//	nameGenerator->SetDirectory( path );
//	
//	const std::vector< std::string > seriesUID = nameGenerator->GetSeriesUIDs();
//	std::string seriesIdentifier = seriesUID.begin()->c_str();
//	
//	std::vector< std::string > names = nameGenerator->GetFileNames( seriesIdentifier );
//
//	*/
//	
//	// Create regex finder to match file names
//	std::cout << inDir << std::endl;
//
//	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
//	
//	fit->SetDirectory( inDir );
//	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
//	fit->SetSubMatch(1);
//
//	names = fit->GetFileNames();
//	
//	if (slice1 > -1 && slice2 > -1 && slice2 > slice1)
//	{
//		names.erase( names.begin(), names.begin()+slice1);
//
//		unsigned int change = slice2 - slice1;
//
//		names.erase( names.begin() + change, names.end() );
//	} else if (slice1 > -1 && slice2 == -1) {
//		names.erase( names.begin(), names.begin()+slice1-1);
//		names.erase( names.begin()+1, names.end() );
//	} else {
//		std::cout << "Error truncating slices. Restoring full set." << std::endl;
//	}
//
//    reader->SetFileNames( names );
//	try
//	{
//		reader->Update();
//	}
//	catch( itk::ExceptionObject & err )
//	{
//		std::cerr << "Error reading dicom: " << err << std::endl;
//		return 0;
//	}
//
//	// Get dictionary array pointer
//	metaDataDictionaryArray = *(reader->GetMetaDataDictionaryArray());
//
//	std::cout << "Number of filenames: " << names.size() << std::endl;
//	std::cout << "Number of dictionary entries: " << metaDataDictionaryArray.size() << std::endl;
//	
//	T::Pointer output = reader->GetOutput();
//
//	/*typedef itk::MetaDataDictionary DictionaryType;
//	DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
//    output->SetMetaDataDictionary(dictionary);*/
//	
//    return output;
//}
//
//ImageType::Pointer ReadDicom2( std::string inDir, std::string outDir)
//{
//	typedef itk::ImageSeriesReader< ImageType >     ReaderType;
//	typedef itk::GDCMImageIO                        ImageIOType;
//	typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
//
//	ImageIOType::Pointer gdcmIO = ImageIOType::New();
//	NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
//	
//	namesGenerator->SetInputDirectory( inDir );
//
//	const ReaderType::FileNamesContainer & filenames = 
//							namesGenerator->GetInputFileNames();
//
//	// Software Guide : EndCodeSnippet
//
//	unsigned int numberOfFilenames =  filenames.size();
//
//	std::cout << numberOfFilenames << std::endl; 
//	for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
//	{
//		std::cout << "filename # " << fni << " = ";
//		std::cout << filenames[fni] << std::endl;
//	}
//
//	ReaderType::Pointer reader = ReaderType::New();
//
//	reader->SetImageIO( gdcmIO );
//	reader->SetFileNames( filenames );
//
//	try
//	{
//		reader->Update();
//	}
//	catch (itk::ExceptionObject &excp)
//	{
//		std::cerr << "Exception thrown while writing the image" << std::endl;
//		std::cerr << excp << std::endl;
//		//return EXIT_FAILURE;
//	}
//
//	return reader->GetOutput();
//	
//}

//ImageType::Pointer ReadDicom3( std::string inDir, std::string outDir, int slice1=-1, int slice2=-1)
//{	
//	typedef itk::ImageSeriesReader< ImageType >     ReaderType;
// 
//  typedef itk::GDCMImageIO                        ImageIOType;
//  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
//
//  ImageIOType::Pointer gdcmIO = ImageIOType::New();
//  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
//  
//  namesGenerator->SetInputDirectory( inDir );
//
//  const ReaderType::FileNamesContainer & filenames = 
//                            namesGenerator->GetInputFileNames();
//  // Software Guide : EndCodeSnippet
//
//  unsigned int numberOfFilenames =  filenames.size();
//  std::cout << numberOfFilenames << std::endl; 
// 
//  // copy names
//  std::vector<std::string> names;
//  names.resize( filenames.size() );
//
//  for (int i=0; i<names.size(); i++)
//  {
//	names[i] = filenames[i];
//  }
//
//  // reverse names
//  reverse(names.begin(),names.end());
//
//   for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
//    {
//    std::cout << "filename # " << fni << " = ";
//    std::cout << names[fni] << std::endl;
//    }
//
//
//
//  ReaderType::Pointer reader = ReaderType::New();
//
//  reader->SetImageIO( gdcmIO );
//  reader->SetFileNames( names );
//  
//
//  try
//    {
//    // Software Guide : BeginCodeSnippet
//    reader->Update();
//    // Software Guide : EndCodeSnippet
//    }
//  catch (itk::ExceptionObject &excp)
//    {
//    std::cerr << "Exception thrown while writing the image" << std::endl;
//    std::cerr << excp << std::endl;
//    system("pause");
//	//return EXIT_FAILURE;
//    }
//
//  Write(reader->GetOutput(),"inputDicom.nii");
//
//  // Get dictionary array pointer
//  metaDataDictionaryArray = *(reader->GetMetaDataDictionaryArray());
//
//  //const char * outputDirectory = outDir.c_str();
// 
//  //itksys::SystemTools::MakeDirectory( outputDirectory );
//  
// // seriesWriter = SeriesWriterType::New();
//
// // //seriesWriter->SetInput( reader->GetOutput() );
// // seriesWriter->SetImageIO( gdcmIO );
// //
// // namesGenerator->SetOutputDirectory( outputDirectory );
//
// // seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
// // 
// // seriesWriter->SetMetaDataDictionaryArray( 
// //                       reader->GetMetaDataDictionaryArray() );
// // 
// //// try
// ////   {
// ////   seriesWriter->Update();
// ////   }
// //// catch( itk::ExceptionObject & excp )
// ////   {
// ////   std::cerr << "Exception thrown while writing the series " << std::endl;
// ////   std::cerr << excp << std::endl;
// ////   system("pause");
//	//////return EXIT_FAILURE;
// ////   }
//
//  return reader->GetOutput();
//
//}