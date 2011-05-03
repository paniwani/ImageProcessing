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
	writer->Update();
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

template <typename T>
typename T::Pointer ReadDicom( std::string path, int slice1=0, int slice2=-1)
{	
	// Create reader
	itk::ImageSeriesReader<T>::Pointer reader = itk::ImageSeriesReader<T>::New();
    itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
    reader->SetImageIO( dicomIO );

	/*
	
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	
	nameGenerator->SetDirectory( path );
	
	const std::vector< std::string > seriesUID = nameGenerator->GetSeriesUIDs();
	std::string seriesIdentifier = seriesUID.begin()->c_str();
	
	std::vector< std::string > names = nameGenerator->GetFileNames( seriesIdentifier );

	*/
	
	// Create regex finder to match file names
	itk::RegularExpressionSeriesFileNames::Pointer fit = itk::RegularExpressionSeriesFileNames::New();
	
	fit->SetDirectory( path );
	fit->SetRegularExpression("[^.]*i([0-9]+).dcm");
	fit->SetSubMatch(1);

	std::vector<std::string> names = fit->GetFileNames();
	
	if (slice2 > 0 && slice2 > slice1)
	{
		names.erase( names.begin(), names.begin()+slice1);
		names.erase( names.begin()+slice2-slice1, names.end() );
	}

    reader->SetFileNames( names );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error reading dicom: " << err << std::endl;
		return 0;
	}
	
	T::Pointer output = reader->GetOutput();
	
	/*
    // Orient all input images into LAI orientation (spine is at top of image)
    itk::OrientImageFilter<T,T>::Pointer orienter = itk::OrientImageFilter<T,T>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI); //LPI
    orienter->SetInput( output );
    orienter->Update();
	output = orienter->GetOutput();
	*/
	
	
    return output;
}