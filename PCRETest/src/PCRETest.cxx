#include "itkImage.h"
#include <iostream>
#include <pcrecpp.h>

#define PCRE_STATIC

int main()
{
	typedef itk::Image< unsigned short, 3 > ImageType;
	ImageType::Pointer image = ImageType::New();
	std::cout << "ITK Hello World !" << std::endl;

	pcrecpp::RE re("h.*o");
	re.FullMatch("hello");


	system("pause");
	return 0;
}