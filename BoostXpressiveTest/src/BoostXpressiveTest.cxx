#include "itkImage.h"
#include <iostream>
#include <boost/xpressive/xpressive.hpp>

using namespace boost::xpressive;

int main()
{
	typedef itk::Image< unsigned short, 3 > ImageType;
	ImageType::Pointer image = ImageType::New();
	std::cout << "ITK Hello World !" << std::endl;

	std::string hello( "222785561" );

    sregex rex = sregex::compile( "^2+[5-7]+1+$" );
    smatch what;

    if( regex_match( hello, what, rex ) )
    {
		std::cout << !what.empty();
    }
	system("pause");
	return 0;
}