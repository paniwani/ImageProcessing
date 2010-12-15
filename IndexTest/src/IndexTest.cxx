#include "itkImage.h"
#include <iostream>
int main()
{
	typedef itk::Image< unsigned short, 3 > ImageType;
	
	ImageType::IndexType a,b;
	a[0]=0;
	a[1]=1;
	a[2]=2;

	std::vector< ImageType::IndexType > vec;
	vec.push_back( a );

	b = vec[0];

	vec.clear();

	std::cout << b[0] << std::endl;
	std::cout << b[1] << std::endl;
	std::cout << b[2] << std::endl;






	system("pause");
	return 0;
}