#include "itkImage.h"
#include <iostream>
int main()
{
typedef itk::Image< unsigned short, 3 > ImageType;
ImageType::Pointer image = ImageType::New();
std::cout << "ITK Hello World !" << std::endl;
system("pause");
return 0;
}