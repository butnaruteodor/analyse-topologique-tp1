#include <DGtal/base/Common.h>
 
int main(int argc, char** argv)
{
  
  DGtal::trace.info() << "Helloworld from DGtal ";
  DGtal::trace.emphase() << "(version "<< DGtal_VERSION << ")"<< std::endl;
  
  return 0;
}