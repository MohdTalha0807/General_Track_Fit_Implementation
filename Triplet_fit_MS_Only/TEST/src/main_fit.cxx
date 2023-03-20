#include <iostream>
#include <string>

int main (int argc, char * argv[])
{

  if (argc < 3)
  {
    std::cout << "[ERROR] Not enough arguments provided!" << std::endl;
    std::cout << "        Call like this:"                << std::endl;
    std::cout << "        ./name_of_script PATH/2/INPUT/file.root PATH/2/OUTPUT/file.root" << std::endl;
    return 0;
  }

  std::cout << "[INFO] Arguments to program:" << std::endl;
  std::cout << "       |- Input file:  " << argv[1] << std::endl;
  std::cout << "       |- Output file: " << argv[2] << std::endl;

  return 0;

}
