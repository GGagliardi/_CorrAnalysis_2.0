#include "../include/init.h"

using namespace std;



int main(int narg, char** argv)
{
  if(narg != 2) crash ("ONLY INPUT FILE IS REQUIRED");

  string InputFile = argv[1];

  MasterClass_analysis Master(InputFile);

  return 1;
}











