#include "../include/init.h"

using namespace std;



int main(int narg, char** argv)
{

  MPI_Init(&narg, &argv);

  debug_loop();
  
  if(narg != 2) crash ("ONLY INPUT FILE IS REQUIRED");

  string InputFile = argv[1];

  MasterClass_analysis Master(InputFile);

  MPI_Finalize();

  return 0;
}











