#include "Timer.h"
#include "BTCCollider.h"
#include "SECP256k1.h"
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdexcept>
#include "hash/sha512.h"
#include "hash/sha256.h"

#define RELEASE "1.2"

using namespace std;

// Function definitions for printUsage(), getInt(), and getInts() remain unchanged

int main(int argc, char* argv[]) {
  // Global Init
  Timer::Init();
  rseed((unsigned long)time(NULL));

  // Init SecpK1
  Secp256K1 *secp = new Secp256K1();
  secp->Init();

  int a = 1;
  int dp = -1;
  bool gpuEnable = false;
  bool stop = false;
  vector<int> gpuId = {0};
  vector<int> gridSize;
  string seed = "";
  vector<string> prefix;
  string outputFile = "";
  int nbCPUThread = Timer::getCoreNumber();
  bool tSpecified = false;
  bool extraPts = false;
  bool checkFlag = false;
  uint32_t cSize = 40;
  uint64_t rekey = 0;
  Point startPuKey;
  startPuKey.Clear();
  string workFile = "";
  string iWorkFile = "";
  uint32_t savePeriod = 60;

  // Command line argument parsing remains unchanged

  printf("BTCCollider v" RELEASE "\n");

  if(gridSize.size()==0) {
    for (int i = 0; i < gpuId.size(); i++) {
      gridSize.push_back(0);
      gridSize.push_back(0);
    }
  } else if(gridSize.size() != gpuId.size()*2) {
    printf("Invalid gridSize or gpuId argument, must have coherent size\n");
    exit(-1);
  }

  // Let one CPU core free per gpu if gpu is enabled
  // It will avoid hanging the system
  if( !tSpecified && nbCPUThread>1 && gpuEnable)
    nbCPUThread-=(int)gpuId.size();
  if(nbCPUThread<0)
    nbCPUThread = 0;

  BTCCollider *v = new BTCCollider(secp, gpuEnable, stop, outputFile, workFile, iWorkFile, savePeriod, cSize, dp, extraPts);

  if(checkFlag)
    v->Check(gpuId, gridSize);
  else
    v->Search(nbCPUThread,gpuId,gridSize);

  delete v;
  delete secp;

  return 0;
}
