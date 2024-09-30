/*
 * This file is part of the BTCCollider distribution (https://github.com/solimanko/coli).
 * Copyright (c) 2024 Abdullah Soliman Dev.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "BTCCollider.h"
#include "Base58.h"
#include "Bech32.h"
#include "hash/sha256.h"
#include "hash/sha512.h"
#include "IntGroup.h"
#include "Timer.h"
#include "hash/ripemd160.h"
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream>

#ifndef WIN64
#include <pthread.h>
#endif

using namespace std;

#ifdef WIN64
DWORD WINAPI _InitKey(LPVOID lpParam) {
#else
void *_InitKey(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->InitKey(p);
  return 0;
}

// ----------------------------------------------------------------------------

BTCCollider::BTCCollider(Secp256K1 *secp, bool useGpu, bool stop, std::string outputFile, std::string workFile, 
                         std::string iWorkFile, uint32_t savePeriod, uint32_t n, int dp,bool extraPoints) {
  this->secp = secp;
  this->useGpu = useGpu;
  this->outputFile = outputFile;
  this->useSSE = useSSE;
  this->nbGPUThread = 0;
  this->colSize = n;
  this->CPU_GRP_SIZE = 128;
  this->initDPSize = dp;
  this->extraPoints = extraPoints;
  this->workFile = workFile;
  this->saveWorkPeriod = savePeriod;

  if (iWorkFile.length() > 0) {
    LoadWork(iWorkFile);
  } else {
    // Seed
    initialSeed = Timer::getSeed(32);
    offsetCount = 0;
    nbLoadedWalk = 0;
    fetchedWalk = 0;
    offsetTime = 0.0;
    loadedX=NULL;
    loadedY=NULL;
  }

  seed.SetInt32(0);
  sha256((uint8_t *)initialSeed.c_str(), (int)initialSeed.length(), (uint8_t *)seed.bits64);

  printf("Collision: %d bits\n", this->colSize);
  printf("Seed: %s\n", seed.GetBase16().c_str());

  // Derived from pairgen (https://github.com/basil00/pairgen.git)
  // Given a hash160 H comprised on {h0, h1, .., h9} 16-bit parts, then H is
  // mapped to a public key P as follows:
  //     P = pub[0][h0] + pub[1][h1] + .. + pub[9][h9]
  // The calculation is truncated according to the length n.  The corresponding
  // private key P' is:
  //     P' = priv[0]+h0 + priv[1]+h1*2^16 + .. + priv[9]+h9*2^144
  // Each base private key is chosen randomly and computed in advanced. 

  TH_PARAM params[10];
  memset(params, 0, sizeof(params));
  printf("Initializing:");
#ifndef WIN64
  fflush(stdout);
#endif
  Point GP = secp->G;
  Int KP;
  KP.SetInt32(1);
  for (int i = 0; i < 10; i++) {
    Gp[i] = GP;
    Kp[i] = KP;
    for (int j = 0; j < 16; j++)
      GP = secp->DoubleDirect(GP);
    KP.ShiftL(16);
  }

  size_t kSize = 10*65536*2*sizeof(Int);
  pub = (Int *)malloc(kSize);
  if (pub == NULL) {
    printf("Cannot allocate memory for input keys !\n");
    exit(0);
  }

  THREAD_HANDLE threadIDs[10];
  for (int i = 0; i < 10; i++) {
    params[i].threadId = i;
    Rand(&seed,&params[i].localSeed);
    threadIDs[i] = LaunchThread(_InitKey,params+i);
  }

  JoinThreads(threadIDs,10);
  FreeHandles(threadIDs,10);
  printf("Done\n");

  nbFull = colSize / 16 ; // Number of full word
  int leftBit = colSize % 16;
  colMask = (1 << (16 - leftBit)) - 1;
  colMask = ~colMask;
#ifdef WIN64
  colMask = _byteswap_ushort(colMask);
#else
  colMask = __builtin_bswap16(colMask);
#endif
  hashTable.SetParam(colSize,nbFull,colMask);
  
  beta1.SetBase16("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee");
  lambda1.SetBase16("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72");
  beta2.SetBase16("851695d49a83f8ef919bb86153cbcb16630fb68aed0a766a3ec693d68e6afa40");
  lambda2.SetBase16("ac9c52b33fa3cf1f5ad9e3fd77ed9ba4a880b9fc8ec739c2e0cfc810b51283ce");

  char *ctimeBuff;
  time_t now = time(NULL);
  ctimeBuff = ctime(&now);
  printf("Start %s", ctimeBuff);
}

// ----------------------------------------------------------------------------

void BTCCollider::LoadWork(string fileName) {
  FILE *f = fopen(fileName.c_str(), "rb");
  if (f == NULL) {
    printf("LoadWork: Cannot open %s for reading\n", fileName.c_str());
    exit(-1);
  }

  // Read global param
  char tmpS[65];
  tmpS[64] = 0;
  fread(tmpS, 64, 1, f);
  initialSeed = string(tmpS);
  fread(&this->colSize, sizeof(uint32_t), 1, f);
  fread(&this->initDPSize, sizeof(uint32_t), 1, f);
  fread(&this->extraPoints, sizeof(bool), 1, f);
  fread(&this->offsetCount, sizeof(uint64_t), 1, f);
  fread(&this->offsetTime, sizeof(double), 1, f);

  // Read random walks state
  fread(&this->nbLoadedWalk, sizeof(uint64_t), 1, f);
  loadedX = (hash160_t *)malloc(nbLoadedWalk * sizeof(hash160_t));
  loadedY = (hash160_t *)malloc(nbLoadedWalk * sizeof(hash160_t));
  if (loadedX == NULL || loadedY == NULL) {
    printf("LoadWork: Failed to allocate memory !\n");
    exit(-1);
  }
  for (uint64_t n = 0; n < nbLoadedWalk; n++) {
    fread(&loadedX[n], sizeof(hash160_t), 1, f);
    fread(&loadedY[n], sizeof(hash160_t), 1, f);
  }

  // Read hash table
  hashTable.LoadTable(f);
  fclose(f);

  fetchedWalk = 0;

  // Need to reseed here otherwise we may perform the same work again
  string newSeed = Timer::getSeed(32);
  seed.SetInt32(0);
  sha256((uint8_t *)newSeed.c_str(), (int)newSeed.length(), (uint8_t *)seed.bits64);

  printf("LoadWork: %s ok [2^%.2f walks]\n", fileName.c_str(), log2((double)nbLoadedWalk));
}

// ----------------------------------------------------------------------------

void BTCCollider::SaveWork(uint64_t totalCount, double totalTime, TH_PARAM *threads, int nbThread) {
  FILE *f = fopen(workFile.c_str(), "wb");
  if (f == NULL) {
    printf("SaveWork: Cannot open %s for writing\n", workFile.c_str());
    return;
  }

  Lock();

  // Wait that all threads blocks before saving works
  saveRequest = true;
  while (!isWaiting(threads))
    Timer::SleepMillis(50);

  // Save global param
  fwrite(initialSeed.c_str(), 64, 1, f);
  fwrite(&colSize, sizeof(uint32_t), 1, f);
  fwrite(&dpSize, sizeof(uint32_t), 1, f);
  fwrite(&extraPoints, sizeof(bool), 1, f);
  fwrite(&totalCount, sizeof(uint64_t), 1, f);
  fwrite(&totalTime, sizeof(double), 1, f);

  // Save random walks
  uint64_t totalWalk = 0;
  for (int i = 0; i < nbThread; i++)
    totalWalk += threads[i].nbWalk;

  fwrite(&totalWalk, sizeof(uint64_t), 1, f);
  for (int i = 0; i < nbThread; i++) {
    for (uint64_t n = 0; n < threads[i].nbWalk; n++) {
      fwrite(&threads[i].x[n], sizeof(hash160_t), 1, f);
      fwrite(&threads[i].y[n], sizeof(hash160_t), 1, f);
    }
  }

  // Save hash table
  hashTable.SaveTable(f);
  fclose(f);

  // Unblock threads
  saveRequest = false;
  Unlock();

  char *ctimeBuff;
  time_t now = time(NULL);
  ctimeBuff = ctime(&now);
  printf("\nWork saved: %s", ctimeBuff);
}

// ----------------------------------------------------------------------------

void BTCCollider::FetchWalks(hash160_t *x, hash160_t *y, uint64_t nbWalk) {
  for (int n = 0; n < nbWalk; n++) {
    if (fetchedWalk < nbLoadedWalk) {
      // Fecth a loaded walk
      x[n] = loadedX[fetchedWalk];
      y[n] = loadedY[fetchedWalk];
      fetchedWalk++;
    } else {
      // New random walk
      Rand(&seed, &x[n]);
      y[n] = x[n];
    }
  }
}

// ----------------------------------------------------------------------------

void BTCCollider::Check(std::vector<int> gpuId, std::vector<int> gridSize) {
  if(initDPSize<0)
    initDPSize = colSize/3;
  SetDP(initDPSize);

  // Check Int lib
  //Int::Check();

  // Check SSE and CPU group
  hash160_t *x = new hash160_t[CPU_GRP_SIZE];
  hash160_t *xc = new hash160_t[CPU_GRP_SIZE];
  for (int i = 0; i < CPU_GRP_SIZE; i++) {
    Rand(&seed, &x[i]);
    xc[i] = x[i];
  }

  for (int i = 0; i < CPU_GRP_SIZE; i++) {
    xc[i] = F(xc[i]);
  }

  IntGroup *grp = new IntGroup(CPU_GRP_SIZE);
  Point *pts = new Point[CPU_GRP_SIZE];
  Int *dInv = new Int[CPU_GRP_SIZE];
  FGroup(grp, pts, dInv, x);

  bool ok = true;
  int i = 0;
  while (ok && i < CPU_GRP_SIZE) {
    ok = (hashTable.compareHash(&x[i], &xc[i]) == 0);
    if (ok) i++;
  }

  if (ok) {
    printf("CPU Group OK!\n");
  } else {
    printf("CPU Group Not OK at %d!\n", i);
    printf("Hg=%s\n", GetHex(x[i]).c_str());
    printf("Hc=%s\n", GetHex(xc[i]).c_str());
  }

#ifdef WITHGPU
  // Check gpu
  if (useGpu) {
    printf("GPU allocate memory:");
    int x = gridSize[0];
    int y = gridSize[1];
    if (!GPUEngine::GetGridSize(gpuId[0], &x, &y)) {
      return;
    }
    GPUEngine h(x,y, gpuId[0], 65536);
    printf(" done\n");
    printf("GPU: %s\n", h.deviceName.c_str());
    printf("GPU: %.1f MB\n", h.GetMemory()/1048576.0);
    int nbH = h.GetNbThread() * GPU_GRP_SIZE;
    hash160_t *iHash = (hash160_t *)malloc(nbH *sizeof(hash160_t));
    for(int i=0;i<nbH;i++)
      Rand(&seed, &iHash[i]);
    h.SetExtraPoint(extraPoints);
    h.SetMasks(colMask,dMask,nbFull);
    printf("GPU SetKeys:");
    h.SetKeys(pub);
    printf(" done\n");
    printf("GPU SetStartingHashes:");
    if (!h.SetStartingHashes((uint64_t *)iHash, (uint64_t *)iHash)) {
      printf(" failed !");
      return;
    }
    printf(" done\n");
    HashTable *h1 = new HashTable();
    HashTable *h2 = new HashTable();
    h1->SetParam(colSize,nbFull,colMask);
    h2->SetParam(colSize, nbFull, colMask);
 
    vector<ITEM> hashFound;
    h.Launch(hashFound);
    for (int i = 0; i < (int)hashFound.size(); i++)
      h1->AddHash((hash160_t *)(hashFound[i].h1), (hash160_t *)(hashFound[i].h2));
    printf("GPU found %d items\n", h1->GetNbItem());
    int nb
    int nbCPU = 0;
    hash160_t *h = (hash160_t *)malloc(nbH * sizeof(hash160_t));
    hash160_t *hc = (hash160_t *)malloc(nbH * sizeof(hash160_t));
    for (int i = 0; i < nbH; i++) {
      h[i] = iHash[i];
      hc[i] = iHash[i];
    }
    for (int i = 0; i < nbH; i++) {
      hc[i] = F(hc[i]);
      if (hashTable.compareHash(&h[i], &hc[i]) == 0) {
        h2->AddHash(&h[i], &hc[i]);
        nbCPU++;
      }
    }
    printf("CPU found %d items\n", nbCPU);

    if (nbCPU != h1->GetNbItem()) {
      printf("CPU/GPU differs !\n");
    } else {
      int nb = h1->GetNbItem();
      bool ok = true;
      int i = 0;
      while (ok && i < nb) {
        vector<Int> hashesGPU = h1->GetItemsInt(i);
        vector<Int> hashesCPU = h2->GetItemsInt(i);
        ok = (hashesGPU[0].IsEqual(&hashesCPU[0]) && hashesGPU[1].IsEqual(&hashesCPU[1]));
        i++;
      }
      if (!ok) {
        printf("CPU/GPU differs !\n");
      } else {
        printf("CPU/GPU OK\n");
      }
    }

    delete h1;
    delete h2;
    free(h);
    free(hc);
    free(iHash);
  }
#endif

  delete[] x;
  delete[] xc;
  delete grp;
  delete[] pts;
  delete[] dInv;
}

// ----------------------------------------------------------------------------

void BTCCollider::InitKey(TH_PARAM *p) {
  Int k;
  Point key;
  int i = p->threadId;
  for (int j = 0; j < 65536; j++) {
    Rand(&p->localSeed, &k);
    key = secp->ComputePublicKey(&k);
    pub[i * 131072 + j * 2] = key.x;
    pub[i * 131072 + j * 2 + 1] = key.y;
    if ((j % 8192) == 0) {
      printf(".");
#ifndef WIN64
      fflush(stdout);
#endif
    }
  }
}

// ----------------------------------------------------------------------------

void BTCCollider::SetDP(int size) {
  dpSize = size;
  dMask = 0xFFFFFFFFFFFFFFFFULL;
  dMask = dMask >> (64 - dpSize);
}

// ----------------------------------------------------------------------------

void BTCCollider::Run(int nbThread) {
  double t0;
  double t1;
  uint64_t count;
  uint64_t lastCount;
  uint64_t nbFoundCollision;
  vector<ITEM> hashFound;
  vector<ITEM> hashFound2;

  TH_PARAM *params = (TH_PARAM *)malloc(nbThread * sizeof(TH_PARAM));
  THREAD_HANDLE *thHandles = (THREAD_HANDLE *)malloc(nbThread * sizeof(THREAD_HANDLE));

  memset(params, 0, nbThread * sizeof(TH_PARAM));

  printf("Number of CPU thread: %d\n", nbThread);

  // Set starting parameters
  for (int i = 0; i < nbThread; i++) {
    params[i].obj = this;
    params[i].threadId = i;
    params[i].isRunning = true;
    params[i].nbit = 0;
    params[i].hStart = (hash160_t *)malloc(CPU_GRP_SIZE * sizeof(hash160_t));
    params[i].hStop = (hash160_t *)malloc(CPU_GRP_SIZE * sizeof(hash160_t));
    params[i].x = (hash160_t *)malloc(CPU_GRP_SIZE * sizeof(hash160_t));
    params[i].y = (hash160_t *)malloc(CPU_GRP_SIZE * sizeof(hash160_t));
    params[i].nbWalk = 0;
    params[i].maxWalk = CPU_GRP_SIZE;
    FetchWalks(params[i].x, params[i].y, params[i].maxWalk);
    params[i].nbWalk = params[i].maxWalk;
  }

  // Launch CPU threads
  for (int i = 0; i < nbThread; i++)
    thHandles[i] = LaunchThread(_FindCollision, params + i);

#ifdef WITHGPU
  // Launch GPU
  if (useGpu) {
    int x = gridSize[0];
    int y = gridSize[1];
    if (!GPUEngine::GetGridSize(gpuId[0], &x, &y)) {
      return;
    }
    GPUEngine *g = new GPUEngine(x, y, gpuId[0], 65536);
    g->SetExtraPoint(extraPoints);
    g->SetMasks(colMask, dMask, nbFull);
    g->SetKeys(pub);
    nbGPUThread = g->GetNbThread();
    printf("GPU: %s\n", g->deviceName.c_str());
    printf("GPU: %.1f MB\n", g->GetMemory() / 1048576.0);
    printf("GPU: %d threads\n", nbGPUThread);
  }
#endif

  // Key rate smoothing filter
  #define FILTER_SIZE 8
  double lastkeyRate[FILTER_SIZE];
  double lastGpukeyRate[FILTER_SIZE];
  uint32_t filterPos = 0;

  double keyRate = 0.0;
  double gpuKeyRate = 0.0;
  char timeStr[256];

  // Wait that all threads have started
  while (!isReady(params)) {
    Timer::SleepMillis(500);
  }

  // Reset timer
  Timer::Init();
  t0 = Timer::get_tick();
  t1 = t0;
  lastCount = 0;
  nbFoundCollision = 0;

  while (isAlive(params)) {
    int delay = 2000;
    while (isAlive(params) && delay > 0) {
      Timer::SleepMillis(500);
      delay -= 500;
    }

    count = getGlobalCount();
    uint64_t totalCount = count + offsetCount;
    double t2 = Timer::get_tick();
    double keyRate = (double)(count - lastCount) / (t2 - t1);
    double totalKeyRate = (double)(totalCount) / (t2 - t0 + offsetTime);

    memcpy(lastkeyRate, lastkeyRate + 1, (FILTER_SIZE - 1) * sizeof(double));
    lastkeyRate[FILTER_SIZE - 1] = keyRate;
    if (filterPos == FILTER_SIZE) {
      double avgKeyRate = 0.0;
      for (int i = 0; i < FILTER_SIZE; i++)
        avgKeyRate += lastkeyRate[i];
      avgKeyRate /= (double)(FILTER_SIZE);
      char *crs = "\n";
      if (rekey > 0) {
        crs = "";
      }
      printf("%s[%.2f Mkey/s][%.2f Mkey/s][2^%.2f][%s][%s]%s",
             crs,
             avgKeyRate / 1000000.0,
             totalKeyRate / 1000000.0,
             log2((double)totalCount),
             GetTimeStr(t2 - t0 + offsetTime).c_str(),
             hashTable.GetSizeInfo().c_str(),
             crs);
    } else {
      filterPos++;
    }

    if (saveWorkPeriod > 0) {
      if ((count % saveWorkPeriod) == 0) {
        SaveWork(totalCount, t2 - t0 + offsetTime, params, nbThread);
      }
    }

    lastCount = count;
    t1 = t2;
  }

  free(params);
  free(thHandles);
}

// ----------------------------------------------------------------------------

void BTCCollider::FindCollision(TH_PARAM *p) {
  // Global init
  TH_PARAM *params = (TH_PARAM *)p;
  int nbThread = params->nbThread;
  int thl = params->threadId;
  counters[thl] = 0;

  // CPU Thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE);
  Point *pts = new Point[CPU_GRP_SIZE];
  Int *dInv = new Int[CPU_GRP_SIZE];

  // Using Affine coord
  Int dy;
  Int dx;
  Int rx;
  Int ry;
  Int _s;
  Int _p;

  uint64_t nbStep = 0;
  hash160_t *x = params->x;
  hash160_t *y = params->y;
  hash160_t *hStart = params->hStart;
  hash160_t *hStop = params->hStop;

  while (params->isRunning) {
    // Random walk
    for (int j = 0; j < 1024 * 16 && params->isRunning; j++) {
      for (int i = 0; i < params->nbWalk; i++) {
        hStart[i] = x[i];
        hStop[i] = y[i];
      }

      FGroup(grp, pts, dInv, x);

      for (int i = 0; i < params->nbWalk; i++) {
        x[i] = F(x[i]);
        y[i] = F(F(y[i]));
      }

      nbStep++;
      counters[thl] += params->nbWalk;
    }

    // Check collision
    if (params->isRunning) {
      for (int i = 0; i < params->nbWalk; i++) {
        if (hashTable.compareHash(&x[i], &y[i]) == 0) {
          // Collision
          Lock();
          vector<Int> col = hashTable.GetCollision(&x[i], &y[i]);
          if (col.size() > 0) {
            Int privKey = col[0];
            Int privKey2 = col[1];
            Point publicKey = secp->ComputePublicKey(&privKey);
            Point publicKey2 = secp->ComputePublicKey(&privKey2);
            printf("\nCollision found:\n");
            printf("PrivKey1: %s\n", privKey.GetBase16().c_str());
            printf("PrivKey2: %s\n", privKey2.GetBase16().c_str());
            printf("PubKey1: %s\n", publicKey.GetBase16Compressed().c_str());
            printf("PubKey2: %s\n", publicKey2.GetBase16Compressed().c_str());
            FILE *f = fopen(outputFile.c_str(), "a");
            if (f) {
              fprintf(f, "PrivKey1: %s\n", privKey.GetBase16().c_str());
              fprintf(f, "PrivKey2: %s\n", privKey2.GetBase16().c_str());
              fprintf(f, "PubKey1: %s\n", publicKey.GetBase16Compressed().c_str());
              fprintf(f, "PubKey2: %s\n", publicKey2.GetBase16Compressed().c_str());
              fclose(f);
            }
          }
          Unlock();
        }
      }
    }

    // Save hash
    if (params->isRunning) {
      Lock();
      for (int i = 0; i < params->nbWalk; i++)
        hashTable.AddHash(&hStart[i], &hStop[i]);
      Unlock();
    }

    // Wait for save request
    if (params->isRunning) {
      while (saveRequest) {
        Timer::SleepMillis(50);
      }
    }
  }

  delete grp;
  delete[] pts;
  delete[] dInv;
}

// ----------------------------------------------------------------------------

hash160_t BTCCollider::F(hash160_t x) {
  uint64_t *x64 = (uint64_t *)&x;
  uint64_t *y64 = (uint64_t *)&x;
  uint32_t *x32 = (uint32_t *)&x;
  uint32_t *y32 = (uint32_t *)&x;

  // Derived from pairgen (https://github.com/basil00/pairgen.git)
  // Given a hash160 H comprised on {h0, h1, .., h9} 16-bit parts, then H is
  // mapped to a public key P as follows:
  //     P = pub[0][h0] + pub[1][h1] + .. + pub[9][h9]
  // The calculation is truncated according to the length n.  The corresponding
  // private key P' is:
  //     P' = priv[0]+h0 + priv[1]+h1*2^16 + .. + priv[9]+h9*2^144
  // Each base private key is chosen randomly and computed in advanced. 

  Int px;
  Int py;
  Int dx;
  Int dy;
  Int rx;
  Int ry;
  Int _s;
  Int _p;

  px.Set32Bytes((unsigned char *)&pub[0 * 131072 + x32[0] * 2]);
  py.Set32Bytes((unsigned char *)&pub[0 * 131072 + x32[0] * 2 + 1]);

  for (int i = 1; i < 5; i++) {
    dx.Set32Bytes((unsigned char *)&pub[i * 131072 + x32[i] * 2]);
    dy.Set32Bytes((unsigned char *)&pub[i * 131072 + x32[i] * 2 + 1]);
    secp->AddDirect(px, py, dx, dy);
  }

  // P = lambda1*(P + beta1.G)
  secp->AddDirect(px, py, beta1);
  secp->MulDirect(px, py, lambda1);

  // P = P + beta2.G
  secp->AddDirect(px, py, beta2);

  // P = lambda2*P
  secp->MulDirect(px, py, lambda2);

  px.Get32Bytes((unsigned char *)y64);
  py.Get32Bytes((unsigned char *)(y64 + 1));

  return x;
}

// ----------------------------------------------------------------------------

void BTCCollider::FGroup(IntGroup *grp, Point *pts, Int *dInv, hash160_t *x) {
  uint64_t *x64 = (uint64_t *)x;
  uint32_t *x32 = (uint32_t *)x;

  Int px[CPU_GRP_SIZE];
  Int py[CPU_GRP_SIZE];
  Int dx[CPU_GRP_SIZE];
  Int dy[CPU_GRP_SIZE];

  for (int i = 0; i < CPU_GRP_SIZE; i++) {
    px[i].Set32Bytes((unsigned char *)&pub[0 * 131072 + x32[i * 5] * 2]);
    py[i].Set32Bytes((unsigned char *)&pub[0 * 131072 + x32[i * 5] * 2 + 1]);
  }

  for (int j = 1; j < 5; j++) {
    for (int i = 0; i < CPU_GRP_SIZE; i++) {
      dx[i].Set32Bytes((unsigned char *)&pub[j * 131072 + x32[i * 5 + j] * 2]);
      dy[i].Set32Bytes((unsigned char *)&pub[j * 131072 + x32[i * 5 + j] * 2 + 1]);
    }
    grp->AddDirect(px, py, dx, dy);
  }

  // P = lambda1*(P + beta1.G)
  grp->AddDirect(px, py, beta1);
  grp->MulDirect(px, py, lambda1);

  // P = P + beta2.G
  grp->AddDirect(px, py, beta2);

  // P = lambda2*P
  grp->MulDirect(px, py, lambda2);

  for (int i = 0; i < CPU_GRP_SIZE; i++) {
    px[i].Get32Bytes((unsigned char *)(x64 + i * 5));
    py[i].Get32Bytes((unsigned char *)(x64 + i * 5 + 1));
  }
}

// ----------------------------------------------------------------------------

string BTCCollider::GetHex(hash160_t h) {
  char tmp[41];
  for (int i = 0; i < 20; i++)
    sprintf(tmp + 2 * i, "%02X", (uint8_t)(h.i8[i]));
  return string(tmp);
}

// ----------------------------------------------------------------------------

bool BTCCollider::isAlive(TH_PARAM *params) {
  for (int i = 0; i < nbThread; i++)
    if (params[i].isRunning)
      return true;
  return false;
}

// ----------------------------------------------------------------------------

bool BTCCollider::isReady(TH_PARAM *params) {
  for (int i = 0; i < nbThread; i++)
    if (!params[i].isReady)
      return false;
  return true;
}

// ----------------------------------------------------------------------------

bool BTCCollider::isWaiting(TH_PARAM *params) {
  for (int i = 0; i < nbThread; i++)
    if (!params[i].isWaiting)
      return false;
  return true;
}

// ----------------------------------------------------------------------------

uint64_t BTCCollider::getGlobalCount() {
  uint64_t count = 0;
  for (int i = 0; i < nbThread; i++)
    count += counters[i];
  return count;
}

// ----------------------------------------------------------------------------

void BTCCollider::SetupRanges(uint32_t totalThreads) {
  // Set ranges
  Int threads;
  threads.SetInt32(totalThreads);
  rangeDiff = secp->order.Div(&threads);
  rangeStart.Set(&secp->order);
  rangeStart.Sub(&rangeDiff);
  rangeEnd.Set(&secp->order);
  rangeStart.Add(1);
}

// ----------------------------------------------------------------------------

void BTCCollider::getCPUStartingKey(int thId, Int *key) {
  Int k(&rangeStart);
  Int r;
  r.Rand(&rangeStart, &rangeEnd);
  k.Add(&r);
  key->Set(&k);
}

// ----------------------------------------------------------------------------

void BTCCollider::SetKeys(Point *p) {
  // Set starting public key for each thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE);
  Point *pts = new Point[CPU_GRP_SIZE];
  Int *dx = new Int[CPU_GRP_SIZE];
  Int *dy = new Int[CPU_GRP_SIZE];

  for (int i = 0; i < nbThread; i++) {
    Int pk;
    getCPUStartingKey(i, &pk);
    Point p = secp->ComputePublicKey(&pk);
    startPubKey[i * 2] = p.x;
    startPubKey[i * 2 + 1] = p.y;
  }

  delete grp;
  delete[] pts;
  delete[] dx;
  delete[] dy;
}

// ----------------------------------------------------------------------------

string BTCCollider::GetTimeStr(double s) {
  char tmp[256];
  int d = (int)(s / 86400.0);
  s -= d * 86400;
  int h = (int)(s / 3600.0);
  s -= h * 3600;
  int m = (int)(s / 60.0);
  s -= m * 60;
  if (d > 0) sprintf(tmp, "%d Day %02d:%02d:%02.0f", d, h, m, s);
  else sprintf(tmp, "%02d:%02d:%02.0f", h, m, s);
  return string(tmp);
}

// ----------------------------------------------------------------------------

void BTCCollider::Rand(Int *seed, hash160_t *h) {
  uint8_t tmp[32];
  sha256((uint8_t *)seed->bits64, 32, tmp);
  memcpy(h->i8, tmp, 20);
  seed->SetInt32(0);
  seed->Set32Bytes(tmp);
}

// ----------------------------------------------------------------------------

void BTCCollider::Rand(Int *seed, Int *i) {
  uint8_t tmp[32];
  sha256((uint8_t *)seed->bits64, 32, tmp);
  i->SetInt32(0);
  i->Set32Bytes(tmp);
  seed->SetInt32(0);
  seed->Set32Bytes(tmp);
}

// ----------------------------------------------------------------------------

#ifdef WIN64
DWORD WINAPI _FindCollision(LPVOID lpParam) {
#else
void *_FindCollision(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->FindCollision(p);
  p->isRunning = false;
  return 0;
}

// ----------------------------------------------------------------------------

BTCCollider::~BTCCollider() {
  if (pub) free(pub);
  if (loadedX) free(loadedX);
  if (loadedY) free(loadedY);
}
