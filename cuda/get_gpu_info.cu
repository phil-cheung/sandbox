/* This code determines the spec of CUDA devices */
#include <stdio.h>

int getSPcores(cudaDeviceProp devProp)
{
    int cores = 0;
    int mp = devProp.multiProcessorCount;
    switch (devProp.major){
     case 2: // Fermi
      if (devProp.minor == 1) cores = mp * 48;
      else cores = mp * 32;
      break;
     case 3: // Kepler
      cores = mp * 192;
      break;
     case 5: // Maxwell
      cores = mp * 128;
      break;
     case 6: // Pascal
      if (devProp.minor == 1) cores = mp * 128;
      else if (devProp.minor == 0) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     default:
      printf("Unknown device type\n");
      break;
      }
    return cores;
}

void get_device_info() {
    int dev_count;
    cudaGetDeviceCount( &dev_count);
    printf("Num of CUDA devices: %i\n", dev_count);

    cudaDeviceProp dev_prop;
    int i;
    for (i=0; i < dev_count; i++) {
        cudaGetDeviceProperties( &dev_prop, i);

        /* Device architecture limits */
        printf("\nDEVICE LIMITS\n");
        printf("Max Threads per block: %i\n", dev_prop.maxThreadsPerBlock);
        printf("Num of SM: %i\n", dev_prop.multiProcessorCount);
        printf("Num of cores: %i\n", getSPcores(dev_prop));
        printf("Clock speed per core (kHz): %i\n", dev_prop.clockRate);
        printf("Max threads per SM: %i\n", dev_prop.maxThreadsPerMultiProcessor);
        printf("Warp size: %i\n", dev_prop.warpSize);

        /* Memory limits */
        printf("\nMEMORY LIMITS\n");
        printf("Total Global Mem (bytes): %zd\n", dev_prop.totalGlobalMem);
        printf("Total num of registers per Multiprocessor: %i\n", dev_prop.regsPerMultiprocessor);
        printf("Total num of registers per block: %i\n", dev_prop.regsPerBlock);
        printf("Shared Mem per MultiProcessor (bytes): %zd\n", dev_prop.sharedMemPerMultiprocessor);
        printf("Shared Mem per Block (bytes): %zd\n", dev_prop.sharedMemPerBlock);
        printf("Total Const Mem (bytes): %zd\n", dev_prop.totalConstMem);
    }
}


int main() {
    get_device_info();

    return EXIT_SUCCESS;
}
