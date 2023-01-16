#include "classify.h"
#include <omp.h>
#include <math.h>

Data classify(Data &D, const Ranges &R, unsigned int numt)
{ // Classify each item in D into intervals (given by R). Finally, produce in D2 data sorted by interval
   assert(numt < MAXTHREADS);
   Counter counts[R.num()]; // I need on counter per interval. Each counter can keep pre-thread subcount.
   #pragma omp parallel num_threads(numt)
   {
      int tid = omp_get_thread_num(); // I am thread number tid
      int len = D.ndata/numt;
      int start = tid*len;
      int end = tid*len+len;
      if (end+len> D.ndata) end = D.ndata;
      for(int i=start; i<end; ++i) { // Threads together share-loop through all of Data
         int v = D.data[i].value = R.range(D.data[i].key);// For each data, find the interval of data's key,
							  // and store the interval id in value. D is changed.
         counts[v].increase(tid); // Found one key in interval v
      }
   }

   // Accumulate all sub-counts (in each interval;'s counter) into rangecount
   unsigned int *rangecount = new unsigned int[R.num()];
   for(int r=0; r<R.num(); r++) { // For all intervals
      rangecount[r] = 0;
      for(int t=0; t<numt; t++) // For all threads
         rangecount[r] += counts[r].get(t);
      // std::cout << rangecount[r] << " elements in Range " << r << "\n"; // Debugging statement
   }

   // Compute prefx sum on rangecount.
   for(int i=1; i<R.num(); i++) {
      rangecount[i] += rangecount[i-1];
   }
   for(int i=R.num()-1; i>0; i--) {
      rangecount[i] = rangecount[i-1];
   }
   rangecount[0] = 0;

   // Now rangecount[i] has the number of elements in intervals before the ith interval.

   Data D2 = Data(D.ndata); // Make a copy
   int cap = ceil((float)R.num()/numt);           // This is the max number of ranges that will be handled by a thread.
   #pragma omp parallel num_threads(numt)
   {
      int tid = omp_get_thread_num();
      int start = tid*cap;
      int end = tid*cap+cap;
      for (int i=0;i<D.ndata; i++){
         int r = D.data[i].value;
         if ( start<=r && r< end ) D2.data[rangecount[r]++] = D.data[i];
      }
   }
   return D2;
}