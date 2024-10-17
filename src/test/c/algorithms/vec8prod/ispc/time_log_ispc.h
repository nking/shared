#ifndef TIME_LOG_ISPC
#define TIME_LOG_ISPC

// process the compile options to create timing logs
#if defined(TIME_D) && defined(TIME_THR)
   #error "only 1 of TIME_D, TIME_THR, or TIME_TOT can be set"
#elif defined(TIME_D) && defined(TIME_TOT)
   #error "only 1 of TIME_D, TIME_THR, or TIME_TOT can be set"
#elif defined(TIME_THR) && defined(TIME_TOT)
   #error "only 1 of TIME_D, TIME_THR, or TIME_TOT can be set"
#endif

// should only be 1 log flag set if any at this point
// TIME_D logs memory data access, such as loading vectors
// TIME_THR logs the individual thread methods
// TIME_TOT logs the entire method
#if defined(TIME_D)
    #define INIT_TIME_TITLE(s) print("TITLE %\n", #s);
    #define INIT_TIME() unsigned uniform int64 start, stop;
    #define START_TOT_TIME()
    #define STOP_TOT_TIME(s)
    #define START_D_TIME() start = clock();
    #define STOP_D_LOAD_TIME() {stop = clock(); print("load data %\n", (stop-start));}
    #define STOP_D_STORE_TIME() {stop = clock(); print("store data %\n", (stop-start));}
    #define START_THR_TIME()
    #define STOP_THR_TIME(s)
#elif defined(TIME_THR)
    #define INIT_TIME_TITLE(s) print("TITLE %\n", #s);
    #define INIT_TIME() unsigned uniform int64 start, stop;
    #define START_TOT_TIME()
    #define STOP_TOT_TIME(s)
    #define START_D_TIME()
    #define STOP_D_TIME(s)
    #define START_THR_TIME() start = clock();
    #define STOP_THR_TIME() {stop = clock(); print("thr thr %\n", #s, (stop-start));}
#elif defined(TIME_TOT)
    #define INIT_TIME_TITLE(s) print("TITLE %\n", #s);
    #define INIT_TIME() unsigned uniform int64 start, stop;
    #define START_TOT_TIME() start = clock();
    #define STOP_TOT_TIME(s) {stop = clock(); print("% tot %\n", #s, (stop-start));}
    #define START_D_TIME()
    #define STOP_D_TIME(s)
    #define START_THR_TIME()
    #define STOP_THR_TIME(s)
#else
   #define INIT_TIME_TITLE(s)
   #define INIT_TIME()
   #define START_TOT_TIME()
   #define STOP_TOT_TIME(s)
   #define START_D_TIME()
   #define STOP_D_TIME(s)
   #define START_THR_TIME()
   #define STOP_THR_TIME(s)
#endif

//LOG(format, ...) fprintf(stdout, format, ##__VA_ARGS__)

#endif