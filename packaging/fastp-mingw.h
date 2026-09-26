// Put ahead of every fastp source when building with MinGW-w64
// (g++ -include packaging/fastp-mingw.h): the POSIX calls mingw-w64 lacks
// or narrows. fastp's own sources are otherwise portable.
//
//  - pwrite: fastp's writer threads write finished blocks concurrently at
//    computed offsets, so it has to be a real positional write, not a
//    seek + write. fastp_mingw_pwrite (fastp-mingw.cpp) is WriteFile with
//    an OVERLAPPED offset, which is Windows' pwrite.
//  - ftell/fseek return and take `long`, 32 bits on Windows; the 64-bit
//    variants keep a >2 GB FASTQ readable.
#pragma once
#ifdef _WIN32
#include <sys/types.h>
#include <stdio.h>
#ifdef __cplusplus
extern "C"
#endif
ssize_t fastp_mingw_pwrite(int fd, const void *buf, size_t count, long long offset);
#define pwrite fastp_mingw_pwrite
#define ftell _ftelli64
#define fseek _fseeki64
#endif
