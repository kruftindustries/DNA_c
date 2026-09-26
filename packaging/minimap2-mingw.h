// Put ahead of every minimap2 source when building with MinGW-w64
// (CPPFLAGS="-include packaging/minimap2-mingw.h"). mingw-w64's ftell and
// fseek work on 32-bit offsets; minimap2 tells whether an index file has
// further parts by comparing ftell() with the index size (index.c,
// mm_idx_reader_eof), which never matched an 8.5 GB whole-genome .mmi, so
// it treated the index as multi-part and wrote SAM without @SQ lines.
// stdio.h is included first so its own prototypes are left alone.
#pragma once
#ifdef _WIN32
#include <stdio.h>
#define ftell _ftelli64
#define fseek _fseeki64
#endif
