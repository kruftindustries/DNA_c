#!/usr/bin/env bash
# Build the sequencing tools the Windows release ships in its tools\ folder:
# samtools, bcftools (both on htslib) and minimap2, from their release
# tarballs, in an MSYS2 MinGW64 shell. This is the environment their own
# projects test Windows builds in (htslib, samtools and bcftools each run a
# mingw64 job upstream); minimap2 has no Windows CI but carries a WIN32
# branch and builds as is. fastp (read QC, optional to the pipeline) has no
# Windows port of its own; it is built best-effort from MSYS2's isa-l,
# libdeflate and highway packages with fastp-mingw.{h,cpp} supplying pwrite
# and 64-bit ftell/fseek, and left out with a warning if that fails.
#
# Packages (pacman -S): mingw-w64-x86_64-toolchain mingw-w64-x86_64-zlib
#   mingw-w64-x86_64-bzip2 mingw-w64-x86_64-xz mingw-w64-x86_64-libdeflate
#   mingw-w64-x86_64-curl-winssl mingw-w64-x86_64-isa-l mingw-w64-x86_64-highway
#   make autoconf automake
# curl-winssl rather than curl: htslib reaches https:// URLs (the Ensembl
# FASTA, the 1000 Genomes VCFs) through libcurl, and the Schannel build uses
# the Windows certificate store, so no CA bundle has to be shipped.
#
#   packaging/windows-tools.sh [outdir]    -> packaging/tools/windows/
#
# The output folder holds the three .exe files, the MinGW DLLs they load
# (found with ldd) and the licences. windows-package.bat copies it into the
# zip as tools\; ToolLocator looks there before PATH.
set -euo pipefail

HTS_VERSION=${HTS_VERSION:-1.24}
MINIMAP2_VERSION=${MINIMAP2_VERSION:-2.31}
FASTP_VERSION=${FASTP_VERSION:-1.3.7}

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/tools/windows}
work=${WORK:-$here/tools/build-windows}
mkdir -p "$out" "$work"
jobs=$(nproc 2>/dev/null || echo 4)

fetch() {  # url [tar flags]
    local f=$work/$(basename "$1")
    [ -s "$f" ] || curl -fsSL --retry 5 --retry-delay 5 -o "$f" "$1"
    tar -x${2:-j}f "$f" -C "$work"
}

echo "== htslib $HTS_VERSION"
fetch "https://github.com/samtools/htslib/releases/download/$HTS_VERSION/htslib-$HTS_VERSION.tar.bz2"
( cd "$work/htslib-$HTS_VERSION"
  ./configure --enable-libcurl --disable-gcs --disable-s3 --disable-plugins
  make -j"$jobs" lib-static )

echo "== samtools $HTS_VERSION"
fetch "https://github.com/samtools/samtools/releases/download/$HTS_VERSION/samtools-$HTS_VERSION.tar.bz2"
( cd "$work/samtools-$HTS_VERSION"
  ./configure --with-htslib="$work/htslib-$HTS_VERSION" --without-curses --disable-ref-cache
  make -j"$jobs" samtools
  [ -f samtools.exe ] || mv samtools samtools.exe )

echo "== bcftools $HTS_VERSION"
fetch "https://github.com/samtools/bcftools/releases/download/$HTS_VERSION/bcftools-$HTS_VERSION.tar.bz2"
( cd "$work/bcftools-$HTS_VERSION"
  ./configure --with-htslib="$work/htslib-$HTS_VERSION" --disable-bcftools-plugins
  make -j"$jobs" bcftools
  [ -f bcftools.exe ] || mv bcftools bcftools.exe )

echo "== minimap2 $MINIMAP2_VERSION"
fetch "https://github.com/lh3/minimap2/releases/download/v$MINIMAP2_VERSION/minimap2-$MINIMAP2_VERSION.tar.bz2"
( cd "$work/minimap2-$MINIMAP2_VERSION"
  # mingw-w64's ftell/fseek are 32-bit. minimap2 tells whether an index has
  # more parts by comparing ftell() with the index size (index.c,
  # mm_idx_reader_eof); with an 8.5 GB whole-genome .mmi that never
  # matched, minimap2 assumed a multi-part index and wrote SAM without @SQ
  # lines; minimap2-mingw.h maps them to the 64-bit variants (after
  # stdio.h, so its prototypes are untouched). HAVE_KALLOC is the
  # Makefile's own CPPFLAGS, kept.
  make -j"$jobs" CPPFLAGS="-DHAVE_KALLOC -include $here/minimap2-mingw.h" minimap2
  [ -f minimap2.exe ] || mv minimap2 minimap2.exe )

# The makefiles name their targets without .exe; MinGW's gcc adds the
# suffix to the file it writes, hence the target names above and the
# renames in case a toolchain does not.
echo "== fastp $FASTP_VERSION (best effort)"
fastp_built=no
if fetch "https://github.com/OpenGene/fastp/archive/refs/tags/v$FASTP_VERSION.tar.gz" z \
   && ( cd "$work/fastp-$FASTP_VERSION" \
        && cp "$here/fastp-mingw.cpp" src/fastp_mingw_compat.cpp \
        && make -j"$jobs" CXX="g++ -include $here/fastp-mingw.h" fastp \
        && { [ -f fastp.exe ] || mv fastp fastp.exe; } ); then
    fastp_built=yes
else
    echo "::warning::fastp $FASTP_VERSION did not build under MSYS2; the Windows package ships without it (read QC is skipped)"
fi

echo "== collecting into $out"
rm -rf "$out"
mkdir -p "$out/licenses"
cp "$work/samtools-$HTS_VERSION/samtools.exe" "$work/bcftools-$HTS_VERSION/bcftools.exe" \
   "$work/minimap2-$MINIMAP2_VERSION/minimap2.exe" "$out/"
cp "$work/htslib-$HTS_VERSION/LICENSE" "$out/licenses/htslib-LICENSE"
cp "$work/samtools-$HTS_VERSION/LICENSE" "$out/licenses/samtools-LICENSE"
cp "$work/bcftools-$HTS_VERSION/LICENSE" "$out/licenses/bcftools-LICENSE"
cp "$work/minimap2-$MINIMAP2_VERSION/LICENSE.txt" "$out/licenses/minimap2-LICENSE"
if [ "$fastp_built" = yes ]; then
    cp "$work/fastp-$FASTP_VERSION/fastp.exe" "$out/"
    cp "$work/fastp-$FASTP_VERSION/LICENSE" "$out/licenses/fastp-LICENSE"
fi
cat > "$out/README.txt" <<TXT
Sequencing tools for the Genetic Health FASTQ pipeline and its 1000 Genomes
test sample, built for Windows x64 by packaging/windows-tools.sh:

  samtools $HTS_VERSION, bcftools $HTS_VERSION (htslib $HTS_VERSION)   https://www.htslib.org
  minimap2 $MINIMAP2_VERSION                                          https://github.com/lh3/minimap2
$([ "$fastp_built" = yes ] && echo "  fastp $FASTP_VERSION                                            https://github.com/OpenGene/fastp" \
                            || echo "  (fastp did not build for this release; the pipeline skips read QC without it)")

Each is MIT-licensed by its authors (see licenses/). The DLLs are the MinGW-w64
runtime libraries they load (zlib, bzip2, xz, libdeflate, libcurl, winpthreads, ...).
TXT

# The DLLs each executable loads from the MinGW prefix (system DLLs are
# left to Windows). ldd understands PE files under MSYS2.
for exe in "$out"/*.exe; do
    ldd "$exe" | awk '/mingw64/ { print $3 }' | while read -r dll; do
        [ -f "$out/$(basename "$dll")" ] || cp "$dll" "$out/"
    done
done
ls -l "$out"
echo "== done: $(ls "$out"/*.exe | wc -l) executables, $(ls "$out"/*.dll | wc -l) DLLs"
