#!/usr/bin/env bash
# Build the sequencing tools the Windows release ships in its tools\ folder:
# samtools, bcftools (both on htslib) and minimap2, from their release
# tarballs, in an MSYS2 MinGW64 shell. This is the environment their own
# projects test Windows builds in (htslib, samtools and bcftools each run a
# mingw64 job upstream); minimap2 has no Windows CI but carries a WIN32
# branch and builds as is. fastp is not built: it needs ISA-L and Highway,
# has no Windows port, and the pipeline skips read QC without it.
#
# Packages (pacman -S): mingw-w64-x86_64-toolchain mingw-w64-x86_64-zlib
#   mingw-w64-x86_64-bzip2 mingw-w64-x86_64-xz mingw-w64-x86_64-libdeflate
#   mingw-w64-x86_64-curl-winssl make autoconf automake
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

here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/tools/windows}
work=${WORK:-$here/tools/build-windows}
mkdir -p "$out" "$work"
jobs=$(nproc 2>/dev/null || echo 4)

fetch() {  # url
    local f=$work/$(basename "$1")
    [ -s "$f" ] || curl -fsSL --retry 5 --retry-delay 5 -o "$f" "$1"
    tar -xjf "$f" -C "$work"
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
  make -j"$jobs" samtools.exe )

echo "== bcftools $HTS_VERSION"
fetch "https://github.com/samtools/bcftools/releases/download/$HTS_VERSION/bcftools-$HTS_VERSION.tar.bz2"
( cd "$work/bcftools-$HTS_VERSION"
  ./configure --with-htslib="$work/htslib-$HTS_VERSION" --disable-bcftools-plugins
  make -j"$jobs" bcftools.exe )

echo "== minimap2 $MINIMAP2_VERSION"
fetch "https://github.com/lh3/minimap2/releases/download/v$MINIMAP2_VERSION/minimap2-$MINIMAP2_VERSION.tar.bz2"
( cd "$work/minimap2-$MINIMAP2_VERSION"
  make -j"$jobs" minimap2.exe )

echo "== collecting into $out"
rm -rf "$out"
mkdir -p "$out/licenses"
cp "$work/samtools-$HTS_VERSION/samtools.exe" "$work/bcftools-$HTS_VERSION/bcftools.exe" \
   "$work/minimap2-$MINIMAP2_VERSION/minimap2.exe" "$out/"
cp "$work/htslib-$HTS_VERSION/LICENSE" "$out/licenses/htslib-LICENSE"
cp "$work/samtools-$HTS_VERSION/LICENSE" "$out/licenses/samtools-LICENSE"
cp "$work/bcftools-$HTS_VERSION/LICENSE" "$out/licenses/bcftools-LICENSE"
cp "$work/minimap2-$MINIMAP2_VERSION/LICENSE.txt" "$out/licenses/minimap2-LICENSE"
cat > "$out/README.txt" <<TXT
Sequencing tools for the Genetic Health FASTQ pipeline and its 1000 Genomes
test sample, built for Windows x64 by packaging/windows-tools.sh:

  samtools $HTS_VERSION, bcftools $HTS_VERSION (htslib $HTS_VERSION)   https://www.htslib.org
  minimap2 $MINIMAP2_VERSION                                          https://github.com/lh3/minimap2

Each is MIT-licensed by its authors (see licenses/). The DLLs are the MinGW-w64
runtime libraries they load (zlib, bzip2, xz, libdeflate, libcurl, winpthreads, ...).
fastp is not included: it has no Windows build, and the pipeline skips read QC
without it.
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
