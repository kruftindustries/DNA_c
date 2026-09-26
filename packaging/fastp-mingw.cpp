// pwrite for MinGW-w64, for fastp (see fastp-mingw.h). Copied into fastp's
// src/ by windows-tools.sh so its Makefile's wildcard compiles it.
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <io.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>

extern "C" ssize_t fastp_mingw_pwrite(int fd, const void *buf, size_t count, long long offset)
{
    HANDLE h = reinterpret_cast<HANDLE>(_get_osfhandle(fd));
    if (h == INVALID_HANDLE_VALUE) {
        errno = EBADF;
        return -1;
    }
    OVERLAPPED ov;
    memset(&ov, 0, sizeof ov);
    ov.Offset = static_cast<DWORD>(offset & 0xffffffffu);
    ov.OffsetHigh = static_cast<DWORD>(static_cast<unsigned long long>(offset) >> 32);
    const DWORD chunk = count > 0x7fffffffu ? 0x7fffffffu : static_cast<DWORD>(count);
    DWORD written = 0;
    // A synchronous handle with an OVERLAPPED offset: the write lands at
    // `offset`, and concurrent calls with their own offsets do not race on
    // a shared file position. The caller loops on short writes.
    if (!WriteFile(h, buf, chunk, &written, &ov)) {
        errno = GetLastError() == ERROR_DISK_FULL ? ENOSPC : EIO;
        return -1;
    }
    return static_cast<ssize_t>(written);
}
#endif
