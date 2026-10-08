"""Read-only identity for local verification workers, including PID reuse."""
import os


def process_birth_token(pid):
    if isinstance(pid,bool) or not isinstance(pid,int) or pid<=0:return None
    if os.name!='nt':
        try:
            from pathlib import Path
            return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]
        except (OSError,IndexError):return None
    import ctypes
    from ctypes import wintypes
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.OpenProcess.argtypes=(wintypes.DWORD,wintypes.BOOL,wintypes.DWORD)
    kernel.OpenProcess.restype=wintypes.HANDLE
    kernel.GetProcessTimes.argtypes=(wintypes.HANDLE,*([ctypes.POINTER(wintypes.FILETIME)]*4))
    kernel.GetProcessTimes.restype=wintypes.BOOL
    kernel.GetExitCodeProcess.argtypes=(wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD))
    kernel.GetExitCodeProcess.restype=wintypes.BOOL
    kernel.CloseHandle.argtypes=(wintypes.HANDLE,)
    handle=kernel.OpenProcess(0x1000,False,pid)
    if not handle:return None
    try:
        code=wintypes.DWORD()
        if not kernel.GetExitCodeProcess(handle,ctypes.byref(code)) or code.value!=259:return None
        times=[wintypes.FILETIME() for _ in range(4)]
        if not kernel.GetProcessTimes(handle,*(ctypes.byref(item) for item in times)):return None
        return str((times[0].dwHighDateTime<<32)|times[0].dwLowDateTime)
    finally:kernel.CloseHandle(handle)


def is_verification_worker_running(job):
    if not job.get('verification_worker'):return False
    identity=job
    if job.get('worker_birth_token') is None:
        # Compatibility for an already running verifier registered before
        # identity metadata was introduced. The sidecar is local JSON only.
        from pathlib import Path
        import json,re
        from workspace_paths import RUNS_DIR
        identifier=job.get('id','')
        if not re.fullmatch('[0-9a-f]{16}',identifier):return False
        try:identity=json.loads((RUNS_DIR/identifier/'worker_identity.json').read_text(encoding='utf-8'))
        except (OSError,ValueError):return False
    return bool(identity.get('worker_birth_token') is not None and
                process_birth_token(identity.get('worker_pid'))==identity['worker_birth_token'])
