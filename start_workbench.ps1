# Use the package virtual environment when it exists. Pass all arguments on.
# Example: .\start_workbench.ps1 --port 8769 --workers 2
$ErrorActionPreference = 'Stop'
$taskPython = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $taskPython -PathType Leaf)) {
    $taskPython = 'python'
}
& $taskPython (Join-Path $PSScriptRoot 'start_workbench.py') @args
exit $LASTEXITCODE
