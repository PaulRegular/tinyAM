param(
  [Parameter(Mandatory)][string] $Rscript,
  [Parameter(Mandatory)][string] $SimSurveySource,
  [string] $Results = 'analysis/formula_expansion/results/pld-cod-zero-sim',
  [string] $ReportResults = '',
  [string] $FisheryCache = '',
  [int] $Replicates = 24,
  [int] $Grid = 40
)
$ErrorActionPreference = 'Stop'
$rscriptPath = (Resolve-Path -LiteralPath $Rscript).Path
$env:TINYAM_PLD_SIMSURVEY = (Resolve-Path -LiteralPath $SimSurveySource).Path
$env:TINYAM_PLD_REPORT_DIR = if ($ReportResults) {
  (New-Item -ItemType Directory -Path $ReportResults -Force).FullName
} else { '' }
$env:TINYAM_PLD_FISHERY_CACHE = if ($FisheryCache) {
  (Resolve-Path -LiteralPath $FisheryCache).Path
} else { '' }
$directory = (New-Item -ItemType Directory -Path $Results -Force).FullName
$archive = Join-Path $directory 'supervisor'
New-Item -ItemType Directory -Force -Path $archive | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$resources = Join-Path $archive "resources-$stamp.csv"
$launch = @{
  FilePath = $rscriptPath
  ArgumentList = @('--vanilla', 'analysis/formula_expansion/simulate_pld.R',
    'sampling_fit', $Replicates, $Grid, ('"{0}"' -f $directory))
  WindowStyle = 'Hidden'
  WorkingDirectory = (Get-Location).Path
  PassThru = $true
  RedirectStandardOutput = Join-Path $archive "run-$stamp.log"
  RedirectStandardError = Join-Path $archive "errors-$stamp.log"
}
$worker = Start-Process @launch
try { $worker.PriorityClass = 'BelowNormal' } catch { Write-Warning $_.Exception.Message }
$clock = [System.Diagnostics.Stopwatch]::StartNew()
try {
  while (!$worker.HasExited) {
    $worker.Refresh()
    if ($worker.HasExited) { break }
    $memory = Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory -ErrorAction SilentlyContinue
    $cpu = Get-CimInstance Win32_PerfFormattedData_PerfOS_Processor -ErrorAction SilentlyContinue |
      Where-Object Name -EQ '_Total' | Select-Object -First 1
    [pscustomobject]@{
      time = Get-Date -Format o
      elapsed_seconds = $clock.Elapsed.TotalSeconds
      system_available_mib = $memory.AvailableMBytes
      system_cpu_percent = $cpu.PercentProcessorTime
      worker_cpu_seconds = $worker.TotalProcessorTime.TotalSeconds
      worker_working_set_mib = [math]::Round($worker.WorkingSet64 / 1MB, 1)
      worker_peak_mib = [math]::Round($worker.PeakWorkingSet64 / 1MB, 1)
    } | Export-Csv -LiteralPath $resources -Append -NoTypeInformation
    Start-Sleep -Seconds 15
  }
} finally {
  $worker.Refresh()
  if (!$worker.HasExited) { Stop-Process -Id $worker.Id }
}
$worker.WaitForExit()
exit $worker.ExitCode
