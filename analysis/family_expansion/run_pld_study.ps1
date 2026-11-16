param(
  [Parameter(Mandatory)][string] $Rscript,
  [string] $SimSurveySource,
  [ValidateSet('audit','spatial','controlled','comparison')]
  [string] $Phase = 'comparison',
  [string] $Results
)
$ErrorActionPreference = 'Stop'
$rscriptPath = (Resolve-Path -LiteralPath $Rscript).Path
$rBin = Split-Path -Parent $rscriptPath
$rtermPath = @((Join-Path $rBin 'Rterm.exe'), (Join-Path $rBin 'x64/Rterm.exe')) |
  Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
if (!$rtermPath) { throw 'Cannot find Rterm.exe beside the supplied Rscript installation.' }
if ($Phase -in 'audit','spatial') {
  if (!$SimSurveySource) { throw 'Specify -SimSurveySource for the spatial simulation.' }
  $env:TINYAM_PLD_SIMSURVEY = (Resolve-Path -LiteralPath $SimSurveySource).Path
}
if (!$Results) {
  $folder = switch ($Phase) {
    'comparison' { 'pld-family-comparison' }
    'controlled' { 'pld-controlled' }
    default { 'pld-whole-area' }
  }
  $Results = Join-Path 'analysis/family_expansion/results' $folder
}
$directory = (New-Item -ItemType Directory -Path $Results -Force).FullName
$archive = Join-Path $directory 'supervisor'
New-Item -ItemType Directory -Force -Path $archive | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$resources = Join-Path $archive "resources-$stamp.csv"
$launch = @{
  FilePath = $rtermPath
  ArgumentList = @('--vanilla', '--slave', '--file=analysis/family_expansion/simulate_pld.R',
    '--args', $Phase, ('"{0}"' -f $directory))
  WindowStyle = 'Hidden'
  WorkingDirectory = (Get-Location).Path
  PassThru = $true
  RedirectStandardOutput = Join-Path $archive "run-$stamp.log"
  RedirectStandardError = Join-Path $archive "errors-$stamp.log"
}
$worker = Start-Process @launch
try { $worker.PriorityClass = 'BelowNormal' } catch { Write-Warning $_.Exception.Message }
$clock = [System.Diagnostics.Stopwatch]::StartNew()
$lowMemoryChecks = 0
$node = Get-Command node.exe -ErrorAction SilentlyContinue
$nodeMetrics = @'
const os = require("node:os");
function ticks() {
  let idle = 0;
  let total = 0;
  for (const cpu of os.cpus()) {
    idle += cpu.times.idle;
    total += Object.values(cpu.times).reduce((sum, ticks) => sum + ticks);
  }
  return {idle, total};
}
const start = ticks();
setTimeout(() => {
  const end = ticks();
  const used = 1 - (end.idle - start.idle) / (end.total - start.total);
  console.log(JSON.stringify({memory: os.freemem() / 2 ** 20, cpu: 100 * used}));
}, 500);
'@
try {
  while (!$worker.HasExited) {
    $worker.Refresh()
    if ($worker.HasExited) { break }
    $available = $null
    $processor = $null
    try {
      $memory = Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory -ErrorAction Stop
      $cpu = Get-CimInstance Win32_PerfFormattedData_PerfOS_Processor -ErrorAction Stop |
        Where-Object Name -EQ '_Total' | Select-Object -First 1
      $available = $memory.AvailableMBytes
      $processor = $cpu.PercentProcessorTime
    } catch {
      if ($node) {
        $system = & $node.Source -e $nodeMetrics | ConvertFrom-Json
        $available = $system.memory
        $processor = $system.cpu
      }
    }
    [pscustomobject]@{
      time = Get-Date -Format o
      elapsed_seconds = $clock.Elapsed.TotalSeconds
      system_available_mib = $available
      system_cpu_percent = $processor
      worker_cpu_seconds = $worker.TotalProcessorTime.TotalSeconds
      worker_working_set_mib = [math]::Round($worker.WorkingSet64 / 1MB, 1)
      worker_peak_mib = [math]::Round($worker.PeakWorkingSet64 / 1MB, 1)
    } | Export-Csv -LiteralPath $resources -Append -NoTypeInformation
    if ($null -ne $available -and $available -lt 2048) {
      $lowMemoryChecks++
    } else {
      $lowMemoryChecks = 0
    }
    if ($lowMemoryChecks -ge 2) {
      throw 'Study paused after two checks below 2 GiB available memory. Resume with the same command when memory is available.'
    }
    Start-Sleep -Seconds 15
  }
} finally {
  $worker.Refresh()
  if (!$worker.HasExited) { Stop-Process -Id $worker.Id }
}
$worker.WaitForExit()
exit $worker.ExitCode
