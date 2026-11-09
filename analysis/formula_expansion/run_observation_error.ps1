param(
  [string]$SimSurveySource,
  [string]$Rscript = "Rscript",
  [ValidateSet("Report", "SamplingReport", "PopulationPilot", "PopulationDiagnostic", "PopulationExtension", "Pilot", "Study", "ZPilot", "ZStudy", "ZReport")] [string]$Phase = "Report",
  [ValidateSet(0, 20, 40)] [int]$Grid = 0,
  [ValidateSet(3, 24)] [int]$PopulationLimit = 24
)

$ErrorActionPreference = "Stop"
$root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "../..")).Path
$source = ""
if ($Phase -notin @("Report", "SamplingReport", "ZReport")) {
  if (!$SimSurveySource) { throw "Set -SimSurveySource to the pinned SimSurvey source folder." }
  $source = (Resolve-Path -LiteralPath $SimSurveySource).Path
}
$rExe = (Get-Command $Rscript -CommandType Application).Source
# Monitor the real R process, rather than the Windows architecture launcher.
$x64 = Join-Path (Split-Path -Parent $rExe) "x64/Rscript.exe"
if (Test-Path -LiteralPath $x64) { $rExe = $x64 }
$populationPhase = $Phase -in @("PopulationPilot", "PopulationDiagnostic", "PopulationExtension", "Report")
$zPhase = $Phase -in @("ZPilot", "ZStudy", "ZReport")
$results = if ($zPhase) {
  Join-Path $PSScriptRoot "results/observation_error_z_covariance"
} elseif ($populationPhase) {
  Join-Path $PSScriptRoot "results/observation_error_population_pilot"
} else {
  Join-Path $PSScriptRoot "results/observation_error_seed21"
}
New-Item -ItemType Directory -Force -Path $results | Out-Null
$run = Get-Date -Format "yyyyMMdd-HHmmss-fff"
$gate = Join-Path $results "gate-$run.txt"
$resourceLog = Join-Path $results "resources-$run.csv"
$lock = [IO.File]::Open((Join-Path $results "launcher.lock"), "OpenOrCreate", "ReadWrite", "None")
$worker = $null

function Get-StudyLoad {
  $memory = Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory
  $cpu = Get-CimInstance Win32_PerfFormattedData_PerfOS_Processor -Filter "Name = '_Total'"
  [pscustomobject]@{
    time = [DateTime]::UtcNow.ToString("o")
    ram_gib = [math]::Round($memory.AvailableBytes / 1GB, 3)
    commit_gib = [math]::Round(($memory.CommitLimit - $memory.CommittedBytes) / 1GB, 3)
    cpu_pct = [double]$cpu.PercentProcessorTime
  }
}

function Get-Headroom($load) {
  [math]::Min($load.ram_gib, $load.commit_gib)
}

function Write-Resource($load, $phaseName, $gridSize, $processId, $peak, $state) {
  [pscustomobject]@{
    time = $load.time
    phase = $phaseName
    grid = $gridSize
    worker_pid = $processId
    ram_gib = $load.ram_gib
    commit_gib = $load.commit_gib
    cpu_pct = $load.cpu_pct
    worker_peak_gib = [math]::Round($peak, 3)
    state = $state
  } | Export-Csv -LiteralPath $resourceLog -Append -NoTypeInformation
}

function Write-GateSignal($state) {
  # The R worker briefly opens this file to read its current signal.
  for ($attempt = 1; $attempt -le 10; $attempt++) {
    try {
      Set-Content -LiteralPath $gate -Value $state -Encoding ascii
      return
    } catch [System.IO.IOException] {
      if ($attempt -eq 10) { throw }
      Start-Sleep -Milliseconds 100
    }
  }
}
function Wait-Headroom($needed, $phaseName, $gridSize) {
  $ready = 0
  $waiting = $false
  while ($ready -lt 2) {
    $load = Get-StudyLoad
    if ((Get-Headroom $load) -ge $needed -and $load.cpu_pct -le 80) {
      $ready++
    } else {
      $ready = 0
      if (!$waiting) {
        Write-Host ("Waiting for {0:F2} GiB of RAM and commit headroom for {1}." -f $needed, $phaseName)
        $waiting = $true
      }
    }
    Write-Resource $load $phaseName $gridSize $null 0 "preflight"
    if ($ready -lt 2) { Start-Sleep -Seconds 10 }
  }
  $load
}

function Invoke-StudyWorker($phaseName, $gridSize, $needed, $knownPeak = 0) {
  $attempt = 0
  $peak = [double]$knownPeak
  while ($true) {
    Wait-Headroom $needed $phaseName $gridSize | Out-Null
    $attempt++
    $logName = "$phaseName-grid$gridSize-$run-$attempt"
    $stdout = Join-Path $results "$logName.log"
    $stderr = Join-Path $results "$logName.err.log"
    $envForStudy = @{
      TINYAM_OE_SIMSURVEY = $source
      TINYAM_OE_PHASE = $phaseName
      TINYAM_OE_GRID = [string]$gridSize
      TINYAM_OE_GATE = $gate
      TINYAM_OE_POPULATION_LIMIT = [string]$PopulationLimit
      OMP_NUM_THREADS = "1"
      OPENBLAS_NUM_THREADS = "1"
      MKL_NUM_THREADS = "1"
      BLIS_NUM_THREADS = "1"
      VECLIB_MAXIMUM_THREADS = "1"
    }
    $scriptPath = "analysis/formula_expansion/simulate_observation_error.R"
    $launch = @{
      FilePath = $rExe
      ArgumentList = @("--vanilla", $scriptPath)
      WorkingDirectory = $root
      Environment = $envForStudy
      WindowStyle = "Hidden"
      PassThru = $true
      RedirectStandardOutput = $stdout
      RedirectStandardError = $stderr
    }
    Write-GateSignal "go"
    $script:worker = Start-Process @launch
    $null = $worker.Handle
    if (!$worker.HasExited) { $worker.PriorityClass = "BelowNormal" }
    Write-Host ("Started {0} on the {1}x{1} grid. Log is {2}." -f $phaseName, $gridSize, $stdout)

    $belowCritical = 0
    $cpuSince = $null
    $stoppedForMemory = $false
    while (!$worker.HasExited) {
      $load = Get-StudyLoad
      $worker.Refresh()
      if ($worker.HasExited) { break }
      $workerMemory = [math]::Max($worker.PeakWorkingSet64, $worker.PrivateMemorySize64) / 1GB
      $peak = [math]::Max($peak, $workerMemory)
      $headroom = Get-Headroom $load
      if ($load.cpu_pct -gt 80) {
        if ($null -eq $cpuSince) { $cpuSince = [DateTime]::UtcNow }
      } else {
        $cpuSince = $null
      }
      $cpuBusy = $null -ne $cpuSince -and ([DateTime]::UtcNow - $cpuSince).TotalSeconds -ge 60
      $state = if ($headroom -lt 3 -or $cpuBusy) { "wait" } else { "go" }
      if ($headroom -lt 2) { $belowCritical++ } else { $belowCritical = 0 }
      if ($belowCritical -ge 2) {
        $state = "stop"
        if (!$worker.HasExited) { Stop-Process -Id $worker.Id }
        $stoppedForMemory = $true
      }
      Write-GateSignal $state
      Write-Resource $load $phaseName $gridSize $worker.Id $peak $state
      if ($stoppedForMemory) { break }
      $null = $worker.WaitForExit(10000)
    }
    $worker.WaitForExit()
    $status = $worker.ExitCode
    $worker.Dispose()
    $script:worker = $null
    if ($stoppedForMemory) {
      if ($phaseName -eq "pilot") { return -$peak }
      $needed = [math]::Max($needed, $peak + 3)
      Write-Host "Paused $phaseName for memory pressure. Completed checkpoints will be reused."
      continue
    }
    if ($status -ne 0) {
      throw "$phaseName exited with code $status. See $stdout and $stderr."
    }
    return $peak
  }
}

try {
  if ($zPhase) {
    if ($Grid -eq 20) { throw "The Z covariance study uses the 40x40 grid." }
    Invoke-StudyWorker ($Phase.ToLowerInvariant()) 40 5 | Out-Null
    return
  }
  if ($populationPhase) {
    if ($Grid -eq 20) { throw "The population pilot uses the 40x40 grid." }
    Invoke-StudyWorker ($Phase.ToLowerInvariant()) 40 5 | Out-Null
    return
  }
  if ($Phase -eq "SamplingReport") {
    $reportGrid = if ($Grid -eq 0) { 40 } else { $Grid }
    Invoke-StudyWorker "samplingreport" $reportGrid 5 | Out-Null
    return
  }
  $requestedGrid = $Grid
  if ($Grid -eq 0) {
    $Grid = if ($Phase -eq "Pilot" -or
      (Test-Path -LiteralPath (Join-Path $results "pilot-grid40.json"))) { 40 } else { 20 }
  }
  $pilotProfile = Join-Path $results "pilot-grid$Grid.json"
  if ($Phase -eq "Pilot") {
    while (!(Test-Path -LiteralPath $pilotProfile)) {
      $peak = Invoke-StudyWorker "pilot" $Grid 5
      if ($peak -lt 0) {
        if ($Grid -ne 40 -or $requestedGrid -ne 0) { throw "Pilot stopped under memory pressure." }
        Write-Host "The 40x40 pilot reached the critical memory limit. Trying the 20x20 design."
        $Grid = 20
        $pilotProfile = Join-Path $results "pilot-grid20.json"
        continue
      }
      @{ peak_gib = $peak; grid = $Grid; resources = $resourceLog } |
        ConvertTo-Json | Set-Content -LiteralPath $pilotProfile -Encoding utf8
    }
  }
  if (!(Test-Path -LiteralPath $pilotProfile)) { throw "Complete the pilot for this grid first." }
  $pilotPeak = (Get-Content -LiteralPath $pilotProfile -Raw | ConvertFrom-Json).peak_gib
  if ($Phase -eq "Study") {
    Invoke-StudyWorker "study" $Grid ([math]::Max(5.0, [double]($pilotPeak + 3))) | Out-Null
  }
} finally {
  if ($null -ne $worker -and !$worker.HasExited) { Stop-Process -Id $worker.Id }
  $lock.Dispose()
  if (Test-Path -LiteralPath $gate) { Remove-Item -LiteralPath $gate }
}
