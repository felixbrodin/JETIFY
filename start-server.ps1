# Starts a tiny local HTTP server for the Jetify web app.
# Usage:  powershell -ExecutionPolicy Bypass -File start-server.ps1
# Then open:  http://localhost:8000
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$port = 8000
if ($args.Count -gt 0) { $port = [int]$args[0] }

$mime = @{
  ".html" = "text/html; charset=utf-8"
  ".js"   = "application/javascript; charset=utf-8"
  ".css"  = "text/css; charset=utf-8"
  ".json" = "application/json; charset=utf-8"
  ".csv"  = "text/csv; charset=utf-8"
  ".tsv"  = "text/tab-separated-values; charset=utf-8"
  ".png"  = "image/png"
  ".ico"  = "image/x-icon"
}

$listener = New-Object System.Net.HttpListener
$listener.Prefixes.Add("http://localhost:$port/")
$listener.Start()
Write-Host "Server at http://localhost:$port/  (Ctrl+C to stop)" -ForegroundColor Green

while ($listener.IsListening) {
  $ctx = $listener.GetContext()
  $req = $ctx.Request; $res = $ctx.Response
  try {
    $rel = [Uri]::UnescapeDataString($req.Url.AbsolutePath).TrimStart('/')
    if ($rel -eq "") { $rel = "index.html" }
    $path = [System.IO.Path]::GetFullPath((Join-Path $root $rel))
    if (-not $path.StartsWith([System.IO.Path]::GetFullPath($root), [System.StringComparison]::OrdinalIgnoreCase)) { throw "forbidden" }
    if (Test-Path -LiteralPath $path -PathType Leaf) {
      $ext = [System.IO.Path]::GetExtension($path).ToLower()
      $res.ContentType = if ($mime.ContainsKey($ext)) { $mime[$ext] } else { "application/octet-stream" }
      $bytes = [System.IO.File]::ReadAllBytes($path)
      $res.ContentLength64 = $bytes.Length
      $res.OutputStream.Write($bytes, 0, $bytes.Length)
      Write-Host ("200  " + $rel) -ForegroundColor Gray
    } else {
      $res.StatusCode = 404
      $res.OutputStream.Write([Text.Encoding]::UTF8.GetBytes("404"), 0, 3)
      Write-Host ("404  " + $rel) -ForegroundColor DarkYellow
    }
  } catch {
    $res.StatusCode = 500
    Write-Host ("500  " + $rel + "  " + $_.Exception.Message) -ForegroundColor Red
  }
  $res.Close()
}
