# https://stackoverflow.com/questions/5648931
Function Get-Anaconda-Path {
    process {
        $bases = "Registry::HKEY_LOCAL_MACHINE\Software\Python\ContinuumAnalytics", "Registry::HKEY_CURRENT_USER\Software\Python\ContinuumAnalytics"
        $found = -1
        for ($i = 0; $i -lt $bases.Count; $i++) {
            if (Test-Path $bases[$i]) {
                $found = $i
                break;
            }
        }

        if ($found -lt 0) {
            Write-Error "Cannot find Anaconda. Is Anaconda for Python 3 installed (https://www.anaconda.com/download)?" -ErrorAction Stop
        }

        $base = $bases[$i]
        $subs = (Get-Item -LiteralPath $base).GetSubKeyNames()
        $addr = ""
        for ($i = 0; $i -lt $subs.Count; $i++) {
            $item = Get-Item -LiteralPath ($base + "\" + $subs[$i])
            $version = $item.GetValue("SysVersion") -as [decimal]
            if ($version -gt 3) {
                $addr = $base + "\" + $subs[$i]
            }
        }

        if ($addr.Length -lt 1) {
            Write-Error "Anaconda for Python 3+ must be installed (Anaconda for Python 2 is not supported)" -ErrorAction Stop
        }
        else {
            $item = Get-Item -LiteralPath "$addr\InstallPath"
            $item = ($item).GetValue("ExecutablePath")
            return $item
        }
    }
}
Function Activate-Anaconda {
	process {
		$p = Get-Anaconda-Path
        $p = (Split-Path -Parent -Path $p)
		# https://github.com/BCSharp/PSCondaEnvs
		$env:Path = "$p;$p\Library\mingw-w64\bin;$p\Library\usr\bin;$p\Library\bin;$p\Scripts;$p\bin;" + $env:Path
		$env:CONDA_DEFAULT_ENV = "root"
        $env:CONDA_PREFIX = $p
	}
}
