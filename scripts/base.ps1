$python_not_found = "Cannot find Python. Is Anaconda for Python 3.7 installed?"
$path1 = "Registry::HKEY_CURRENT_USER\Software\Python\ContinuumAnalytics\Anaconda37-64\InstallPath"
$path2 = "Registry::HKEY_LOCAL_MACHINE\Software\Python\ContinuumAnalytics\Anaconda37-64\InstallPath"
$value = "ExecutablePath"

# https://stackoverflow.com/questions/5648931
Function Test-RegistryValue {
    param(
            [Alias("PSPath")]
            [Parameter(Position = 0, Mandatory = $true, ValueFromPipeline = $true, ValueFromPipelineByPropertyName = $true)]
            [String]$Path
            ,
            [Parameter(Position = 1, Mandatory = $true)]
            [String]$Name
            ,
            [Switch]$PassThru
         ) 

        process {
            if (Test-Path $Path) {
                $Key = Get-Item -LiteralPath $Path
                    if ($Key.GetValue($Name, $null) -ne $null) {
                        if ($PassThru) {
                            Get-ItemProperty $Path $Name
                        } else {
                            $true
                        }
                    } else {
                        $false
                    }
            } else {
                $false
            }
        }
}
Function Activate-Anaconda {
	process {
		$p = $false
		if (Test-RegistryValue -Path $path1 -Name $value) {
				$p = (Test-RegistryValue -PassThru -Path $path1 -Name $value).$value
		} elseif (Test-RegistryValue -Path $path2 -Name $value) {
				$p = (Test-RegistryValue -PassThru -Path $path2 -Name $value).$value
		}

		if ($p) {
				$p = (Split-Path -Parent -Path $p)
				# https://github.com/BCSharp/PSCondaEnvs
				$env:Path = "$p;$p\Library\mingw-w64\bin;$p\Library\usr\bin;$p\Library\bin;$p\Scripts;$p\bin;" + $env:Path
				$env:CONDA_DEFAULT_ENV = "root"
				$env:CONDA_PREFIX = $p
		} else {
				Write-Error $python_not_found -ErrorAction Stop
		}
	}
}
