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
Function Run-Python {
    param(
            [Parameter(Position = 0, Mandatory = $true, ValueFromPipeline = $true, ValueFromPipelineByPropertyName = $true)]
            [String]$Cmd
         ) 
	process {
		$path1 = "Registry::HKEY_CURRENT_USER\Software\Python\ContinuumAnalytics\Anaconda36-64\InstallPath"
		$value1 = "ExecutablePath"
		# $path2 = "Registry::HKEY_LOCAL_MACHINE\Software\Python\ContinuumAnalytics\Anaconda36-64\InstallPath"
		# $value2 = "ExecutablePath"
		
		if (Test-RegistryValue -Path $path1 -Name $value1) {
				$p = (Test-RegistryValue -PassThru -Path $path1 -Name $value1).$value1
				iex "& $p $cmd"
		} else {
				Write-Host -ForegroundColor Red "Cannot find Python. Is Anaconda for Python 3 installed?"
		}
    }
}
Function Activate-Anaconda {
	process {
		$path1 = "Registry::HKEY_CURRENT_USER\Software\Python\ContinuumAnalytics\Anaconda36-64\InstallPath"
		$value1 = "ExecutablePath"
		# $path2 = "Registry::HKEY_LOCAL_MACHINE\Software\Python\ContinuumAnalytics\Anaconda36-64\InstallPath"
		# $value2 = "ExecutablePath"
		
		if (Test-RegistryValue -Path $path1 -Name $value1) {
				$p = (Test-RegistryValue -PassThru -Path $path1 -Name $value1).$value1
                $p = (Split-Path -Parent -Path $p)
                $env:Path += "$p;$p\Scripts"
		} else {
				Write-Host -ForegroundColor Red "Cannot find Python. Is Anaconda for Python 3 installed?"
		}
    }
}
