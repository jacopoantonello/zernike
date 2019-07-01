Activate-Anaconda
if (-not $?) {
	exit
}
pip uninstall -y zernike
if (-not $?) {
	exit
}
rm dist\*.whl -ErrorAction SilentlyContinue
python setup.py bdist_wheel
if (-not $?) {
	exit
}
pip install (get-item .\dist\*.whl)
