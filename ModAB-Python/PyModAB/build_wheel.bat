@echo off
REM Build wheel for PyPI distribution
pip install build
py -m build
echo.
echo Wheel created in dist/ directory
echo Upload to PyPI with: twine upload dist/*
