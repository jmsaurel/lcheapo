rm -r dist
rm -r build
rm -r *.egg-info
mkdir dist
rm lcheapo_obspy/_examples/.DS_Store
rm lcheapo_obspy/data/.DS_Store
python setup.py sdist
python setup.py bdist_wheel
twine check dist/*
#twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose
twine upload dist/*
