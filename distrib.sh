rm -r dist
rm -r build
rm -r *.egg-info
mkdir dist
rm lcheapo/_examples/.DS_Store
rm lcheapo/data/.DS_Store
python setup.py sdist
python setup.py bdist_wheel
twine check dist/*
#twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose
twine upload dist/*
