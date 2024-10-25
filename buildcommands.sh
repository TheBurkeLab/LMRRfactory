################################## USER INPUTS #########################################
VERSION=0.0.39
COMMITMSG="update to version 0.0.39"
APIKEY=pypi-AgENdGVzdC5weXBpLm9yZwIkODkxMDdkYTgtMzlkOC00MjU4LWI2OGUtOWVjMmU2NjE1OTFiAAIqWzMsImViNWY0ZmQwLTkzNjItNDExZS04ZWViLTI0ZThkNDk1ODcyNyJdAAAGIEFgJrbxSO9oefkf0Z4A3rnmt_ScMPKs4-UQwaV2REkv
########################################################################################

export TWINE_USERNAME=__token__
export TWINE_PASSWORD=$APIKEY

sed -i 's/version = .*/version = '\'$VERSION\''/' pyproject.toml
rm -rf dist
rm -rf src/LMRRfactory.egg-info

python3 -m build
python setup.py sdist bdist_wheel

git add .
git commit -m "$COMMITMSG"
git push origin main

python3 -m twine upload --repository testpypi dist/*
# python3 -m twine upload --repository pypi dist/* # uncomment to upload to real PYPI

pip uninstall -y LMRRfactory
pip install -i https://test.pypi.org/simple/ LMRRfactory==$VERSION
pip install -i https://test.pypi.org/simple/ LMRRfactory==$VERSION


# # chmod +x buildcommands.sh python setup.py sdist bdist_wheel
# # # scons test-kinetics toolchain=msvc verbose_tests=y -j4 googletest=submodule > testlog.txt 2>&1

