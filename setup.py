from setuptools import setup, find_packages

VERSION = '0.0.4'
DESCRIPTION = 'Auto-generate LMR-R reactions and mechanisms'
LONG_DESCRIPTION = "Auto-generates LMR-R reactions (and mechanisms) according to the user's choice of Plog, Troe, or Chebyshev sub-formats"

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="LMRRfactory", 
        version=VERSION,
        author="Patrick Singal",
        author_email="p.singal@columbia.edu",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        include_package_data=True,
        package_data={"LMRRfactory":['data','ext/cantera/build/python']},
        # install_requires=['cantera','numpy','pyyaml','scipy'],
        install_requires=['numpy','pyyaml','scipy'],

        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)