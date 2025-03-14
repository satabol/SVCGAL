from setuptools import setup, find_packages
from pathlib import Path, os

NAME = 'pySVCGAL'
VERSION = '0.0.25'
# - First Wrapper
DESCRIPTION = 'Python pySVCGAL.'

#LONG_DESCRIPTION = os.path.join( Path(__file__).parent, "README.md").read_text()
this_directory = Path(__file__).parent
LONG_DESCRIPTION = (this_directory / "README.md").read_text()

setup(
    name = NAME,
    version = VERSION,
    author = 'Satabol',
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    packages = find_packages(),
    include_package_data=True,
    package_data={f'{NAME}.clib':['*.so', '*.dll', '*.dylib']}, # https://stackoverflow.com/questions/70334648/how-to-correctly-install-data-files-with-setup-py
    install_requires = [],    
    keywords = ["straight skeleton 2D offset", "straight skeleton", "skeleton"],
    classifiers= [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",        
        "Programming Language :: Python :: 3",        
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ]
)