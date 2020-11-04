
from setuptools import find_packages, setup
import versioneer

setup(
    name='dlp_qc_utils',
    description='dlp qc',
    packages=find_packages(),
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
