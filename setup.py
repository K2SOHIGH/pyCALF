from setuptools import setup, find_packages


setup(
    name='pyCALF',
    version='0.1',
    description='search calcyanin within a set of amino acid sequences',
    url='',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        "pyhmmer",
        "biopython",
        "numpy==1.22.4",
        "pyyaml==6.0",
        "pandas",
        # "gffutils",
    ],    
    include_package_data=True,
    entry_points = {
        'console_scripts': ['pycalf = pycalf.main:main'],
    },
    zip_safe=False
)
