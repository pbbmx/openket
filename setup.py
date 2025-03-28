from setuptools import setup, find_packages

setup(
    name="openket",
    version="0.1.0",
    author="IIMAS UNAM",
    description="Free software for manipulating quantum objects in Dirac notation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pbbmx/openket.git",
    packages=find_packages(include=["openket", "openket.*"]),
    install_requires=[
        "sympy>=1.0",
        "scipy>=1.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)