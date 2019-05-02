from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "asqcan",
    version = "0.3",
    author = "Daniel Bogema",
    author_email = "daniel.bogema@dpi.nsw.gov.au",
    description = ("An ASsembly, Quality Control and ANnotation pipeline"),
    license = "GPL-3.0",
    keywords = "genomics assembly annotation quality control",
    url = "https://github.com/bogemad/asqcan",
    py_modules=['asqcan'],
    scripts=['bin/asqcan'],
    long_description=read('README.md'),
)
