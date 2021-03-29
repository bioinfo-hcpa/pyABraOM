from setuptools import setup, find_packages


setup(
  name= "pyABraOM",
  version = '0.0.1',
  description= "'A Python API to communicate with ABraOM database",
  author ="Paola Carneiro and Felipe Colombeli",
  author_email= "pa0la_barcellosca@gmail.com and fcolombelli@inf.ufrgs.br ",
  url="abraompackage.com",
  license= "MIT License",
  install_requires=["requests>=2.23.0",
                    "pandas>=1.1.5",
                    "times"],
  py_modules=["Search_gene",
              "Search_Region",
              "Variant_ID", 
              'Searches'],
  packages =find_packages()
)
