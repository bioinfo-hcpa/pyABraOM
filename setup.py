from setuptools import setup, find_packages


setup(
  name= "pyabraom",
  version="1.0",
  description= "A package to search variants in ABraOM database",
  author ="Paola Carneiro and Felipe Colombeli",
  author_email= "pa0la_barcellosca@gmail.com and fcolombelli@inf.ufrgs.br ",
  url="abraompackage.com",
  license= "MIT License",
  install_requires=["requests","pandas","times"],
  py_modules=["Search_gene","Search_Region","Variant_ID"],
  packages =find_packages()
)
