import os
import subprocess
import pathlib
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


extension = Extension(
    name="pyoorb",
    sources=["pyoorb.f90", "pyoorb.pyf"],
    include_dirs=["../build"]
)

class PyoorbBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        try:
            os.chdir("../")
            self.spawn(["make", "pyoorb", "-j4"])
        finally:
            os.chdir("./python")

        src = "../lib/" + self.get_ext_filename(ext.name)
        dst = self.get_ext_fullpath(ext.name)
        self.mkpath(os.path.dirname(dst))
        self.move_file(src, dst)


with open("../VERSION", "r") as f:
    version = f.read().strip()


setup(
    name="pyoorb",
    version=version,
    ext_modules=[extension],
    cmdclass={
        "build_ext": PyoorbBuild,
    }
)
