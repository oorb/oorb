import os
import subprocess
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
            self.spawn(["make", "-j4"])            
            self.spawn(["make", "pyoorb", "-j4"])
        finally:
            os.chdir("./python")

        src = "../lib/" + self.get_ext_filename(ext.name)
        dst = self.get_ext_fullpath(ext.name)
        self.mkpath(os.path.dirname(dst))
        self.move_file(src, dst)


def deduce_version():
    # This is a gnarly hack, but it ensures consistency.
    stdout = subprocess.PIPE
    cmd_output = subprocess.run(
        ["./build-tools/compute-version.sh",  "-u"],
        stdout=stdout,
    )
    cmd_output.check_returncode()
    return cmd_output.stdout.decode("utf8").strip()


setup(
    name="pyoorb",
    version=deduce_version(),
    ext_modules=[extension],
    install_requires=["numpy"],
    license="GPL3",
    cmdclass={
        "build_ext": PyoorbBuild,
    }
)
