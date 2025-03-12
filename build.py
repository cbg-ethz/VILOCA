from glob import glob
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess

def build(setup_kwargs):
    # Build the Rust component using maturin
    subprocess.run(["maturin", "develop", "--release", "--bindings", "pyo3"], check=True)

    # Define the C++ extension module using pybind11
    ext_modules = [
        Pybind11Extension(
            "libshorah",
            sources=sorted(glob("src/cpp_module/*.cpp")),
            include_dirs=["src/cpp_module"],
            libraries=["hts"],  # Link external libraries if needed
            undef_macros=["HAVE_POPCNT"],
            extra_compile_args=["-std=c++11"]
        ),
    ]

    # Update setup kwargs for setuptools
    setup_kwargs.update({
        "ext_modules": ext_modules,
        "cmdclass": {"build_ext": build_ext},
        "zip_safe": False,
    })
