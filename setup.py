import os
from setuptools import setup

ext_modules = []
cmdclass = {}

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext as pybind11_build_ext

    class optional_build_ext(pybind11_build_ext):
        """Compiles C++ extension with fallback if compiler or build fails."""
        def build_extension(self, ext):
            try:
                super().build_extension(ext)
            except Exception as exc:
                print(f"\n[WARNING] C++ extension compilation failed for {ext.name}: {exc}")
                print("[WARNING] lr-rna-qc will operate using the pure-Python fallback engine.\n")

    ext_modules = [
        Pybind11Extension(
            "lr_rna_qc._fast_qc",
            [
                os.path.join("src", "cpp", "fast_qc.cpp"),
                os.path.join("src", "cpp", "bindings.cpp"),
            ],
            cxx_std=17,
        )
    ]
    cmdclass = {"build_ext": optional_build_ext}
except ImportError:
    pass

setup(
    ext_modules=ext_modules,
    cmdclass=cmdclass,
)
