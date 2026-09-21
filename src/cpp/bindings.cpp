#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "fast_qc.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_fast_qc, m) {
    m.doc() = "C++ accelerated FASTQ quality control engine for Oxford Nanopore RNA-seq";

    py::class_<lr_rna_qc::ReadQC>(m, "ReadQC")
        .def_readonly("mean_q", &lr_rna_qc::ReadQC::mean_q)
        .def_readonly("length", &lr_rna_qc::ReadQC::length)
        .def("__eq__", &lr_rna_qc::ReadQC::operator==)
        .def("__repr__", [](const lr_rna_qc::ReadQC& qc) {
            return "<ReadQC length=" + std::to_string(qc.length) +
                   " mean_q=" + std::to_string(qc.mean_q) + ">";
        });

    m.def("process_quality_string", &lr_rna_qc::process_quality_string,
          "Compute exact mean Phred quality score and read length using C++ LUT",
          py::arg("qual"));

    m.def("process_fastq_chunk", &lr_rna_qc::process_fastq_chunk,
          "Batch process FASTQ chunk vectors in C++",
          py::arg("headers"), py::arg("seqs"), py::arg("quals"));
}
