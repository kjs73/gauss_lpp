from __future__ import division
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.special import hyp2f1, gamma
import subprocess
import os
import sys

def to_string(inp, digits_after_point = 16):
    format_string = "{0:."
    format_string += str(digits_after_point)
    format_string += "f}"
    return format_string.format(inp)

def save_pdf(plt, file_name):
    pdf = PdfPages(file_name)
    plt.savefig(pdf, format="pdf")
    pdf.close()
    plt.close()
    
class compare_and_plot(object):
    def __init__(self, file_base = "gauss_lpp", nr_points = 200):
        self.file_base = file_base
        self.nr_points = nr_points
        self.phi_min = (1 / self.nr_points)
        self.phi_max = (1 - self.phi_min)
        self.phi_min *= np.pi
        self.phi_max *= np.pi
        self.set_up_evaluations()
        self.get_scipy_data()
        self.get_cpp_data()
        self.generate_plots()
    def set_up_evaluations(self):
        self.kappa = range(1, 9)
        self.scipy_data = []
        self.cpp_data = []
        self.phi = np.linspace(self.phi_min, self.phi_max, self.nr_points)
    def get_scipy_data(self):
        for kappa in self.kappa:
            print "gsl, kappa:", kappa
            self.scipy_data.append([self.get_lpp(kappa, phi) for phi in self.phi])
    def get_cpp_data(self):
        for kappa in self.kappa:
            print "integral, kappa:", kappa
            self.cpp_data.append([self.get_cpp_lpp(kappa, phi) for phi in self.phi])
    def generate_plots(self):
        plt.rcParams.update({'font.size': 17})
        plt.rcParams.update({'figure.autolayout': True})
        plt.xlabel(r"$\phi / \pi$")
        plt.ylabel(r"Probability $P_{\kappa}(\phi)$")
        for kidx in xrange(len(self.kappa)):
            plt.plot(self.phi / np.pi, self.scipy_data[kidx])
        save_pdf(plt, self.file_base + "_scipy_plot.pdf")
        plt.xlabel(r"$\phi / \pi$")
        plt.ylabel(r"Probability $P_{\kappa}(\phi)$")
        for kidx in xrange(len(self.kappa)):
            plt.plot(self.phi / np.pi, self.cpp_data[kidx])
        save_pdf(plt, self.file_base + "_cpp_plot.pdf")
        for kidx in xrange(len(self.kappa)):
            plt.xlabel(r"$\phi / \pi$")
            plt.ylabel(r"Absolute error $| P_{\kappa, \mathrm{C++}}(\phi) - P_{\kappa, \mathrm{Scipy}}(\phi) |$")
            plt.plot(self.phi / np.pi, [np.abs(ci - si) for (ci, si) in zip(self.cpp_data[kidx], self.scipy_data[kidx])], "k")
            save_pdf(plt, self.file_base + "_error_kappa_" + to_string(self.kappa[kidx], 0) + ".pdf")
        self.test_cpp_at_phi_0p5_pi()
        self.test_cpp_at_kappa_2()
    def test_cpp_at_phi_0p5_pi(self):
        plt.xlabel(r"$\kappa$")
        plt.ylabel(r"Absolute error $| P_{\kappa, \mathrm{C++}}(\pi/2) - 1/2 |$")
        kappa = np.linspace(1, 8, 200)
        plt.plot(kappa, [np.abs(self.get_cpp_lpp(kappai, np.pi / 2) - 1 / 2) for kappai in kappa], "k")
        save_pdf(plt, self.file_base + "_error_phi_0p5_pi.pdf")
    def test_cpp_at_kappa_2(self):
        def analytic_result(phi):
            return 1 + np.cos(phi) * np.sin(phi) / np.pi - phi / np.pi
        plt.xlabel(r"$\phi / \pi$")
        plt.ylabel(r"Absolute error $| P_{\kappa, \mathrm{C++}}(\phi) - (1 + \cos(\phi)\sin(\phi)/\pi - \phi/\pi)|$")
        plt.plot(self.phi / np.pi, [np.abs(self.get_cpp_lpp(2, phi) - analytic_result(phi)) for phi in self.phi], "k")
        save_pdf(plt, self.file_base + "_error_vs_analytical_kappa_2.pdf")
    def get_lpp(self, kappa, phi):
        def cot(x):
            return 1 / np.tan(x)
        return 1 / 2 + gamma(4 / kappa) / np.sqrt(np.pi) / gamma((8 - kappa) / (2 * kappa)) * cot(phi) * hyp2f1(1 / 2, 4 / kappa, 3 / 2, - cot(phi) ** 2)
    def get_cpp_lpp(self, kappa, phi):
        cwd = os.path.abspath(".")
        out = subprocess.Popen(["./get_lpp", to_string(kappa), to_string(phi)], stdout = subprocess.PIPE, cwd = cwd).communicate()[0]
        return float(out)
    
if __name__ == "__main__":
    compare_and_plot()