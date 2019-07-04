#!/usr/bin/env python
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##################################################################################################
# weighted NEW CDF
from scipy.stats import truncnorm
import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import scipy.integrate as integrate
from scipy import special
import scipy.optimize as opt
from sympy.solvers import solve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate g(r)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NM pdf/r**2


def Weightedpdf(fx, fmu, fsigma, fa, fb):
    # return normpdf(x,mu,sigma)/(x**2)
    return truncnorm.pdf(fx, fa, fb, loc=fmu, scale=fsigma)/(fx**2)

# integral NUMBER(0,2*MU)


def InteWeightedpdf(fmu, fsigma, xl, xr, summation, fa, fb):
    rr = np.linspace(xl, xr, summation)
    r = rr[1:]
    resultMatrix = Weightedpdf(r, fmu, fsigma, fa, fb)
    return np.sum(resultMatrix)*((xr-xl)/summation)


def awPDFDistribution(fx, fmu, fsigma, summation, fa, fb):
    fintegralPDF = InteWeightedpdf(fmu, fsigma, 0, 2*fmu, summation, fa, fb)
    # print(fintegralPDF)
    areaweightedPro = Weightedpdf(fx, fmu, fsigma, fa, fb)/fintegralPDF
    return areaweightedPro
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ending g(r)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate Phi(r)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# integrate from a to b


def InteAwpdf(fx, fmu, fsigma, xl, xr, summation, fa, fb):
    xx = np.linspace(xl, xr, summation)
    s = xx[1:]
    r = s.transpose()
    size = np.shape(r)
    resu = []
    dimension = r.ndim
    if dimension == 1:
        resultMatrix = awPDFDistribution(r, fmu, fsigma, summation, fa, fb)
        inn = np.true_divide((np.subtract(xr, xl)), summation)
        eachResult = np.sum(resultMatrix)*inn
        resu.append(eachResult)
    else:
        for ii in range(0, size[0]):
            resultMatrix = awPDFDistribution(
                r[ii, :], fmu, fsigma, summation, fa, fb)
            inn = np.true_divide((np.subtract(xr[ii], xl)), summation)
            eachResult = np.sum(resultMatrix)*inn
            resu.append(eachResult)
    return resu


def Phi(fx, fmu, fsigma, xr, summation, fa, fb):
    # fy= InteAwpdf(x,mu,sigma,0,2*mu,summation)  #1
    numerator = InteAwpdf(fx, fmu, fsigma, 0, xr, summation, fa, fb)
    #awCDF= fx/fy
    return numerator
# print(Phi(x,mu,sigma,x,interval,a,b))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end Phi(r)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AW Percent Point Function


def awGrainSize(fprobability, fx, mu, sigma, summation, a, b):
    def devi(fx):
        numerator = Phi(fx, mu, sigma, fx, summation, a, b)
        # fy= Phi(x,mu,sigma,2*mu,interval,a,b)
        #awCDF= fx/fy
        # deviValue = (awCDF-fprobability)**2
        deviValue = (np.subtract(numerator, fprobability))**2
        return deviValue
    fgs = opt.fmin(devi, mu-2*sigma)
    fgrainsize = fgs[0]
    if fgrainsize < mu-2*sigma:
        fgrainsize = mu-2*sigma
    elif fgrainsize > mu+2*sigma:
        fgrainsize = mu+2*sigma
    return fgrainsize


# run in this file
if __name__ == "__main__":
    # initial parameters
    mu = 20
    sigma = 8
    a1, a2 = mu-2*sigma, mu+2*sigma

    x1, x2, interval = a1, a2, 1000
    x = np.linspace(a1, a2, interval)
    # print(x)
    a, b = (a1-mu)/sigma, (a2-mu)/sigma

    print(awGrainSize(0.0001, x, mu, sigma, interval, a, b))
    # Plots
    # plt.figure(1)
    # plt.plot(x,truncnorm.pdf(x, a, b, loc=mu, scale=sigma),'b',label='normpdf')
    # plt.legend()

    # plt.figure(2)
    # plt.plot(x,Weightedpdf(x,mu,sigma,a,b),'r',label='pdf/x**2')
    # plt.legend()

    # plt.figure(3)
    # plt.plot(x,awPDFDistribution(x,mu,sigma,interval,a,b),'r',label='awPDF')
    # plt.legend()

    plt.figure(4)
    plt.plot(x, truncnorm.cdf(x, a, b, loc=mu,scale=sigma), 'b-', label='Normal Distribution')
    plt.plot(x, Phi(x, mu, sigma, x, interval, a, b), 'b--', label='Area Weighted Normal Distribution' )
    # ''' subplots
    f,(ax1,ax2) = plt.subplots(1,2)
    
    ax1.plot(x, truncnorm.cdf(x, a, b, loc=mu,scale=sigma), 'b-', label='ND of (20,8)')
    ax1.plot(x, Phi(x, mu, sigma, x, interval, a, b), 'b--', label='AW ND of (20,8)' )
    ax1.legend(loc='lower center',fontsize=14)
    ax1.set_xlabel('Radius of Grain Size',fontsize=14)
    ax1.set_ylabel('Probability',fontsize=14)
    ax1.set_title('CDF Curves of Two Distribution with (mu,sigma) being (20,8)',fontsize=14)
    
    
    ax2.plot(x, truncnorm.cdf(x, -2, 2, loc=20,scale=2), 'g-', label='ND of (20,2)')
    ax2.plot(x, Phi(x, 20, 2, x, interval, -2, 2), 'g--', label='AW ND of (20,2)' )
    ax2.set_xlim([16, 24])

    plt.legend(loc='lower center',fontsize=14)
    plt.xlabel('Radius of Grain Size',fontsize=14)
    plt.ylabel('Probability',fontsize=14)
    plt.title('CDF Curves of Two Distribution with (mu,sigma) being (20,2)',fontsize=14)
    plt.show()
