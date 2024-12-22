###############################################################
#
# A python script to generate hydrostatic white dwarfs and put
# them into readable table to perform stellar simulations
#
###############################################################

#import packages#
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

###############################################################

#adiabatic constants#
ye = 0.5
gamma = 4/3

#constants#
gconst = 6.67430e-8
clight = 2.99792458e10
solar = 1.98847e33
rsolar = 6.96342e10

# Conversion between units #
lencgs2code = (clight**2)/(solar*gconst)
masscgs2code = (1.0/solar)
tcgs2code = (clight**3)/(solar*gconst)

# Derived conversion #
rhocgs2code = (masscgs2code/lencgs2code**3)
# taucgs2code = (masscgs2code*lencgs2code**2/tcgs2code**2)

energycgs2code = (1/clight**2)
taucgs2code = (rhocgs2code*energycgs2code)

# convert polytropic constants#
me2 = 9.1093837015E-28*masscgs2code
mb2 = 1.66053906660E-24*masscgs2code
h_bar = (1.054571817E-27)*(lencgs2code**2*masscgs2code/tcgs2code)
# kpoly = (3**(2/3)*math.pi**(4/3)*h_bar**2*ye**(5/3))/(5*me2*mb2**(5/3))
kpoly = 0.75*(3/8/math.pi)**(1/3)*h_bar*2*np.pi*(ye/mb2)**(4/3)
print(kpoly)
###############################################################
# other paramteres

#factor of atmpsheric density#
rhofac = 1.0e-6

# number of grid#
n_grid = 1200

# start and ending coordiante #
r_start = 1e-50
r_end = 1200

# central density #
rhoc = 5.0e10*rhocgs2code
rhoa = rhoc*rhofac

###############################################################
# set array #

#coordinate#
rface = np.zeros(n_grid+7)
rgrid = np.zeros(n_grid+6)
rho = np.zeros(n_grid+6)
m = np.zeros(n_grid+6)
p = np.zeros(n_grid+6)

###############################################################
# initialize #

#grid size#
dr = (r_end - r_start)/(n_grid)

# face coordinate #
for i in range (0, len(rface)):
  rface[i] =  r_start + (i-3)*dr

# cell center coordinate #
for i in range (0, len(rgrid)):
  rgrid[i] =  0.5*(rface[i+1] + rface[i])

###############################################################
# define function #

def get_drhodr (m_in, r_in, rho_in):
  drhodr = - m_in*rho_in**(2 - gamma)/(kpoly*gamma*r_in**2)
  return drhodr

def get_dmdr (r_in, rho_in):
  dmdr = 4*math.pi*r_in**2*rho_in
  return dmdr

###############################################################
# loop over #

m[3] = (4/3)*math.pi*rhoc*rface[4]**3
rho[3] = rhoc

# do the rk-4 integration#
for n in range (3, n_grid+3):
  
  # get initial input #
  r_0 = rgrid[n]
  m_0 = m[n]
  rho_0 = rho[n]

  #atmosphere#
  if(rho_0 < rhoa): rho_0 = rhoa

  # get k1 #
  k1_rho = get_drhodr(m_0, r_0, rho_0)
  k1_m = get_dmdr(r_0, rho_0)

  #atmosphere#
  if(rho_0 < rhoa): 
    k1_rho = 0
    k1_m = 0

  # update #
  r_1 = rgrid[n] + dr/2
  m_1 = m[n] + dr*k1_m/2
  rho_1 = rho[n] + dr*k1_rho/2

  #atmosphere#
  if(rho_1 < rhoa): rho_1 = rhoa

  # get k2 #
  k2_rho = get_drhodr(m_1, r_1, rho_1)
  k2_m = get_dmdr(r_1, rho_1)

  #atmosphere#
  if(rho_1 < rhoa): 
    k2_rho = 0
    k2_m = 0

  # update #
  r_2 = rgrid[n] + dr/2
  m_2 = m[n] + dr*k2_m/2
  rho_2 = rho[n] + dr*k2_rho/2

  #atmosphere#
  if(rho_2 < rhoa): rho_2 = rhoa

  # get k2 #
  k3_rho = get_drhodr(m_2, r_2, rho_2)
  k3_m = get_dmdr(r_2, rho_2)

  #atmosphere#
  if(rho_2 < rhoa): 
    k3_rho = 0
    k3_m = 0

  # update #
  r_3 = rgrid[n] + dr
  m_3 = m[n] + dr*k3_m
  rho_3 = rho[n] + dr*k3_rho

  #atmosphere#
  if(rho_3 < rhoa): rho_3 = rhoa

  # get k2 #
  k4_rho = get_drhodr(m_3, r_3, rho_3)
  k4_m = get_dmdr(r_3, rho_3)

  #atmosphere#
  if(rho_3 < rhoa): 
    k4_rho = 0
    k4_m = 0

  # next step #
  m[n+1] = m[n] + dr*(k1_m + 2*k2_m + 2*k3_m + k4_m)/6
  rho[n+1] = rho[n] + dr*(k1_rho + 2*k2_rho + 2*k3_rho + k4_rho)/6

  #atmosphere#
  if(rho[n+1] < rhoa): rho[n+1] = rhoa

###############################################################
# boundary condition #
rho[0] = rho[5]
rho[1] = rho[4]
rho[2] = rho[3]
rho[n_grid+3] = rho[n_grid+2]
rho[n_grid+4] = rho[n_grid+2]
rho[n_grid+5] = rho[n_grid+2]

###############################################################
# find pressure #
p[:] = kpoly*rho[:]**(gamma)

###############################################################
#output#

#convert face coordinate#
rface = pd.DataFrame(rface)
rface.to_csv('grid.dat', index=False, header=False, sep=' ')

# Dataframe #
hydro = pd.DataFrame((rho, p)).T
hydro.to_csv('hydro.dat', index=False, header=False, sep=' ')

###############################################################
