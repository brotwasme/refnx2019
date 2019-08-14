import numpy as np
import matplotlib.pyplot as plt
import matplotlib
 
from scipy import special
from refnx.dataset import ReflectDataset
from refnx.analysis import CurveFitter, Objective, Transform, parameter, possibly_create_parameter, Parameters
from refnx.reflect import ReflectModel, SLD, Erf, Component
from refnx.refnx import Spline, structure

# |       .        o           .       |
# air-pro, phi-pho, center, mirror (i-o), pro-wat 
#       a-w
#   (x/a)^2 + (y/b)^2 = 1 <-ellipse
# air-pro to cetre = a, b is from area taken up

class param():
    def __init__(real=0, image=0, name='unknown', thickness=0,) 
        self.real = real
        
        

class eggs(Component):#
    def __init__(self, air, hydropho, hydrophil, water, tilt):
        super(eggs)
        names = ['air', 'hydropho', 'hydrophil', 'water']
        more = ['tilt']
        pass

class bsla_thesis(Component): #as in thesis
    def __init__(self, i_as=11.1, sig_as=0.1, i_ap=24.1, th_p= 25.1, i_ps=24.1, n_w=500.1, d_sld=2.01, a=13, b=16.5, c=27.5):
        #interface air solvent, interface width air solvent, interface air protien, interface protien solvent, number of water molecules, sld of deuterated protein 
        super(bsla_thesis).__init__()
        name='bsla'
        self.name = name
        self.interface_air_protein = possibly_create_parameter(0.0,
                name='%s - interface_air_protein' % name)
        self.interface_protein_solvent = possibly_create_parameter(11.1,
                name='%s - interface_protein_solvent' % name)
        self.protein_length = possibly_create_parameter(25.1,
                name='%s - protein_length' % name)
        self.number_of_water_molecules = possibly_create_parameter(500.1,
                name='%s - number_of_water_molecules' % name)
        self.interface_width_air_solvent = possibly_create_parameter(0.1,
                name='%s - interface_width_air_solvent' % name)
        self.interface_width_protein_solvent = possibly_create_parameter(0.1,
                name='%s - interface_width_protein_solvent' % name)
        self.sld_of_deuterated_protein = possibly_create_parameter(3.23*(10**(-6)),
                name='%s - sld_of_deuterated_protein' % name)
#         self.semi_axis_a = possibly_create_parameter(a,
#                 name='%s - semi_axis_a' % name)
#         self.semi_axis_b = possibly_create_parameter(b,
#                 name='%s - semi_axis_b' % name)
#         self.semi_axis_c = possibly_create_parameter(c,
#                 name='%s - semi_axis_c' % name)
        self.a = 13
        self.b = 16.5
        self.c = 27.5
        self.volume_of_water_molecule = 30
        self.major_axis_length = 55
        self.d = (self.c - self.a)/self.volume_of_water_molecule

    def __call__(self, z, structure):
        self.calculations()
        pass

    @property
    def  parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.interface_air_protein,
                  self.interface_protein_solvent,
                  self.protein_length,
                  self.number_of_water_molecules,
                  self.interface_width_air_solvent,
                  self.interface_width_protein_solvent,
                  self.sld_of_deuterated_protein])
        return p

    def logp(self):
        full = np.arange(150)
        a_p = np.array(list(map(self.area_protein, full)))
        a_w = np.array(list(map(self.area_water, full)))
        max_w = max(a_w)
        max_p = max(a_p)
        max_a = max_w if max_w > max_p else max_p
        return 0

    def calculations(self):
#         a = 13
#         b = 16.5
#         c = 27.5
#         volume_of_water_molecule = 30
        self.e = self.major_axis_length - self.protein_length
        self.delta_a = self.a + (self.d*self.e)
        self.delta_c = self.c - (self.d*self.e)

    def area_protein(self, z):
        first_bit = np.pi*self.delta_a*self.b/(self.delta_c**2)
        second_bit = z - self.interface_air_protein
        final_bit = 2*self.delta_c - second_bit
        return first_bit*second_bit*final_bit

    def area_water(self, z):
        water_length = (self.interface_air_protein 
                        + self.protein_length - self.interface_protein_solvent)
        z_s = self.interface_protein_solvent*0.5*water_length # calc
        first_bit = self.number_of_water_molecules*self.volume_of_water_molecule/2*self.d
        second_bit = (Erf((z, - z_s + 0.5*self.d)/self.interface_width_air_solvent)
                      - Erf((z, - z_s - 0.5*self.d)/self.interface_width_air_solvent))
        return first_bit*second_bit


#
#air = SLD(value=0+0j, name='air')
#----------
#or 
#re =Parameter(0)
#im=Parameter(0)
#air = SLD([re,im])
#---------
#polymer = SLD(1,'polymer')
#silicon = SLD(2,'silicon')
#structure = air(thick=0,rough=0) | polymer(200,4) | silicon(0,3)
# air-polymer roughness of 4, polymer size of 200
#  #Erf()  <-error function
# structure[1].interfaces = Erf() # air-polymer interface
# structure[2].interfaces = Erf()
#need lnprob
#sort vol fractions 
#add areas etc.
#possibly using possibly` creat parameter and parameter p.extend(
#tilt and other values not used either


class bsla_non(Component): #not as in thesis
    def __init__(self, area_per_molecule=1, angle=1, length_of_molecule=1):
        #area per molecule, 
        super(bsla).__init__()
        name='bsla'

    
    #for calculating unphysical thing
    def volume_fraction_protien(self, z):
        pass
    
    def volume_fraction_air(self, z):
        pass
    
    def volume_fraction_water(self, z):
        pass
    
    
class bsla_equation(Componemt): # equation based
    def __init__(self, extent, vs, dz,name='bsla', microslab_max_thickness):
        pass
    
    def __call__(self, z):
        pass
    
    