import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

class MonAlphaDetectorGeometricalEfficiency:
    def __init__(self, ion_config, detector_settings):
        # Unpack the configuration dictionary
        energies = ion_config['energies']
        branching_ratios = ion_config['branching_ratios']
        srim_losses = ion_config['srim_losses']

        # Assign the parameters from the configuration
        self.ENERGY_1_226RA = energies['energy_1_226ra']
        self.ENERGY_2_226RA = energies['energy_2_226ra']
        self.BR_RATIO_1 = branching_ratios['br_ratio_1']
        self.BR_RATIO_2 = branching_ratios['br_ratio_2']
        self.SRIM_LOSS_226RA = srim_losses['srim_loss_226ra']
        self.SRIM_LOSS_226RA_SI = srim_losses['srim_loss_226ra_si']
        self.SRIM_LOSS_226RA_AIR = srim_losses['srim_loss_226ra_air']

        # Assign default settings
        environment = detector_settings['environment']
        geometry = detector_settings['geometry']
        implantation = detector_settings['implantation']
        diode = detector_settings['diode']
        simulation = detector_settings['simulation']
        self.PRESSURE_ALPHA_POT = environment['pressure_alpha_pot']
        self.DIST_SAMP_DIODE = geometry['dist_samp_diode']
        self.ACT_SHAPE = implantation['act_shape']
        self.ACT_EXT = implantation['act_ext']
        self.IMPL_DEPTH = implantation['impl_depth']
        self.IMPL_SPREAD = implantation['impl_spread']
        self.IMP_DISTR = implantation['imp_distr']
        self.DIODE_SHAPE = diode['diode_shape']
        self.DIODE_SIZE = diode['diode_size']
        self.DIODE_RESOLUTION = diode['diode_resolution']
        self.DIODE_DEAD_LAYER = diode['diode_dead_layer']
        self.N_IONS = simulation['n_ions']

    def generate_positions(self):
        if self.IMP_DISTR == 'point':
            pos_z = np.ones((self.N_IONS)) * -1*(self.IMPL_DEPTH * 1E-6 + self.DIST_SAMP_DIODE)
        elif self.IMP_DISTR == 'gauss':
            pos_z = np.random.normal(loc=-1*(self.IMPL_DEPTH * 1E-6 + self.DIST_SAMP_DIODE), scale=self.IMPL_SPREAD*1E-6, size=self.N_IONS)
            pos_z[pos_z > -1*self.DIST_SAMP_DIODE] = -1*self.DIST_SAMP_DIODE

        if self.ACT_SHAPE == 'point':
            pos_x = np.zeros((self.N_IONS))
            pos_y = np.zeros((self.N_IONS))
        elif self.ACT_SHAPE == 'round':
            pos_phi = np.random.uniform(0, 2*np.pi, self.N_IONS)
            pos_r2 = np.random.uniform(0, self.ACT_EXT**2, self.N_IONS)
            pos_r = np.sqrt(pos_r2)
            pos_x = pos_r*np.cos(pos_phi)
            pos_y = pos_r*np.sin(pos_phi)

        positions = np.stack((pos_x, pos_y, pos_z), axis=1)
        return positions

    def generate_directions(self):
        phis = np.random.uniform(0, 2*np.pi, self.N_IONS)
        thetas = np.arccos(np.random.uniform(0, 1, self.N_IONS))

        dir_x = np.sin(thetas) * np.cos(phis)
        dir_y = np.sin(thetas) * np.sin(phis)
        dir_z = np.cos(thetas)
        directions = np.stack((dir_x, dir_y, dir_z), axis=1)
        return directions

    def calculate_intersections(self, positions, directions):
        vec_a = positions
        vec_b = directions

        vec_c = np.zeros((self.N_IONS, 3))
        vec_d = np.stack((np.ones((self.N_IONS)), np.zeros((self.N_IONS)), np.zeros((self.N_IONS))), axis=1)
        vec_e = np.stack((np.zeros((self.N_IONS)), np.ones((self.N_IONS)), np.zeros((self.N_IONS))), axis=1)

        rhs = vec_c - vec_a
        coefficients = np.stack((vec_b, -vec_d, -vec_e), axis=2)
        solutions = np.linalg.solve(coefficients, rhs)

        intersections = vec_a + np.array((solutions[:,0], solutions[:,0], solutions[:,0])).T * vec_b
        return intersections

    def check_diode_hit(self, intersections):
        if self.DIODE_SHAPE == 'square':
            bool_diode_hit = np.logical_and(np.less(np.abs(intersections[:,0]), self.DIODE_SIZE/2.), np.less(np.abs(intersections[:,1]), self.DIODE_SIZE/2.))
        elif self.DIODE_SHAPE == 'round':
            bool_diode_hit = np.less(np.sqrt(intersections[:,0]**2 + intersections[:,1]**2), self.DIODE_SIZE/2.)
        return bool_diode_hit

    def run_simulation(self):
        positions = self.generate_positions()
        directions = self.generate_directions()
        intersections = self.calculate_intersections(positions, directions)
        bool_diode_hit = self.check_diode_hit(intersections)

        # Additional calculations for distance, energy, etc., can be added here following the same pattern

        print("Simulation complete.")