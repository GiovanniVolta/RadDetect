import numpy as np

class MonAlphaDetectorGeometricalEfficiency:
    """
    A class to simulate the geometrical efficiency of a monolithic alpha detector.
    This class simulates the detection efficiency of alpha particles emitted from a radioactive source,
    considering the geometrical setup of the source and the detector, as well as the physical properties
    of the emitted alpha particles.

    Attributes:
        PRESSURE_ALPHA_POT (float): Pressure in the chamber.
        DIST_SAMP_DIODE (float): Distance from the sample to the diode.
        ACT_SHAPE (str): Shape of the activity distribution ('point' or 'round').
        ACT_EXT (float): Lateral extent of the activity.
        IMPL_DEPTH (float): Depth of implantation.
        IMPL_SPREAD (float): Spread of the implantation depth.
        IMP_DISTR (str): Distribution of the implantation ('point' or 'gauss').
        DIODE_SHAPE (str): Shape of the diode ('square' or 'round').
        DIODE_SIZE (float): Size of the diode.
        DIODE_RESOLUTION (float): Energy resolution of the diode.
        DIODE_DEAD_LAYER (float): Thickness of the diode's dead layer.
        N_IONS (int): Number of ions to simulate.

    Methods:
        generate_positions(): Generates the initial positions of the alpha particles.
        generate_directions(): Generates the initial directions of the alpha particles.
        calculate_intersections(positions, directions): Calculates the intersections of alpha particles with the detection plane.
        check_diode_hit(intersections): Checks if the alpha particles hit the diode.
        run_simulation(): Runs the simulation and prints the completion message.
    """

    def __init__(self, detector_settings):
        """
        Initializes the detector with configurations for ions and detector settings.
        """
    
        # Assign detector settings
        environment = detector_settings['environment']
        geometry = detector_settings['geometry']
        implantation = detector_settings['implantation']
        diode = detector_settings['diode']
        simulation = detector_settings['simulation']
        
        self.PRESSURE_ALPHA_POT = environment['pressure_alpha_pot']
        self.DIST_SAMP_DIODE = geometry['dist_diode_flange'] - geometry['dist_flange_samp'] + geometry['uncertainty_dist_samp_diode']
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
        """
        Generates the initial positions of the alpha particles based on the implantation depth and distribution.
        """

        # Generate z positions based on the implantation distribution
        if self.IMP_DISTR == 'point':
            pos_z = np.ones((self.N_IONS)) * -1 * (self.IMPL_DEPTH * 1E-6 + self.DIST_SAMP_DIODE)
        elif self.IMP_DISTR == 'gauss':
            pos_z = np.random.normal(loc=-1 * (self.IMPL_DEPTH * 1E-6 + self.DIST_SAMP_DIODE), scale=self.IMPL_SPREAD*1E-6, size=self.N_IONS)
            pos_z[pos_z > -1*self.DIST_SAMP_DIODE] = -1 * self.DIST_SAMP_DIODE

        # Generate (x,y) positions based on the activity shape
        if self.ACT_SHAPE == 'point':
            pos_x = np.zeros((self.N_IONS))
            pos_y = np.zeros((self.N_IONS))
        elif self.ACT_SHAPE == 'round':
            pos_phi = np.random.uniform(0, 2 * np.pi, self.N_IONS)
            pos_r2 = np.random.uniform(0, self.ACT_EXT**2, self.N_IONS)
            pos_r = np.sqrt(pos_r2)
            pos_x = pos_r * np.cos(pos_phi)
            pos_y = pos_r * np.sin(pos_phi)
        elif self.ACT_SHAPE == 'square':
            pos_x = np.random.uniform(-self.ACT_EXT, self.ACT_EXT, self.N_IONS)
            pos_y = np.random.uniform(-self.ACT_EXT, self.ACT_EXT, self.N_IONS)

        positions = np.stack((pos_x, pos_y, pos_z), axis=1)
        return positions

    def generate_directions(self):
        """
        Generates the initial directions of the alpha particles assuming isotropic emission.
        """
        
        # Generate random isotropic directions
        phis = np.random.uniform(0, 2*np.pi, self.N_IONS)
        thetas = np.arccos(np.random.uniform(0, 1, self.N_IONS))

        dir_x = np.sin(thetas) * np.cos(phis)
        dir_y = np.sin(thetas) * np.sin(phis)
        dir_z = np.cos(thetas)
        directions = np.stack((dir_x, dir_y, dir_z), axis=1)
        return directions
    
    def calculate_plane_intersection(self, vec_a, vec_b, vec_c, vec_d, vec_e):
        """
        Calculates the intersection of particle trajectories with a given plane.
        Get the line of flight here \vec{a} is a point on the line (starting points)
        and \vec{b} is the vector giving the emmission direction
        1) Get the equation \vec{x} = \vec{a} + g*\vec{b} describing the line of flight of this ion
        2) Find the intersection of this line with the "ceiling plane" \vec{y} = \vec{c}+h*\vec{d}+i*\vec{e}
        To find the intersection, we have to set \vec{x} = \vec{y} and find the values
        for parameters g, h and i. This is done by solving the system of linear equations
        """

        rhs = vec_c - vec_a
        coefficients = np.stack((vec_b, -vec_d, -vec_e), axis=2)
        solutions = np.linalg.solve(coefficients, rhs)
        return solutions
    
    def calculate_intersections(self, positions, directions):
        """
        Calculates the intersections of the alpha particles with the detection plane.
        Look if the intersection point is inside the diode
        """
        # Ceiling plane center of the plane is cener of diode
        vec_c = np.zeros((self.N_IONS, 3))
        vec_d = np.stack((np.ones((self.N_IONS)), np.zeros((self.N_IONS)), np.zeros((self.N_IONS))), axis=1)
        vec_e = np.stack((np.zeros((self.N_IONS)), np.ones((self.N_IONS)), np.zeros((self.N_IONS))), axis=1)

        solutions = self.calculate_plane_intersection(positions, directions, vec_c, vec_d, vec_e)
        
        # Get the intersection
        intersections = positions + np.array((solutions[:,0], solutions[:,0], solutions[:,0])).T * directions
        return intersections

    def check_diode_hit(self, intersections):
        """
        Checks if the alpha particles hit the diode based on the diode shape and size.
        """
            
        if self.DIODE_SHAPE == 'square':
            bool_diode_hit = np.logical_and(np.less(np.abs(intersections[:,0]), self.DIODE_SIZE/2.), np.less(np.abs(intersections[:,1]), self.DIODE_SIZE/2.))
        elif self.DIODE_SHAPE == 'round':
            bool_diode_hit = np.less(np.sqrt(intersections[:,0]**2 + intersections[:,1]**2), self.DIODE_SIZE/2.)
        return bool_diode_hit

    def calculate_detection_efficiency(self, bool_diode_hit):
        """
        Calculates the detection efficiency based on the number of hits.
        """
        detection_efficiency = np.sum(bool_diode_hit) / self.N_IONS
        return detection_efficiency

    def calculate_detection_efficiency_stat_uncertainty(self, bool_diode_hit):
        """
        Calculates the statistical uncertainty of the 
        detection efficiency based on the number of hits.
        """
        detection_efficiency_stat_uncert = np.sqrt(np.sum(bool_diode_hit)) / self.N_IONS
        return detection_efficiency_stat_uncert

    def run_simulation(self,):
        self.positions = self.generate_positions()
        self.directions = self.generate_directions()
        self.intersections = self.calculate_intersections(self.positions, self.directions)
        self.bool_diode_hit = self.check_diode_hit(self.intersections) 
        self.detection_efficiency = self.calculate_detection_efficiency(self.bool_diode_hit)
        self.detection_efficiency_stat_uncert = self.calculate_detection_efficiency_stat_uncertainty(self.bool_diode_hit)
        print(f"Disance Sample - Diode: {self.DIST_SAMP_DIODE} mm")
        print(f"Detection Efficiency: {self.detection_efficiency * 100:.2f}%")
        print(f"Statistical Uncertainty: {self.detection_efficiency_stat_uncert * 100:.2f}%")
        print("Simulation complete.")