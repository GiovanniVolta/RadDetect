# Example of loading the YAML configurations and creating an instance of the class
with open('simulation_parameters.yml', 'r') as file:
    config = yaml.safe_load(file)

with open('default_simulation_settings.yml', 'r') as file:
    default_settings = yaml.safe_load(file)

detector_efficiency = MonAlphaDetectorGeometricalEfficiency(config, default_settings)

