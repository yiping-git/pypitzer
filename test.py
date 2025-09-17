from src.Pitzer.models import FluidPitzer

# aqueous species determined in LA-ICP-MS analysis
species = {
    'Na+': 1, # always be 1 if Na is the internal standard
    'K+': 2,  # K/Na = 2
}

# create a fluid object with information from microthemometric and LA-ICP-MS data
fluid = FluidPitzer(
    # the initial guess
    x0=(3, 3),
  
    # species defined before
    species=species,
  
    # the last melting solid
    solids=['KCl'],
  
    # melting temperature of the last solid, °C
    t = 25,
)

result = fluid.optimize()

print(result)