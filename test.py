from src.Pitzer.models import FluidPitzer


species = {
    'Na+': 1, # always be 1 if Na is the internal standard
    'K+': 2 ,  # K/Na = 2
    # 'Mg+2':0.1
}

# create a fluid object with information from microthemometric and LA-ICP-MS data
fluid = FluidPitzer(
    # the initial guess
    x0=(1,1),
  
    # species defined before
    species=species,
  
    # the last melting solid phase
    equilibrium='KCl(s) = K+(aq) + Cl-(aq)',
  
    # melting temperature of the last-melting solid phase, [°C]
    t = 25,
)

result = fluid.optimize()

print(result)

