import numpy as np
import matplotlib.pyplot as plt
from detran import *
import sys
import time
# Input
inp = InputDB.Create()
inp.put_str("problem_type",    "fixed")
inp.put_str("inner_solver",    "SI")
inp.put_int("inner_max_iters", 1000)
inp.put_dbl("inner_tolerance", 1e-12)
inp.put_int("inner_print_out", 1)
inp.put_str("outer_solver",    "GS")
inp.put_int("outer_max_iters", 100)
inp.put_dbl("outer_tolerance", 1e-8)
inp.put_int("outer_print_out", 0)
inp.put_str("bc_right",         "reflect")
inp.put_int("quad_order",      32)
inp.put_int("store_angular_flux",       1)
# Material
mat = Material.Create(1, 1, False)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.finalize()
# Mesh
fm = [6]
cm = [0, 4]
mesh = Mesh1D.Create(fm, cm, [0])
# Execute 0.99999945  0.99999288  0.99990083  0.99861878  0.98076211  0.73205081
#         1.99927621  1.9975391   1.98989584  1.95800257  1.82531549  1.27338559
execute = Execute1D(sys.argv)
execute.initialize(inp, mat, mesh)
execute.set_external_source(ConstantSource.Create(mesh, execute.get_quadrature(), 1, 1.0))
execute.solve()
# Postprocess
state = execute.get_state()
print np.asarray(state.phi(0)) 
#print np.asarray(state.psi(0, 0, 0))
#print np.asarray(state.psi(0, 1, 0))
plt.plot(np.asarray(state.phi(0)), linewidth=3)
plt.grid(True)
#plt.show()
execute.finalize()
