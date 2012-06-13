def fill_input(inp) :
  inp.put_int("number_groups",      7)
  inp.put_str("equation",           "dd")
  inp.put_int("inner_max_iters",    1)
  inp.put_dbl("inner_tolerance",    1e-7)
  inp.put_int("inner_print_out",    0)
  inp.put_int("outer_max_iters",    1)
  inp.put_dbl("outer_tolerance",    1e-7)
  inp.put_int("outer_print_out",    0)
  inp.put_int("eigen_max_iters",    10000)
  inp.put_dbl("eigen_tolerance",    1e-12)
  inp.put_str("bc_left",            "vacuum")
  inp.put_str("bc_right",           "reflect")
  inp.put_str("bc_bottom",          "reflect")
  inp.put_str("bc_top",             "reflect")
  inp.put_int("store_angular_flux", 1)

def get_geom() :
  cm = [0.0, 0.5, 1.5,  2.0]
  fm = [5, 10, 5]
  mt = [1,1,1,1,0,1,1,1,1]
  return cm, fm, mt
