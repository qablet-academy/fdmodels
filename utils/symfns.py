def find_coefficients(expr, terms):
  """For each of the terms, find its coefficient in expr,
  by setting all other terms to zero."""

  coefs = {}
  for symbol in terms:
    coef = expr.subs(symbol, 1)
    for other in terms:
      if other != symbol:
        coef = coef.subs(other, 0)
    coefs[symbol] = coef
  return coefs


def discretize_crank_nicolson(model_expr, V, dV, Vx, Vxx, rowtype='normal'):
  """Substitute V, Vx, Vxx by crank nicolson discretization."""

    # Define the gridpoints
    # Vu_l, Vm_l, Vd_l are the upper, middle and lower gridpoints on the left hand side.
    # and Vu_r, Vm_r, Vd_r are the corresponding gridpoints on the right hand side. 
    Vu_l, Vm_l, Vd_l, Vu_r, Vm_r, Vd_r = symbols("V^u_l V^m_l, V^d_l, V^u_r, V^m_r, V^d_r")

    # Substitute V, dV 
    expr = model_expr.subs(V, (Vm_l + Vm_r)/2)
    expr = expr.subs(dV, (Vm_r - Vm_l))

    # Substitute Vx, Vxx
    if rowtype == "up":  # There is no Vu at upper boundary
        expr = expr.subs(Vx, ((Vm_l - Vd_l)/2 + (Vm_r - Vd_r)/2)/(dx))
        expr = expr.subs(Vxx, 0)
    if rowtype == "dn":  # There is no Vd at lower boundary
        expr = expr.subs(Vx, ((Vu_l - Vm_l)/2 + (Vu_r - Vm_r)/2)/(dx))
        expr = expr.subs(Vxx, 0)
    else:  # normal case
        expr = expr.subs(Vx, ((Vu_l - Vd_l)/2 + (Vu_r - Vd_r)/2)/(2 * dx))
        expr = expr.subs(Vxx, ((Vu_l + Vd_l - 2*Vm_l) + (Vu_r + Vd_r - 2 * Vm_r))/(2 * dx ** 2))

    return expr, [Vu_l, Vm_l, Vd_l], [Vu_r, Vm_r, Vd_r]
