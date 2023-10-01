def find_coefficients(expr, terms):
  """For each of the terms, find its coefficient in expr,
  by setting all other terms to zero."""

  coefs = {}
  for symbol in terms:
    coef = expr.subs(symbol, 1)
    for other in terms:
      if symbol != coef:
        coef = coef.subs(other, 0)
    coefs[symbol] = coef
  return coefs
