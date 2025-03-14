from QV import calculateUtilities, formatList
import sympy

n = 3
p = [sympy.symbols(f'p{i}',negative=False,real=True) for i in range(n)]
u = [sympy.symbols(f'u{i}',negative=False,real=True) for i in range(n)]
possible_votes = list(range(n+6))

pv = [0.9,0.8,0.6]
uv = [400,350,50]

d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv + uv, texFile = 'tex/LVworse.tex' )

pv = [0.979,0.965,0.953]
uv = [250,65,98]
#{p0: 0.9794544064507311, p1: 0.9654849465351896, p2: 0.9530368061677184, u0: 249, u1: 65, u2: 98}

d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv, texFile = 'tex/QVworse.tex')

pv = [.78,.74,.72]
uv = 3*[78]
d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv )

pv = [.7,.67,.6]
uv = 3*[100]
d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv )

pv = [.9,.8,.6]
uv = 3*[100]
d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv )

pv = [.9,.78,.51]
uv = 3*[100]
d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv )

pv = [.97,.95,.88]
uv = 3*[125]
d = calculateUtilities( p, u, pv, uv )
formatList( d, params = pv+uv, texFile='tex/QVworse_uniform.tex' )