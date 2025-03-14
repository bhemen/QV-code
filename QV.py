#!/usr/bin/env python

import sympy
import itertools
import json
import numpy as np
import time


fLV = lambda x: x      #Linear voting cost function
fQV = lambda x: x**2   #Quadratic voting cost function

def zeros( shape ):
    """
        shape - list of positive integers
        Returns a *mutable* sympy array of 0s with the given shape
    """
    return sympy.MutableDenseNDimArray( [0 for _ in range( sympy.prod(shape) ) ], shape=shape )

def indices(A):
    """
        A - n dimensional array
        Returns a list of all the indices in A
    """
    return list( itertools.product( *[list( range(l)) for l in A.shape] ) )

def winning_sets(v):
    """
        Calculate the winning sets for a vote vector, v
    """
    F = []
    signs = list(itertools.product( [-1,1], repeat=len(v) ) )
    for s in signs:
        if sum( a*b for a,b in zip(v,s) ) > 0:
            F.append( tuple(i for i,e in enumerate(s) if e == 1 ) )
    return F

def success_prob(p, v ):
    """
        Given a set of probabilities p (of length n) and a set of votes, v, calculate the probability of success
        Return a scalar winning probability
    """
    prob = 0
    ws = winning_sets(v)
    for w in ws:
        prob += sympy.prod( [ pp if i in w else 1-pp for i,pp in enumerate(p) ] )
    return prob

def winning_probs( p, possible_votes ):
    """
        p - vector of success probabilities
        possible_votes - set of possible votes
        Returns an n-dimensional array, indexed by possible votes, whose entries are the success probabilities
    """
    n = len(p)
    P = zeros( n*[len(possible_votes)] )
    for v in indices(P):
        P[v] = sympy.simplify( success_prob( p, v ) )
    return P
    
def payoff_matrices( p, possible_votes, u, f, W=None):
    """
        p - vector of success probabilities
        possible_votes - set of possible votes
        u - vector of utilities
        f - cost function (e.g. LV or QV)
        Returns a list of n n-dimensional arrays.  The i'th n-dimensional array gives the payoff for voter i for all possible combinations of votes
    """
    n = len(p)
    if W is None:   #If you're doing repeated calculations with the same set of parametrs, it can save time calculate W once, rather than calculate it every time you call this function
        W = winning_probs( p, possible_votes )
    PMs = [0 for _ in range(n)]
    for i in range(n):
        P = u[i]*W
        for v in indices(P):
            P[v] -= f(v[i])
        PMs[i] = P
    return PMs   

def best_response( pM, i ):
    """
        pM = payoff matrix for player i
        i = player index
        returns a best_response matrix br
        br.shape = pM.shape
        br[x] = 1 if and only if x is player i's best response
    """
    C = indices(pM)
    br = zeros( shape=pM.shape )
    for c in C:
        sl = tuple( list(c[:i]) + [j] + list(c[i+1:]) for j in range( pM.shape[i] ) )
        #print( f"sl = {sl}" )
        try:
            m = max( pM[cc] for cc in sl )
        except Exception as e:
            print( e )
            print( list( pM[cc] for cc in sl ) )
        if pM[c] == m:
            br[c] = 1
    return br

def equilibria( BRs ):
    """
        BRs - list of n best-response matrices (these must be numeric, not symbolic)
        returns a list of equilibria
    """
    n = len(BRs)
    assert all( len( BR.shape ) == n for BR in BRs ), 'Error, BRs must have the same number of dimensions'
    C = indices(BRs[0])
    eqs = set()
    for c in C:
        if all( BRs[i][c] == 1 for i in range(n) ):
            eqs.add(c)
    return eqs

def best_responses( PMs, value_dict ):
    n = len(PMs)
    assert all( len( PM.shape ) == n for PM in PMs ), 'Error, PMs must have the same number of dimensions'
    BRs = [0 for _ in range(n)]
    for i in range(n):
        C = indices( PMs[i] )
        P = zeros( PMs[i].shape )
        for c in C:
            P[c] = PMs[i][c].subs(value_dict)
        BRs[i] = best_response( P, i )
    return BRs

def netU( p, u, e, value_dict, W=None):
    """
        p - list of success probabilities (symbolic)
        u - list of utilities (symbolic)
        e - list of equilibria
        value_dict - instantiations for the p_i and u_i
        Returns - a list of net utilities (of same length as e)
    """
    if W is None:
        W = winning_probs( p, possible_votes )
    netu = 0
    for i,ei in enumerate(e):
        netu += (u[i]*W[e]).subs(value_dict) #Note that net utility is independent of cost function because voting costs are redistributed
    return netu
    
def searchSpace(n,outfile,uniformUtilities=False):
    """
        Search the space for equilibria in an n-voter game, 
        where the utilities in LV and QV are different.

        n - number of players
        uniformUtilities - if True, utilities are uniform
    """
    assert isinstance(n,int), 'Error, n must be an integer'
    assert n > 0, 'Error, n must be positive'

    p = [sympy.symbols(f'p{i}',negative=False,real=True) for i in range(n)]
    u = [sympy.symbols(f'u{i}',negative=False,real=True) for i in range(n)]
    possible_votes = list(range(n+6))
    while True:
        """
            Search for examples where LV and QV have equilibria with different utilities
        """

        n = len(p)
        pv = np.random.uniform( low=.5, high=1, size=n )    #Success probabilities are uniform [.5,1]
        pv = sorted( pv, reverse=True )                     #Success probabilities are sorted in decreasing order
        pv = list( np.round( pv, decimals=3 ) )
#        pv[1] = pv[0] - .01
#        pv = n*[0]
#        pv[0] = np.random.uniform( low=.5, high=.6 )
#        pv[1] = np.random.uniform( low=.5, high=pv[0] )
#        pv[2] = pv[1]
        if uniformUtilities:
            uv = n*[np.random.randint(low=1,high=10)]
        else:
            uv = np.random.randint( low=50, high=500, size=n )  #Utilities are random 50 - 100
        assert len(uv) == n
        value_dict = { k:v for k,v in zip(p,pv) }
        value_dict.update( {k:v for k,v in zip(u,uv) } )
        str_dict = { k: str(v) for k,v in value_dict.items() }
        print( str_dict )                                 #Print something so we know the for loop is running
        W = winning_probs( p, possible_votes )
        LVPMs = payoff_matrices( p, possible_votes, u, fLV, W=W)
        QVPMs = payoff_matrices( p, possible_votes, u, fQV, W=W)
        LVBRs = best_responses( LVPMs, value_dict )
        QVBRs = best_responses( QVPMs, value_dict )
        LVE = equilibria( LVBRs )
        QVE = equilibria( QVBRs )
        CommonE = LVE.intersection( LVE, QVE )                 #Equilibria under both LV and QV
        OnlyLVE = LVE.difference( CommonE )                    #Equilibria that are only under LV
        OnlyQVE = QVE.difference( CommonE )                    #Equilibria that are only under QV
        LVUs = [ netU( p,u, eq, value_dict, W=W) for eq in LVE ]   #LV equilibrium Utilities
        QVUs = [ netU( p,u, eq, value_dict, W=W) for eq in QVE ]   #QV equilibrium utilities
        maxQVU = max(QVUs,default=0)
        maxLVU = max(LVUs,default=0)
        minQVU = min(QVUs,default=0)
        minLVU = min(LVUs,default=0)
        
        if len(OnlyLVE) > 0 and len(OnlyQVE) > 0:
            print( "=================" )
            print( "Both different" )
            if maxQVU > maxLVU:
                print( f"Better QVE" )
                writeParams( pv, uv, "Better QVE" )
            if maxLVU > maxQVU:
                print( f"Better LVE" )
                writeParams( pv, uv, "Better LVE" )
            if minQVU < minLVU:
                print( f"Worse QVE" )
                writeParams( pv, uv, "Worse QVE" )
            if minLVU < minQVU:
                print( f"Worse LVE" )
                writeParams( pv, uv, "Worse LVE" )
            print( f"max(QVE) = {maxQVU}" )
            print( f"min(QVE) = {minQVU}" )
            print( f"max(LVE) = {maxLVU}" )
            print( f"min(LVE) = {minLVU}" )
            print( "=================" )
            continue
        if len(OnlyLVE) > 0:
            print( "=================" )
            print( "Extra LVE" )
            if maxLVU > maxQVU:
                print( f"Better LVE" )
                writeParams( pv, uv, "Better LVE" )
            if minLVU < minQVU:
                print( f"Worse LVE" )
                writeParams( pv, uv, "Worse LVE" )
            print( f"max(QVE) = {maxQVU}" )
            print( f"min(QVE) = {minQVU}" )
            print( f"max(LVE) = {maxLVU}" )
            print( f"min(LVE) = {minLVU}" )
            print( "=================" )
            continue
        if len(OnlyQVE) > 0:
            print( "=================" )
            print( "Extra QVE" )
            if maxQVU > maxLVU:
                print( f"Better QVE" )
                writeParams( pv, uv, "Better QVE" )
            if minQVU < minQVU:
                print( f"Worse QVE" )
                writeParams( pv, uv, "Worse QVE" )
            print( f"max(QVE) = {maxQVU}" )
            print( f"min(QVE) = {minQVU}" )
            print( f"max(LVE) = {maxLVU}" )
            print( f"min(LVE) = {minLVU}" )
            print( "=================" )
            continue

def writeParams( pv, uv, tag ):
    with open( outfile, 'a') as f:
        f.write( ",".join( [str(x) for x in pv+uv] + [tag] ) + "\n" )


def calculateUtilities( p, u, pv, uv ):
    """
        p - vector of symbols representing success probabilities (e.g. p0,p1,...)
        u - vector of symbols representing utilities (e.g. u0,u1,...)
        pv - probability values
        uv - utility values
        Returns Dictionary with equilibria and utilities under QV and LV
    """
    assert len(p) == len(u), 'Error, len(p) != len(u)'
    assert len(pv) == len(p), 'Error, len(p) != len(pv)'
    assert len(uv) == len(u), 'Error, len(u) != len(uv)'
    value_dict = { k:v for k,v in zip(p,pv) }
    value_dict.update( {k:v for k,v in zip(u,uv) } )
    start_time = time.time()
    W = winning_probs( p, possible_votes )
    end_time = time.time()
    print( f"Calculated winning probabilities ({round(end_time-start_time)}s)" )
    start_time = time.time()
    LVPMs = payoff_matrices( p, possible_votes, u, fLV, W=W )
    end_time = time.time()
    print( f"Calculated LV payoff matrices ({round(end_time-start_time)}s)" )
    start_time = time.time()
    QVPMs = payoff_matrices( p, possible_votes, u, fQV, W=W)
    end_time = time.time()
    print( f"Calculated QV payoff matrices ({round(end_time-start_time)}s)" )
    start_time = time.time()
    LVBRs = best_responses( LVPMs, value_dict )
    end_time = time.time()
    print( f"Calculated LV BR matrices ({round(end_time-start_time)}s)" )
    start_time = time.time()
    QVBRs = best_responses( QVPMs, value_dict )
    end_time = time.time()
    print( f"Calculated QV BR matrices ({round(end_time-start_time)}s)" )
    start_time = time.time()
    LVE = equilibria( LVBRs )
    end_time = time.time()
    print( f"Found {len(LVE)} LV Equilibria ({round(end_time-start_time)}s)" )
    start_time = time.time()
    QVE = equilibria( QVBRs )
    end_time = time.time()
    print( f"Calculated {len(QVE)} QV Equilibria ({round(end_time-start_time)}s)" )
    allE = LVE.union( QVE )
    CommonE = LVE.intersection( LVE, QVE )
    OnlyLVE = LVE.difference( CommonE )
    OnlyQVE = QVE.difference( CommonE )

    print( f"{len(CommonE)} common equilibria" )
    print( f"{len(OnlyLVE)} LV-only equilibria" )
    print( f"{len(OnlyQVE)} QV-only equilibria" )

    start_time = time.time()
    EUs = [ (list(eq), 'Both', float(netU( p,u, eq, value_dict, W=W ))) for eq in CommonE] #Calculate net utility for equilibria
    EUs += [ (list(eq), 'LV', float(netU( p,u, eq, value_dict, W=W ))) for eq in OnlyLVE]
    EUs += [ (list(eq), 'QV', float(netU( p,u, eq, value_dict, W=W ))) for eq in OnlyQVE]
    end_time = time.time()
    print( f"Calculated net utilities ({round(end_time-start_time)}s)" )

    EUs = sorted( EUs, key=lambda x: -x[2] ) #Sort equilibria in decreasing order based on net Utility
    print( value_dict )
    return EUs

def formatList( L, params=None, texFile=None ):
    """
        L - list of equilibria
        params - list of length 2n [p0,p1,p2,u0,u1,u2] 
        texFile - filename

        Take a list of equilibria and display them in a printable format
        Then format them for LaTeX
    """
    for x in L:
        print( f"{x[0]} : {x[1]} : {round(x[2],2)}" )

    if texFile:
        print( f"Writing to {texFile}" )
        with open( texFile, 'w' ) as f:
            f.write( '\\begin{center}\n' )
            f.write( '\\begin{tabular}{c|c|c}\n' )
            f.write( 'Strategy & Mechanism & Net Utility \\\\\n' )
            f.write( '\\hline\n' )
            for x in L:
                f.write( f'${str(x[0]).replace("[","(").replace("]",")")}$ & {x[1]} & ${round(x[2],2)}$ \\\\\n' )
            f.write( '\\end{tabular}\n' )
            f.write( '\n\n' )
            if params:
                f.write ( '\\smallskip\n' )
                f.write( "$" + ",".join( [f'p_{i+1} = {params[i]}' for i in range(n)] ) + "$\n" )
                f.write( "$" + ",".join( [f'u_{i-n+1} = {params[i]}' for i in range(n,2*n)] ) + "$\n" )
            f.write( '\\end{center}\n' )
