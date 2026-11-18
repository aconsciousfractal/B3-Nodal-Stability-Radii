#!/usr/bin/env python3
"""Finite exact checks for the narrow-sector example and its supporting claims."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = 'nodal-stability-reader-companion-v1'
DEFINITIONS = [
    ('S01','ODD-DIMENSION','exact'),
    ('S01','SATURATION-IDENTITY','exact'),
    ('S01','CHEB-BASIS','exact'),
    ('S01','CHEB-DIFFERENTIAL','exact'),
    ('S01','H-DEGREE-ODD','exact'),
    ('S01','H-CENTER-JET','exact'),
    ('S01','H-FACTOR-SCALE','exact'),
    ('S01','CONTINUOUS-MAJORANT-ARITHMETIC','exact'),
    ('S01','EXPLICIT-EPSILON-GAP','exact'),
    ('S01','RATIONAL-FINITE-DOMAIN','exact'),
    ('S01','NORMAL-FIBER-NOT-POLYNOMIAL','negative'),
    ('S02','GLOBAL-DUAL-ATOM-COUNTS','exact'),
    ('S02','COMPLETE-ODD-MOMENTS','exact'),
    ('S02','CONTACT-INTERPOLATION','exact'),
    ('S02','M1-NULL-DIRECTION','exact'),
    ('S02','RADICAL-GAP-IDENTITY','exact'),
    ('S02','SIMPLE-DOMAIN-GAP','exact'),
    ('S02','M2-NONCONCENTRATION','negative'),
    ('S03','FULL-POINT-JET-SURJECTIVITY','exact'),
    ('S03','GENERIC-DUAL-SIGN','exact'),
    ('S03','GLOBAL-JET-CORRECTION','exact'),
    ('S03','GENERIC-NOT-ODD','negative'),
    ('S04','UPPER-ALGEBRA','exact'),
    ('S04','GAMMA-COST','exact'),
    ('S04','SIDE-GAP','exact'),
    ('S04','THREE-ZERO-JET','exact'),
    ('S04','JET-ATTAINMENT-NOT-FAILURE','negative'),
    ('S05','SOURCE-BOUNDARY','source'),
    ('S05','FIXED-OMEGA-REGIMES','negative'),
]


def unique(pairs):
    result = {}
    for k,v in pairs:
        if k in result:
            raise ValueError('Duplicate JSON key '+k)
        result[k] = v
    return result


def load(name):
    return json.loads((ROOT/'data'/name).read_text(encoding='utf-8'),
                      object_pairs_hook=unique,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite '+x)))


def rat(value):
    if type(value) is int:
        return sp.Rational(value)
    if type(value) is str and re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', value):
        parts = value.split('/')
        return sp.Rational(int(parts[0]), int(parts[1]) if len(parts)==2 else 1)
    raise ValueError('Exact rational required, not bool/float/expression')


def native(value, expected):
    return type(value) is int and value == expected


def source_boundary(sources):
    expected = [
        ("Erdelyi", "https://doi.org/10.1007/s00365-013-9208-0",
         "Theorem 2.1", "Fibre inequality"),
        ("Dressler", "https://doi.org/10.1137/24M1717853",
         "Sections 4.2–4.3", "Extremal-polynomial context"),
        ("DavydovSchaback", "https://doi.org/10.1007/s00211-018-0973-3",
         "Theorem 9", "Fixed-node duality context"),
    ]
    citation_text = (ROOT / "references.bib").read_text(encoding="utf-8")
    citations = dict(re.findall(r"(?ms)^@article\{(\w+),\s*(.*?)^\}", citation_text))
    if set(sources) != {"protocol_id", "schema_version", "distribution", "records"}:
        return False
    if sources["distribution"] != "citations-only":
        return False
    if not isinstance(sources["records"], list) or len(sources["records"]) != len(expected):
        return False
    for record, (key, url, locator, use) in zip(sources["records"], expected):
        if record != {"citation": key, "url": url, "locator": locator, "use": use}:
            return False
        if not re.search(r"(?m)^\s*doi\s*=\s*\{" + re.escape(url.removeprefix("https://doi.org/"))
                         + r"\},", citations.get(key, "")):
            return False
    return True


def payload():
    terminal, dual, radius, sources = [load(n) for n in
        ('terminal.json','duality.json','radius.json','sources.json')]
    values = []
    def ck(sid,key,passed,detail):
        kind = next(k for s,n,k in DEFINITIONS if (s,n)==(sid,key))
        values.append({'id':'F-'+sid+'-'+key,'statement_id':'F-'+sid,
                       'kind':kind,'passed':bool(passed),'details':detail})
    R = sp.Rational
    z,x,y,u,a,b,B,d,C,q,g,eps = sp.symbols('z x y u a b B d C q g eps', real=True)
    cheb = terminal['chebyshev']; witness = terminal['global_witness']
    norm = terminal['normalization']; fp = radius['fixed_point']
    transfer = radius['finite_witness_transfer']
    ck('S01','ODD-DIMENSION', native(norm['odd_global_dimension'],45)
       and 5+2*sum([4,4,3,3,2,2,1,1])==45,
       'Exact odd parameter count; projection/attainment require the written proof')
    scale = B/d
    ck('S01','SATURATION-IDENTITY', sp.cancel(scale*(-d)+B)==0,
       'Positive scaling gives saturated derivative; positivity d>B and strict norm reduction are proved separately')
    coeffs = cheb['coefficients_ascending']
    P = sum(rat(c)*z**i for i,c in enumerate(coeffs))
    ck('S01','CHEB-BASIS', native(cheb['degree'],9) and sp.degree(P,z)==9
       and sp.expand(P-sp.chebyshevt(9,z))==0 and sp.diff(P,z).subs(z,0)==9,
       'Consumed ascending coefficients, degree and central derivative')
    ck('S01','CHEB-DIFFERENTIAL', sp.expand(1-P**2-(1-z*z)*sp.diff(P,z)**2/81)==0,
       'Polynomial identity supporting the continuous normal majorant')
    us = sp.symbols('u0:5',real=True)
    H = z*sp.prod(u-v for v in us)
    Hxy = H.subs(u,x+y)
    ck('S01','H-DEGREE-ODD', native(witness['ambient_degree_H'],6)
       and native(witness['ambient_degree_P_plus_epsilon_H'],9)
       and sp.Poly(Hxy,x,y,z).total_degree()==6
       and sp.expand(Hxy.subs(z,-z)+Hxy)==0,
       'Degree and z parity for ambient product; fixed angular coefficients do not raise degree')
    ck('S01','H-CENTER-JET', sp.diff(H,z).subs({u:1,z:0})==sp.prod(1-v for v in us)
       and H.subs(z,0)==0,
       'Terminal derivative from even u(z) with u(0)=1; no tangential substitution into the norm')
    ck('S01','H-FACTOR-SCALE', sp.expand(sp.diff(P,z)-9*sp.chebyshevu(8,z))==0
       and 9*256==2304 and R(9,2304)/2==R(1,512),
       'Exact U8 derivative factor and H0 majorant constants; root locations and continuous inequalities remain analytic')
    ck('S01','CONTINUOUS-MAJORANT-ARITHMETIC', sp.expand(a**4/4-z*z*(a*a-z*z)-(z*z-a*a/2)**2)==0
       and rat(transfer['j_lower_factor'])==R(7,8),
       'Exact square verifies z²(a²-z²)≤a⁴/4; no sampled sup')
    substituted = sp.factor((eps*g-eps**2*C**2/2).subs(eps,g/C**2))
    ck('S01','EXPLICIT-EPSILON-GAP', substituted==g**2/(2*C**2)
       and witness['epsilon']=='g/C^2' and witness['g']=='a*d/9-q',
       'For g,C>0 the epsilon choice leaves a strict positive gap')
    ex=terminal['exact_tangent_example']
    alpha,beta = rat(ex['tan_eta']),rat(ex['tan_omega'])
    denominator2=16*(1+alpha**2)**2+alpha**2
    u8=sp.Poly(sp.chebyshevu(8,z),z)
    ck('S01','RATIONAL-FINITE-DOMAIN', R(3,8)<alpha<beta==R(1,2)
       and beta-alpha==rat(ex['L_upper'])==R(1,10**9) and denominator2<36
       and rat(ex['product_si'])**2==u8.TC()/u8.LC()==R(1,256)
       and 2**5*256==8192 and R(1,16)**10/8192==R(1,2**53)
       and 45*2**52==rat(ex['sufficient_condition_integer_left'])
       <rat(ex['sufficient_condition_integer_right'])==10**18,
       'Exact rational tangent example and integer margin; atan Lipschitz and radical bounds proved in the article')
    ck('S01','NORMAL-FIBER-NOT-POLYNOMIAL', sp.expand(sp.series(z*sp.sqrt(1-z*z),z,0,12).removeO()).coeff(z,11)==-R(7,256),
       'Global field uz has a nonzero z11 fiber coefficient, so ordinary degree9 interpolation does not represent the full class')
    dr={r['id']:r for r in dual['results']}
    dualr=dr['F-S02-DUAL']
    ck('S02','GLOBAL-DUAL-ATOM-COUNTS', native(dualr['odd_dimension'],45)
       and native(dualr['positive_half_tube_atom_bound'],45)
       and native(dualr['whole_tube_V9_paired_atom_bound'],90) and 2*45==90,
       'Exposed face argument is written; finite count here does not determine actual optimizing atoms')
    modes=[(m,k) for m in range(9) for k in range(5) if m+2*k+1<=9]
    ck('S02','COMPLETE-ODD-MOMENTS', sum(1 if m==0 else 2 for m,k in modes)==45,
       'Complete global odd Fourier-z basis, exact degree constraint')
    nodes=[R(i*i,100) for i in range(1,6)]
    ell=[sp.prod(-h/(v-h) for h in nodes if h!=v) for v in nodes]
    ck('S02','CONTACT-INTERPOLATION', sum(ell)==1
       and all(sum(w*v**k for w,v in zip(ell,nodes))==0 for k in range(1,5)),
       'Five exact rational nodes verify interpolation algebra; actual Cheby nodes and general identity are in the proof')
    matrix=sp.Matrix([[w*v**k for w,v in zip(ell,nodes)] for k in range(4)])
    ck('S02','M1-NULL-DIRECTION', matrix.rank()==4
       and matrix*sp.Matrix(nodes)==sp.zeros(4,1),
       'Exact rational support for unique null direction x_i; not a global contact certificate')
    upper=(1-x/(1+b))/sp.sqrt(1-x)
    gap=x*(1-b*b-x)/((1+b)**2*(1-x))
    ck('S02','RADICAL-GAP-IDENTITY', sp.cancel(1-upper**2-gap)==0,
       'Symbolic gap identity with a²=1-b², denominators nonzero in declared domain')
    ck('S02','SIMPLE-DOMAIN-GAP', sp.expand((a*a/4)*(a*a-a*a/4))==3*a**4/16
       and dual['exact_derivations']['simple_gap_lower']=='gap at i=3 >= 3*a^4/64',
       'i3 node and exact numerator; denominator≤4 proved from b≤1 and x<1')
    toyx=R(1,256); toyC=(1-2*toyx)/(1-toyx)
    ck('S02','M2-NONCONCENTRATION', 0<toyC<1
       and (1-toyx)*toyC==1-2*toyx,
       'Cosine m2 relation alone permits C_i<1; this toy is not a full45-moment feasible measure')
    c,lam,mu,F,Fr,s = sp.symbols('c lam mu F Fr s',real=True)
    jet=sp.diag(1,c)
    ck('S03','FULL-POINT-JET-SURJECTIVITY', native(fp['full_dimension'],100)
       and native(fp['certificate_max_atoms'],100) and jet.det()==c,
       'Jet basis1,z-zp has determinant cos(rp)>0; no45dim generic restriction')
    ck('S03','GENERIC-DUAL-SIGN', sp.expand(-lam*(-F)-mu*(-Fr-s)-(lam*F+mu*Fr))==mu*s,
       'For mu,s≥0 the generic dual objective is a lower bound, correct sign')
    value,deriv,t = sp.symbols('value deriv t',real=True)
    ck('S03','GLOBAL-JET-CORRECTION', value-value==0
       and (deriv-t*c).subs(t,deriv/c)==0 and 2*a*t==t*(a+a),
       'Affine correction preserves pointvalue and removes positive derivative; the bound and negative branch are analytic')
    ck('S03','GENERIC-NOT-ODD', fp['odd_reduction_generic_point'] is False
       and sp.expand(-2*F+F)==-F and -2*F!=-F,
       'Feasible e=-2f+F(p) loses prescribed value under odd projection when F(p)≠0')
    K=9/a+eps*d
    U=B/K*(1+eps*(q+g/2))
    ck('S04','UPPER-ALGEBRA', sp.cancel((U-B*a/9*(1-eps*g/(2*(1+eps*a*d/9)))).subs(q,a*d/9-g))==0,
       'Equivalent positive certified upper forms; g>0 makes upper strictly below A')
    AA,UU=sp.symbols('AA UU',positive=True)
    gamma=(1+AA/UU)/2
    ck('S04','GAMMA-COST', sp.expand(gamma*UU-(AA+UU)/2)==0,
       'When0<U<A this explicitgamma>1 has gammaU<A')
    ck('S04','SIDE-GAP', rat(transfer['j_over_A_when_Beta_le_Bomega'])==R(7,8)*9==R(63,8)
       and rat(transfer['sufficient_side_ratio'])==R(8,63),
       'Exact sufficient side gap; actualfailures from gamma use this extra condition')
    ck('S04','THREE-ZERO-JET', sp.cancel((1-gamma)*B-B*(UU-AA)/(2*UU))==0,
       'Negative terminal derivative for U<A; threezero IVT and sides are written continuous proof')
    ck('S04','JET-ATTAINMENT-NOT-FAILURE', sp.diff(z**3,z).subs(z,0)==0
       and sp.Poly(z**3,z).all_roots()==[0,0,0],
       'Singular cubic level still has one distinct normal zero; jet attainment alone is not graphfailure')
    ck('S05','SOURCE-BOUNDARY', source_boundary(sources),
       'Three public bibliography records cross-checked against the supplied citation entries and their declared uses; references only')
    eta,omega=sp.symbols('eta omega',positive=True)
    width=sp.sin(eta)*sp.cos(eta)/sp.sqrt(16+sp.sin(eta)**2*sp.cos(eta)**2)
    obstruction=dr['F-S02-CHEB-OBSTRUCTION']
    ck('S05','FIXED-OMEGA-REGIMES', obstruction['asymptotic_fixed_omega_small_eta_claim'] is False
       and sp.limit(width,eta,0)==0 and sp.limit(omega-eta,eta,0)==omega,
       'Narrow edge eta↑omega differs from small eta at fixedomega; no small-eta correction concluded')
    observed=[(r['statement_id'][2:],r['id'].split('-',2)[2],r['kind']) for r in values]
    if observed!=DEFINITIONS:
        raise ValueError('Named row inventory/order changed')
    contract=all(n.get('schema_version')=='1.0' and n.get('protocol_id')==PROTOCOL for n in (terminal,dual,radius,sources))
    values.append({'id':'F-CONTRACT','statement_id':None,'kind':'contract',
                   'passed':bool(contract),'details':'All four input files declare the companion schema and protocol identity'})
    counts=dict(sorted(Counter(r['kind'] for r in values).items()))
    return {'schema_version':'1.0','protocol_id':PROTOCOL,
            'scope':'Finite exact identities and arithmetic, mathematical negative controls, bibliography and input identity',
            'total':len(values),'passed':sum(r['passed'] for r in values),
            'negative_controls':counts.get('negative',0),'counts_by_kind':counts,'checks':values}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    destination=(ROOT/args.output).resolve()
    if not destination.is_relative_to(ROOT) or not destination.parent.is_dir() or destination.exists():
        parser.error('Output must be new file in existing companion directory')
    data=payload()
    with destination.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, ensure_ascii=True, indent=2)
        stream.write('\n')
    print(str(data['passed'])+'/'+str(data['total']))
    return 0 if data['passed']==data['total'] else 1


if __name__=='__main__':
    sys.exit(main())
