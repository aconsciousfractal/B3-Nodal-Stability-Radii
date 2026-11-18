"""Exact scalar and symbolic support for the written uniform support continuum proof."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import sympy as sp


def support():
    rows = []

    def record(identifier, passed, scope, kind="exact"):
        rows.append({"id": identifier, "kind": kind,
                     "status": "PASS" if bool(passed) else "FAIL", "scope": scope})

    a2 = F(625, 100020626)
    beta = F(10000, 10001**2)
    B = F(1249875000, 10001**4)
    record("K01", a2 == beta / (16 + beta), "Exact fixed normal scale.")
    record("K02", 0 < a2 < F(1, 160000), "a<1/400, fixed geometry only.")
    record("K03", 0 < B < F(1, 8000000) < F(3, 625), "Strict endpoint B ordering; angular minimum requires written unimodality proof.")
    record("K04", F(1, 400)*F(11, 10)/(1-F(1, 1600000)) < F(1, 300), "Uniform reach-ratio scalar envelope for |s|<=1,|t|<=1/10.")
    record("K05", F(49, 500)/2 == F(49, 1000), "sin(L/4)>L/8>=49/1000 uses the article's analytic sine and angular bounds.")
    record("K06", F(10, 147)**2 < F(1, 100), "Uniform discriminant margin after dividing the reach ratio by sin(L/4).")
    record("K07", 1-F(1, 100) == F(99, 100), "Strict discriminant>99b^2/100; real support follows from the written chord construction.")
    t = sp.symbols("t", real=True)
    U = 1-40*t**2+240*t**4-448*t**6+256*t**8
    T = 9*t-120*t**3+432*t**5-576*t**7+256*t**9
    N = (t**2-1)*U
    record("K08", sp.expand(sp.diff(T,t)-9*U) == 0, "Exact T9'=9U8 polynomial identity.")
    record("K09", sp.expand(N-(-1+41*t**2-280*t**4+688*t**6-704*t**8+256*t**10)) == 0, "Exact nodal polynomial; its nonvanishing on the continuum uses the U bound.")
    record("K10", F(3,5)-F(7,15625)>F(1,2), "U lower scalar envelope on t^2<=1/100.")
    record("K11", -40+F(12,5)+F(4,15625)<0, "Upper coefficient envelope implies U<1 outside zero.")
    record("K12", -80+F(48,5)+F(32,15625)<0, "U'/t<0 via the written polynomial envelope.")
    record("K13", F(1,65)/(1-F(1,65)) == F(1,64), "Tangent correction magnitude bound; a^2<1/65 is strict.")
    record("K14", 2-F(1,64)==F(127,64)>0, "Positive lambda/t lower envelope, conditional on K10,K12.")
    record("K15", F(39,5)-F(9,15625)>0, "Positive T9/t envelope outside zero.")
    record("K16", -120+F(108,25)+F(4,15625)<0, "Strict T9/t<9 envelope.")
    record("K17", 80+F(168,625)<81, "-tU'<81t^2 envelope.")
    record("K18", F(2,1)/(1-F(1,100))<3 and 3+162==165, "lambda*t<165t^2 envelope, using U>1/2 and negative tangent correction.")
    record("K19", 9+F(1485,100)==F(477,20)<24, "M<24 envelope, with central case M=9 separate.")
    record("K20", F(317,8)-F(12,5)-F(4,15625)>37, "S-U>37t^2 coefficient envelope.")
    record("K21", F(951,8)-F(108,25)-F(4,15625)>114, "9R-T9/t>114t^2 envelope.")
    record("K22", F(333,24)==F(111,8), "Conservative strict objective/mass gap; continuum proof is in the article.")
    record("K23", F(8,9)-F(1,8)==F(55,72)>0, "Critical graph sign coefficient for the critical field.")
    record("K24", F(7,8)*9==F(63,8)>1, "Side threshold j>A; graph/event lemma and uniqueness remain written proofs.")

    # A free symbolic target t is retained: these are identities, not a grid.
    # x=s-t, b=b0+b1*x, q=q2*x^2. The two-root recurrence is independent
    # of the direct jet of the original rotated global monomial.
    x,b0,b1,q2,w = sp.symbols("x b0 b1 q2 w")
    b=b0+b1*x
    q=q2*x**2
    P=[b,2*q]
    for degree in range(2,10):
        P.append(sp.expand(b*P[-1]-q*P[-2]))
    s=t+x
    nodal=(s**2-1)*(1-40*s**2+240*s**4-448*s**6+256*s**8)
    nodal=sp.expand(nodal)
    nt=sp.expand(nodal.subs(x,0))
    nd=sp.expand(sp.diff(nodal,x).subs(x,0))
    lam=-nd/nt+b1/b0
    cancellation=sp.cancel(nd*b0-nt*b1+lam*nt*b0)
    real_components=0
    for m in range(10):
        numerator=sp.expand(b**(m+1)+sum(sp.binomial(m,n)*b**(m-n)*w**n*P[n] for n in range(1,m+1)))
        for k in range(10-m):
            H=sp.Poly(sp.expand(s**k*numerator),x)
            value=sp.cancel(H.eval(0)/b0-t**k*b0**m)
            direct=(k*t**(k-1)*b0**m if k else 0)+(m*t**k*b0**(m-1)*b1 if m else 0)
            jet=sp.cancel((H.diff().eval(0)*b0-H.eval(0)*b1)/b0**2-direct)
            quotient=H.nth(10)/256 if H.degree()==10 else 0
            remainder=sp.Poly(sp.expand(H.as_expr()-quotient*nodal),x)
            components=1 if m==0 else 2
            real_components+=components
            record(f"KM-{m}-{k}", H.degree()<=10 and remainder.degree()<=9 and value==0 and jet==0 and cancellation==0,
                   f"z^{k}(u+iv)^{m}: {components} real component(s); formal degree10/b, original value/first derivative and interpolation annihilator. Nonzero physical a^k and rotated unit phase multiply both sides equally.")
    wrong=sp.cancel(nd*b0-nt*b1+(-nd/nt)*nt*b0)
    record("KN01", sp.expand(wrong+nt*b1)==0 and wrong!=0,
           "Omitting b'/b leaves -N(t)b1; nonzero for original t!=0. At t=0 the correction vanishes, handled separately.","negative")
    record("KN02", sp.cancel(2*q2/b0)!=0,
           "Replacing the double zero of q with a simple zero produces formal first-jet defect 2q2/b0. No alternative target family is executed.","negative")
    return {"schema_version":"1.0","scope":"uniform support exact scalar envelopes and symbolic identities supporting written continuum proof; no grid, numerical solver, novelty or independent acceptance.",
            "normalization":"Fixed taneta=1/100,tanomega=1/2, full100 original jet, entire closed tube norm; t symbolic.",
            "basis_coverage":{"complex_monomial_identities":55,"real_components":real_components},
            "counts":{kind:sum(row["kind"]==kind for row in rows) for kind in ("exact","negative")},
            "checks":rows,"status":"PASS" if all(row["status"]=="PASS" for row in rows) else "FAIL","scientific_acceptance":False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        parser.error("Output must be new and its parent must already exist.")
    result=support()
    with args.output.open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2)
        stream.write("\n")
    print(json.dumps({"status":result["status"],"checks":len(result["checks"]),"counts":result["counts"],"basis_coverage":result["basis_coverage"]},ensure_ascii=False))
    return 0 if result["status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())
