"""Exact support for Example M patch arithmetic and identities; proofs are separate."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import sympy as sp
from uniform_support import support as prior_support


def support():
    prior = prior_support()
    # K08 onward is used in the adapted criterion for free formal parameters
    # or the original scalar domain a²<1/65. K's fixed geometry is excluded.
    rows = [dict(row, id="M-import-" + row["id"])
            for row in prior["checks"]
            if row["id"] not in {f"K{i:02d}" for i in range(1, 8)}]

    def record(identifier, passed, scope, kind="exact"):
        rows.append({"id": identifier, "kind": kind,
                     "status": "PASS" if bool(passed) else "FAIL", "scope": scope})

    v = Q(25, 10001)
    a2 = v*v/(1+v*v)
    record("M01", a2 == Q(625, 100020626),
           "Exact normal scale for taneta=1/100; a=25/sqrt100020626.")
    record("M02", 0 < a2 < Q(1, 160000) < Q(1, 65),
           "Normal scale and inherited original scalar domain; no float.")
    record("M03", v/(1+v*v) == Q(250025, 100020626) > Q(3, 1250),
           "Rational lower envelope for delta; atan integral is proved in the article.")
    record("M04", Q(3, 1250) == Q(6, 25)*Q(1, 100),
           "eta<1/100 makes delta>6eta/25 strict.")
    record("M05", Q(1, 100)+Q(21, 5)*Q(1, 400)
           == Q(41, 2000) < Q(1, 40),
           "omega upper envelope; angular ordering and pi comparison are analytic.")
    record("M06", Q(21, 5) > 4 and Q(21, 5)-Q(3, 2) == Q(27, 10),
           "L>4delta and adapted side alpha>=27delta/10.")
    record("M07", 1+Q(3, 2)*Q(6, 25) == Q(34, 25) > Q(4, 3),
           "theta_c/eta lower envelope; strict by M03-M04 and the proof.")
    record("M08", 1-Q(1, 400)**2/6 > Q(99, 100),
           "sin(delta)/delta lower coefficient; signed Taylor proof is separate.")
    record("M09", 11*Q(99, 100)/(10+Q(1, 10))
           == Q(1089, 1010) > Q(21, 20),
           "W/delta lower beats sin(L/4)/delta upper: candidate beyond L frontier.")
    record("M10", Q(27, 8000) < Q(1, 100)
           and 1-Q(1, 100)**2/2 > Q(9, 10),
           "Cosine envelope used for the entire adapted support, not sampled angles.")
    record("M11", Q(9, 10)*Q(27, 20) == Q(243, 200)
           and Q(11, 10)/Q(243, 200) == Q(220, 243) < 1,
           "Strict adapted support reach ratio; physical sine hypothesis proved in the article.")
    record("M12", 1-Q(220, 243)**2 == Q(10649, 59049) > 0,
           "Uniform discriminant coefficient on near-terminal patch.")

    a, s, t = sp.symbols("a s t", real=True)
    h = sp.symbols("h", positive=True)
    c2 = 1-a*a*t*t
    b2 = (1-a*a*t*s)**2/c2
    q = a*a*(s-t)**2/(c2*h)
    record("M13", sp.cancel(b2-(1-a*a*s*s)-a*a*(s-t)**2/c2) == 0,
           "Exact tangent affine defect for free theta-dependent h.")
    record("M14", sp.cancel(4*q/b2-4*a*a*(s-t)**2/(h*(1-a*a*t*s)**2)) == 0,
           "Exact discriminant ratio with h=4sin²(alpha/2); support is analytic.")
    reach = a*(1+t)/(1+a*a*t)
    record("M15", sp.cancel(sp.diff(reach, t)
           -a*(1-a*a)/(1+a*a*t)**2) == 0
           and sp.cancel(reach.subs(t, Q(1, 10))-11*a/(10+a*a)) == 0,
           "Uniform endpoint reach threshold; maximum over s proved in the article.")
    record("M16", 1-Q(1, 20)**2/6 == Q(2399, 2400)
           and 1-Q(1, 20)**2/2 == Q(799, 800),
           "Target sine/cosine coefficient bounds on theta<=omega<1/40.")
    record("M17", Q(2399, 2400)**3 > Q(799, 800),
           "Strict Bernoulli scalar factor in angular target ratio.")
    record("M18", Q(799, 800)**2 > Q(399, 400),
           "Product envelope for target; no norm or coefficient restriction.")
    record("M19", Q(4, 3)**3*Q(399, 400) == Q(532, 225) > Q(9, 4),
           "Angular target margin on the entire far patch follows from written bounds.")

    u = sp.symbols("u", real=True)
    barrier = 81*u*(1-u/8)**2+sp.Rational(16, 65)*(1-u)*(1-3*u/8)**2
    derivative = (81*(1-u/8)*(1-3*u/8)
                  +sp.Rational(16, 65)*(-(1-3*u/8)**2
                    -sp.Rational(3, 4)*(1-u)*(1-3*u/8)))
    record("M20", sp.expand(sp.diff(barrier, u)-derivative) == 0
           and Q(2835, 64)-Q(28, 65) > 0,
           "H barrier derivative identity and envelope; continuum monotonicity in the article.")
    record("M21", barrier.subs(u, 0) == Q(16, 65)
           and barrier.subs(u, 1) == 81*Q(7, 8)**2,
           "Global coarse H lower including value-only normal sides.")
    record("M22", barrier.subs(u, Q(1, 100))
           > Q(67127049, 65000000) > 1,
           "Original H outer-band exclusion, no new fiber model.")
    record("M23", Q(9, 4)**2 > Q(65, 16) and 81 > 65,
           "B margin compensates global coarse H factor: 9/sqrt65>1.")
    b_eta = (2*Q(1, 100)/(1+Q(1, 100)**2))**3 \
        * (1-Q(1, 100)**2)/(1+Q(1, 100)**2)/64
    record("M24", b_eta == Q(1249875000, 10001**4) > 0,
           "Exact endpoint target coefficient; rho=aB_eta/9 only after written coverage.")
    record("M25", Q(7, 8)*9 == Q(63, 8) > 7,
           "Side margin and gamma in(1,7) true-failure family: IVT/graph proof separate.")
    record("MN01", 1-(Q(1089, 1010)/Q(21, 20))**2 < 0,
           "Uniform prescribed L support has negative extreme discriminant in candidate; no radius counterexample.",
           "negative")
    record("MN02", barrier.subs(u, 0) < 1,
           "Coarse H alone at theta=eta cannot certify target A_eta; target margin is required where this patch is used.",
           "negative")
    return {"schema_version": "1.0",
            "scope": "Example M exact scalar envelopes and formal moment identities; no grid, solver, continuum or novelty certification.",
            "normalization": "Original real full100 and continuous whole closed-tube sup. Near-angle support and far-angle margin are proved in the article.",
            "basis_coverage": prior["basis_coverage"],
            "counts": {kind: sum(row["kind"] == kind for row in rows)
                       for kind in ("exact", "negative")},
            "checks": rows,
            "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
            "scientific_acceptance": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        parser.error("Output must be new and its parent must already exist.")
    result = support()
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "checks": len(result["checks"]),
                      "counts": result["counts"], "basis_coverage": result["basis_coverage"]}))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
