"""La copie JavaScript de la loi (app/static/js/loi.js) rend les mêmes nombres que noyau/loi.py."""
import json
import os
import shutil
import subprocess

import numpy as np
import pytest

from noyau.loi import Loi

JS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "app", "static", "js", "loi.js"))
CAS = [dict(pj=3e13, ft=1.3e6, coh=16.4e6, phi=23.0, GfI=3.8, GfII=84.0),
       dict(pj=3e14, ft=34e6, coh=13.6e6, phi=13.4, GfI=70.0, GfII=700.0, c=3.0, montee="parabolic"),
       dict(pj=3e13, ft=1.3e6, coh=16.4e6, phi=23.0, GfI=3.8, GfII=84.0, adoucissement="linear")]


@pytest.mark.skipif(shutil.which("node") is None, reason="node absent")
@pytest.mark.parametrize("p", CAS)
def test_js_egal_python(p):
    L = Loi(p["pj"], p["ft"], p["coh"], p["phi"], p["GfI"], p["GfII"], c=p.get("c", 6.0),
            adoucissement=p.get("adoucissement", "yan"), montee=p.get("montee", "linear"))
    dn = list(np.linspace(0, 1.3 * (L.dnE + L.ot), 300))
    ds = list(np.linspace(0, 2.0 * L.st, 600))
    script = ("import { Loi } from %s; const L = new Loi(%s); "
              "console.log(JSON.stringify({I: L.I, s1: %s.map((x) => L.sigma(x)), t2: L.mode2(%s, -2e7)}));"
              % (json.dumps("file:///" + JS.replace("\\", "/")), json.dumps(p), json.dumps(dn), json.dumps(ds)))
    r = json.loads(subprocess.run(["node", "--input-type=module", "-e", script], capture_output=True, text=True, check=True).stdout)
    assert r["I"] == pytest.approx(L.I, rel=1e-14)
    assert np.allclose(r["s1"], [L.sigma(x) for x in dn], rtol=1e-12, atol=0)
    assert np.allclose(r["t2"], L.mode2(ds, -2e7), rtol=1e-12, atol=0)
