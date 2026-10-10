// UAV Kit calculation engine (SI units throughout)
const G0 = 9.80665, R = 287.05, GAMMA = 1.4;

function isa(h) {
  // International Standard Atmosphere, valid 0 to 20 km
  let T, p;
  if (h <= 11000) {
    T = 288.15 - 0.0065 * h;
    p = 101325 * Math.pow(T / 288.15, 5.25588);
  } else {
    T = 216.65;
    p = 22632.06 * Math.exp(-G0 / (R * T) * (h - 11000));
  }
  const rho = p / (R * T);
  const a = Math.sqrt(GAMMA * R * T);
  const mu = 1.458e-6 * Math.pow(T, 1.5) / (T + 110.4); // Sutherland
  return { T, p, rho, a, mu, sigma: rho / 1.225 };
}

function wingSizing({ mass, vStall, clMax, ar, vCruise, h }) {
  const { rho, mu } = isa(h);
  const W = mass * G0;
  const S = 2 * W / (rho * vStall * vStall * clMax);
  const b = Math.sqrt(ar * S);
  const c = S / b;
  const clCruise = 2 * W / (rho * vCruise * vCruise * S);
  const Re = rho * vCruise * c / mu;
  return {
    W, S, b, c, clCruise, Re,
    wingLoading: W / S,                       // N/m2
    wingLoadingG: (mass * 1000) / (S * 100),  // g/dm2 (hobby unit)
  };
}

function fixedWingPerf({ mass, S, ar, cd0, e, h, capAh, volts, dod, eta, clMax }) {
  const { rho } = isa(h);
  const W = mass * G0;
  const vStall = Math.sqrt(2 * W / (rho * S * clMax));
  const k = 1 / (Math.PI * e * ar);
  const LDmax = 1 / (2 * Math.sqrt(cd0 * k));
  const clMd = Math.sqrt(cd0 / k);
  const vMd = Math.sqrt(2 * W / (rho * S * clMd));      // best range (electric)
  const vMp = vMd / Math.pow(3, 0.25);                  // best endurance
  const pReq = v => 0.5 * rho * v ** 3 * S * cd0 + 2 * k * W * W / (rho * v * S);
  const energyWh = capAh * volts * dod;
  const endur = v => energyWh / (pReq(v) / eta);        // hours
  const curve = [];
  for (let v = 0.6 * vMp; v <= 2.2 * vMd; v += (2.2 * vMd - 0.6 * vMp) / 60)
    curve.push([v, pReq(v)]);
  return {
    k, LDmax, clMd, vMd, vMp, vStall,
    pMd: pReq(vMd), pMp: pReq(vMp),
    energyWh,
    endurMp: endur(vMp) * 60,               // min
    rangeMd: vMd * endur(vMd) * 3.6,        // km
    curve,
  };
}

function multirotor({ mass, n, propIn, fm, etaElec, h, capAh, volts, dod, maxThrustG }) {
  const { rho } = isa(h);
  const W = mass * G0;
  const T = W / n;
  const D = propIn * 0.0254;
  const A = Math.PI * D * D / 4;
  const pIdeal = Math.pow(T, 1.5) / Math.sqrt(2 * rho * A);
  const pShaft = pIdeal / fm;
  const pElec = n * pShaft / etaElec;
  const energyWh = capAh * volts * dod;
  const tw = (maxThrustG * n / 1000) / mass;
  return {
    T, A, diskLoading: T / A, vInduced: Math.sqrt(T / (2 * rho * A)),
    pIdeal, pShaft, pElec,
    current: pElec / volts,
    endurance: energyWh / pElec * 60,     // min
    tw, energyWh,
  };
}

function battery({ cells, chem, capmAh, cRate, avgA, massG }) {
  const vCell = chem === 'liion' ? 3.6 : 3.7;
  const V = cells * vCell;
  const Ah = capmAh / 1000;
  const Wh = V * Ah;
  return {
    V, Wh,
    maxA: cRate * Ah,
    time80: avgA > 0 ? 0.8 * Ah / avgA * 60 : NaN,
    specE: massG > 0 ? Wh / (massG / 1000) : NaN,
    vFull: cells * 4.2, vEmpty: cells * (chem === 'liion' ? 3.0 : 3.5),
  };
}

function cg(items, macLE, mac) {
  let m = 0, mx = 0;
  for (const it of items) { m += it.mass; mx += it.mass * it.x; }
  const x = m > 0 ? mx / m : NaN;
  return { total: m, x, pctMac: mac > 0 ? (x - macLE) / mac * 100 : NaN };
}

/* ---------------- shared helpers ---------------- */
const D2R = Math.PI / 180;
function polar({ mass, S, ar, cd0, e, h }) {
  const { rho } = isa(h);
  const W = mass * G0, k = 1 / (Math.PI * e * ar);
  const pReq = (v, n = 1) => 0.5 * rho * v ** 3 * S * cd0 + 2 * k * n * n * W * W / (rho * v * S);
  return { rho, W, k, pReq };
}
const helmbold = A => 2 * Math.PI * A / (2 + Math.sqrt(A * A + 4)); // 3D lift slope per rad

/* ---------------- NACA 4-digit airfoil + thin airfoil theory ---------------- */
function naca4({ code, alpha, ar, e }) {
  const s = String(code).padStart(4, '0').slice(0, 4);
  const m = +s[0] / 100, p = +s[1] / 10, t = +s.slice(2) / 100;
  const cam = x => (m === 0 || p === 0) ? [0, 0] : x < p
    ? [m / (p * p) * (2 * p * x - x * x), 2 * m / (p * p) * (p - x)]
    : [m / ((1 - p) ** 2) * (1 - 2 * p + 2 * p * x - x * x), 2 * m / ((1 - p) ** 2) * (p - x)];
  const yt = x => 5 * t * (0.2969 * Math.sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1036 * x ** 4);
  const up = [], lo = [], cl = [];
  let tMax = 0, xT = 0;
  for (let i = 0; i <= 100; i++) {
    const x = 0.5 * (1 - Math.cos(Math.PI * i / 100));
    const [yc, dy] = cam(x), th = Math.atan(dy), y = yt(x);
    up.push([x - y * Math.sin(th), yc + y * Math.cos(th)]);
    lo.push([x + y * Math.sin(th), yc - y * Math.cos(th)]);
    cl.push([x, yc]);
    if (2 * y > tMax) { tMax = 2 * y; xT = x; }
  }
  // thin airfoil integrals, x = (1 - cos th)/2
  let I0 = 0, A1 = 0, A2 = 0; const N = 4000;
  for (let i = 0; i < N; i++) {
    const th = (i + 0.5) * Math.PI / N, x = 0.5 * (1 - Math.cos(th)), dz = cam(x)[1], d = Math.PI / N;
    I0 += dz * (Math.cos(th) - 1) * d; A1 += dz * Math.cos(th) * d; A2 += dz * Math.cos(2 * th) * d;
  }
  const aL0 = -I0 / Math.PI; A1 *= 2 / Math.PI; A2 *= 2 / Math.PI;
  const a0 = 2 * Math.PI, a3 = a0 / (1 + a0 / (Math.PI * e * ar));
  const a = alpha * D2R;
  return {
    m, p, t, tMax, xT, aL0: aL0 / D2R, cm4: Math.PI / 4 * (A2 - A1),
    cl2d: a0 * (a - aL0), a3d: a3 * D2R, CL3d: a3 * (a - aL0),
    shape: up.slice().reverse().concat(lo.slice(1)), camber: cl,
  };
}

/* ---------------- Climb & glide ---------------- */
function climbGlide({ mass, S, ar, cd0, e, clMax, h, pMax, eta, hGlide }) {
  const { rho, W, k, pReq } = polar({ mass, S, ar, cd0, e, h });
  const vS = Math.sqrt(2 * W / (rho * S * clMax)), vMin = 1.2 * vS;
  const pAv = eta * pMax;
  let best = { roc: -1e9 }, bestAng = { g: -1e9 }, vMax = NaN; const curve = [];
  for (let v = vMin; v <= 120; v += 0.05) {
    const roc = (pAv - pReq(v)) / W;
    if (roc >= 0) vMax = v;
    if (roc > best.roc) best = { roc, v };
    const g = Math.asin(Math.max(-1, Math.min(1, roc / v)));
    if (g > bestAng.g) bestAng = { g, v };
    if (v < 4 * vS + 25 && Math.round(v * 20) % 10 === 0) curve.push([v, roc]);
  }
  const vMd = Math.sqrt(2 * W / (rho * S * Math.sqrt(cd0 / k)));
  const vMp = vMd / Math.pow(3, 0.25);
  const vBg = Math.max(vMd, vMin), vMs = Math.max(vMp, vMin);
  const ld = v => W * v / pReq(v);
  return {
    vS, pAv, rocMax: best.roc, vY: best.v, gammaMax: bestAng.g / D2R, vX: bestAng.v, vMax,
    LDmax: 1 / (2 * Math.sqrt(cd0 * k)), vMd, vMp, vBg, vMs, ldBg: ld(vBg),
    sinkMin: pReq(vMs) / W, glideRange: hGlide * ld(vBg) / 1000,
    glideTime: hGlide / (pReq(vMs) / W) / 60, curve,
  };
}

/* ---------------- Turn ---------------- */
function turnPerf({ mass, S, ar, cd0, e, clMax, h, V, bank }) {
  const { rho, W, pReq } = polar({ mass, S, ar, cd0, e, h });
  const phi = bank * D2R, n = 1 / Math.cos(phi);
  const vS = Math.sqrt(2 * W / (rho * S * clMax));
  return {
    n, R: V * V / (G0 * Math.tan(phi)), omega: G0 * Math.tan(phi) / V / D2R,
    t360: 2 * Math.PI * V / (G0 * Math.tan(phi)), vS, vSturn: vS * Math.sqrt(n),
    clReq: n * 2 * W / (rho * V * V * S), pLevel: pReq(V), pTurn: pReq(V, n),
  };
}

/* ---------------- V-n diagram ---------------- */
function vnDiagram({ mass, S, clMax, clMin, nPos, nNeg, vC, vdF, h }) {
  const { rho } = isa(h), W = mass * G0;
  const vS = Math.sqrt(2 * W / (rho * S * clMax));
  const vSn = Math.sqrt(2 * W / (rho * S * Math.abs(clMin)));
  const vA = vS * Math.sqrt(nPos), vG = vSn * Math.sqrt(Math.abs(nNeg)), vD = vdF * vC;
  const env = [];
  const vTop = Math.min(vA, vD);
  for (let v = 0; v <= vTop; v += vTop / 40) env.push([v, (v / vS) ** 2]);
  env.push([vTop, (vTop / vS) ** 2]);
  if (vA < vD) env.push([vA, nPos], [vD, nPos]);
  env.push([vD, 0]);
  if (vG < vC) env.push([vC, nNeg], [vG, nNeg]); else env.push([vG, nNeg]);
  for (let v = vG; v >= 0; v -= vG / 40) env.push([v, -((v / vSn) ** 2)]);
  env.push([0, 0]);
  return { vS, vSn, vA, vG, vD, vC, env };
}

/* ---------------- Takeoff & landing (numerical integration) ---------------- */
function takeoffLanding({ mass, S, ar, cd0, e, clMax, clG, h, T0, vPitch, mu, muB }) {
  const { rho } = isa(h), W = mass * G0, k = 1 / (Math.PI * e * ar);
  const vS = Math.sqrt(2 * W / (rho * S * clMax)), vLO = 1.1 * vS, vTD = 1.15 * vS;
  const cd = cd0 + k * clG * clG, dt = 0.005;
  let v = 0, s = 0, t = 0, ok = true;
  while (v < vLO) {
    const q = 0.5 * rho * v * v, T = Math.max(0, T0 * (1 - v / vPitch));
    const a = (T - q * S * cd - mu * Math.max(0, W - q * S * clG)) / mass;
    if (a <= 0 && t > 0.5) { ok = false; break; }
    v += a * dt; s += v * dt; t += dt; if (t > 300) { ok = false; break; }
  }
  let vl = vTD, sl = 0, tl = 0;
  while (vl > 0) {
    const q = 0.5 * rho * vl * vl;
    const a = (q * S * cd + muB * Math.max(0, W - q * S * clG)) / mass;
    vl -= a * dt; sl += Math.max(vl, 0) * dt; tl += dt; if (tl > 300) break;
  }
  return { ok, vS, vLO, vTD, vApp: 1.3 * vS, sTO: s, tTO: t, vReached: v, sL: sl, tL: tl, tw: T0 / W };
}

/* ---------------- Tail sizing & static margin ---------------- */
function tailSM({ S, c, b, ar, lh, lv, vh, vv, arH, arV, xcg, eta }) {
  const Sh = vh * c * S / lh, bh = Math.sqrt(arH * Sh), ch = Sh / bh;
  const Sv = vv * b * S / lv, hv = Math.sqrt(arV * Sv), cv = Sv / hv;
  const aw = helmbold(ar), at = helmbold(arH), deda = 2 * aw / (Math.PI * ar);
  const hnp = 0.25 + eta * vh * (at / aw) * (1 - deda);
  return { Sh, bh, ch, Sv, hv, cv, aw, at, deda, hnp: hnp * 100, sm: hnp * 100 - xcg };
}

/* ---------------- Spar bending ---------------- */
const MATERIALS = {
  carbon: { name: 'Carbon tube (pultruded)', sig: 600e6, E: 120e9, rho: 1600 },
  al6061: { name: 'Aluminium 6061-T6', sig: 276e6, E: 69e9, rho: 2700 },
  al7075: { name: 'Aluminium 7075-T6', sig: 503e6, E: 72e9, rho: 2810 },
  glass:  { name: 'Glass-fibre tube', sig: 300e6, E: 25e9, rho: 1900 },
  spruce: { name: 'Spruce (wood)', sig: 65e6, E: 9e9, rho: 450 },
};
function spar({ mass, n, b, sf, od, wall, mat }) {
  const M = MATERIALS[mat], W = mass * G0;
  const Do = od / 1000, Di = Math.max(0, od - 2 * wall) / 1000;
  const I = Math.PI * (Do ** 4 - Di ** 4) / 64, A = Math.PI * (Do * Do - Di * Di) / 4;
  const Mroot = n * W * b / (3 * Math.PI);          // elliptic lift, cantilever half-span
  const sigma = Mroot * (Do / 2) / I;
  const L = b / 2, w = (n * W / 2) / L;
  return {
    Mroot, sigma, allow: M.sig, ms: M.sig / (sf * sigma) - 1,
    defl: w * L ** 4 / (8 * M.E * I), deflPct: w * L ** 4 / (8 * M.E * I) / L * 100,
    sparMass: A * b * M.rho, I,
  };
}

/* ---------------- Motor & propeller ---------------- */
function motorProp({ kv, volts, amps, propIn, pitchIn, etaM, fm, h }) {
  const { rho, a } = isa(h);
  const rpmNL = kv * volts, rpm = 0.8 * rpmNL;
  const D = propIn * 0.0254, A = Math.PI * D * D / 4;
  const pE = volts * amps, pS = etaM * pE;
  const T = Math.pow(fm * pS, 2 / 3) * Math.pow(2 * rho * A, 1 / 3);
  return {
    rpmNL, rpm, pitchSpeed: rpm * pitchIn * 0.0254 / 60, pE, pS, T,
    Tg: T / G0 * 1000, gPerW: T / G0 * 1000 / pE,
    tipMach: Math.PI * D * rpm / 60 / a, pitchRatio: pitchIn / propIn,
  };
}

/* ---------------- Mass budget (fixed-wing sizing loop) ---------------- */
function massBudget({ payload, avionics, fS, fP, V, LD, eta, endur, eSpec, reserve }) {
  const tH = endur / 60;
  const fB = G0 * V * tH / (LD * eta * eSpec * (1 - reserve));
  const den = 1 - fS - fP - fB;
  const m = den > 0 ? (payload + avionics) / den : NaN;
  return { fB, den, m, mS: fS * m, mP: fP * m, mB: fB * m, Wh: fB * m * eSpec, pCruise: m * G0 * V / (LD * eta) };
}

/* ---------------- Mapping / photogrammetry ---------------- */
function survey({ sw, f, imgW, imgH, alt, fo, so, V, L, Wd }) {
  const gsd = sw * alt * 100 / (f * imgW);                // cm/px
  const fw = gsd * imgW / 100, fh = gsd * imgH / 100;     // m
  const spacing = fh * (1 - fo / 100), line = fw * (1 - so / 100);
  const lines = Math.ceil(Wd / line) + 1;
  const path = lines * L + (lines - 1) * line;
  return {
    gsd, fw, fh, spacing, interval: spacing / V, line, lines,
    photos: lines * (Math.ceil(L / spacing) + 1), path, time: path / V / 60,
  };
}

/* ---------------- Radio link budget ---------------- */
function linkBudget({ fMHz, ptx, gt, gr, sens, margin, hDrone, hGnd }) {
  const budget = ptx + gt + gr - sens - margin;
  const dLink = Math.pow(10, (budget - 32.44 - 20 * Math.log10(fMHz)) / 20); // km
  const dHor = 4.12 * (Math.sqrt(hDrone) + Math.sqrt(hGnd));
  const d = Math.min(dLink, dHor);
  const fres = 17.32 * Math.sqrt(d / (4 * fMHz / 1000));
  return { budget, dLink, dHor, d, fres, limit: dLink < dHor ? 'link' : 'horizon' };
}

/* ---------------- Wind triangle ---------------- */
function windTri({ tas, ws, wdir, track, dist }) {
  const leg = trk => {
    const b = (wdir - trk) * D2R;
    const head = ws * Math.cos(b), cross = ws * Math.sin(b);
    if (Math.abs(cross) >= tas) return null;
    const wca = Math.asin(cross / tas);
    const gs = tas * Math.cos(wca) - head;
    return { head, cross, wca: wca / D2R, hdg: ((trk + wca / D2R) % 360 + 360) % 360, gs };
  };
  const out = leg(track), back = leg((track + 180) % 360);
  const tOut = out && out.gs > 0 ? dist * 1000 / out.gs / 60 : NaN;
  const tBack = back && back.gs > 0 ? dist * 1000 / back.gs / 60 : NaN;
  return { out, back, tOut, tBack, tCalm: 2 * dist * 1000 / tas / 60 };
}

/* ---------------- Unit conversion tables (factor to SI base) ---------------- */
const UNITS = {
  Length: { m: 1, mm: 1e-3, cm: 1e-2, km: 1e3, in: 0.0254, ft: 0.3048, mi: 1609.344, nmi: 1852 },
  Speed: { 'm/s': 1, 'km/h': 1 / 3.6, kt: 0.514444, mph: 0.44704, 'ft/min': 0.00508 },
  Mass: { kg: 1, g: 1e-3, lb: 0.45359237, oz: 0.028349523 },
  Force: { N: 1, kgf: 9.80665, gf: 0.00980665, lbf: 4.448222 },
  Power: { W: 1, kW: 1e3, hp: 745.6999 },
  Energy: { J: 1, kJ: 1e3, Wh: 3600, kWh: 3.6e6 },
  Pressure: { Pa: 1, hPa: 100, kPa: 1e3, bar: 1e5, psi: 6894.757, inHg: 3386.389 },
  Area: { 'm²': 1, 'dm²': 1e-2, 'cm²': 1e-4, 'in²': 6.4516e-4, 'ft²': 0.09290304 },
  'Wing loading': { 'kg/m²': 1, 'g/dm²': 0.1, 'N/m²': 1 / 9.80665, 'oz/ft²': 0.305152, 'lb/ft²': 4.882428 },
  Angle: { deg: 1, rad: 180 / Math.PI },
};

if (typeof module !== 'undefined') module.exports = {
  isa, wingSizing, fixedWingPerf, multirotor, battery, cg, naca4, climbGlide, turnPerf,
  vnDiagram, takeoffLanding, tailSM, spar, motorProp, massBudget, survey, linkBudget, windTri, UNITS, MATERIALS };
