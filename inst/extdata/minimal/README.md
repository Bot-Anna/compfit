# Minimal fixture scenario

A tiny **2-compartment linear model** used by the test suite so the integration
tests can run without any private data.

```
dX1/dt = -k*X1            (X1: fitted initial state)
dX2/dt =  k*X1 - m*X2     (X2: fixed initial 0; k, m: fitted params)
```

- `modelParams.csv` — the model sheet (`_Level1` compartment indices, `Others`
  time grid 2015–2020 / partition 4, `States`, `Parameters`, then the coefficient
  columns interleaved per compartment: `Linear1`, `Quadratic1`, `Linear2`,
  `Quadratic2`. `Linear<j>` is the first-order column for compartment j;
  `Quadratic*` = zero — no infection term).
- `dataCombined.csv` — two observed streams (`X1`, `X2`) over 2015–2020.
- `dataDummy.csv` — one display-only series (`X1`, the susceptibles).

`test-unit-modelchain.R` builds this through the pure-R chain (no Julia). The
integration tests default to this folder when `FITCM_SCENARIO_DIR` is unset, so
they run end-to-end on any machine with a working Julia bridge.
