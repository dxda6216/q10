#!/usr/bin/env python3
"""
================================================================================
  Q10 calculator for circadian period - desktop (Python/Tkinter) version
================================================================================
  Fits   period(T) = tau_b * Q10 ** ((T_b - T) / 10)   to temperature/period data
  and shows the fit as Period-vs-Temperature, Frequency-vs-Temperature and
  Arrhenius plots, with an optional Arrhenius regression line.

  Requirements:  Python 3.9+, numpy, scipy, matplotlib   (tkinter ships with Python)
  Run:           python q10_calculator.py

  A browser version of this calculator (single HTML file, no installation) is
  available from  https://github.com/dxda6216/q10/blob/main/index.html

  ==============================================================================
  1. Theory
  ==============================================================================
  [Q10 (temperature coefficient)]
    The factor by which a rate R changes for a 10 °C rise in temperature.
    From rates R1, R2 measured at temperatures T1, T2:
        Q10 = (R2 / R1) ^ (10 / (T2 - T1))
    Typical enzyme-catalysed reactions have Q10 of about 2 to 3.

  [Application to circadian rhythms]
    The "rate" of a circadian clock is taken to be its frequency, f = 1/period.
    Circadian clocks are temperature-compensated: their period depends only
    weakly on temperature, so Q10 is typically close to 1 (roughly 0.8 to 1.2).
      Q10 > 1 : period shortens as temperature rises (the clock runs faster)
      Q10 = 1 : perfect temperature compensation
      Q10 < 1 : period lengthens as temperature rises (overcompensation)

  [Model that is fitted]
    With tau_b the period at the base temperature T_b:
        tau(T) = tau_b * Q10 ^ ((T_b - T) / 10)                         ... (1)
    Frequency in cycles per day (day^-1):
        f(T) = 24 / tau(T)
    Taking logarithms:
        ln f(T) = ln(24 / tau_b) + ((T - T_b) / 10) * ln Q10            ... (2)
    So ln f is a straight line in T with slope ln(Q10) / 10.

  [Meaning of the base temperature T_b]
    T_b is only a reparameterisation: it decides at which temperature the
    fitted period tau_b is reported. Q10, its standard error, R^2 and the
    fitted curve itself are identical for every choice of T_b; only tau_b and
    its error change. Choosing T_b far outside the data (e.g. absolute zero,
    -273.15 °C) makes tau_b a large extrapolation with a large error.

  [Arrhenius plot]
    The Arrhenius equation k = A * exp(-Ea / (R * T_K)) gives a straight line
    when ln k is plotted against 1/T_K (T_K = absolute temperature), with
    slope -Ea/R (Ea: activation energy, R = 8.314 J mol^-1 K^-1).
    This program plots ln f against 1000/T_K (K^-1). Model (2) is linear in T,
    not in 1/T, so the fitted Q10 curve is slightly bent on the Arrhenius plot.
    The two descriptions are approximately related by
        Ea ~= R * T1 * T2 * ln(Q10) / 10     (T1, T2 in kelvin)

  [Arrhenius regression line]
    Optionally (checkbox, on by default) the Arrhenius plot also shows the
    ordinary least-squares straight line of ln f against 1000/T_K, fitted to
    the same data used for Q10:
        ln f = intercept + slope * (1000 / T_K),    Ea = -slope * 1000 * R
    It is drawn over the same temperature span as the Q10 curve, and its Ea,
    slope, intercept, standard errors and R^2 appear in the Results. Because
    the two models differ in form (linear in T vs. linear in 1/T) and in the
    quantity whose residuals are minimised (period vs. ln f), the two lines
    are close but not identical; they diverge slightly outside the data range.

  ==============================================================================
  2. Calculation method
  ==============================================================================
    (a) Non-linear least squares (scipy.optimize.curve_fit) minimising the sum
        of squared residuals of the period itself (not its logarithm):
            SSE = sum_i (tau_i - tau_hat(T_i))^2
    (b) Starting values: linear regression of ln(tau) on T, converted to
        tau_b and Q10 through relation (2) (falls back to (24, 1)).
    (c) Standard errors: curve_fit's covariance matrix (absolute_sigma=False),
        i.e. Cov = s^2 * (J^T J)^-1 with s^2 = SSE / (n - 2), where J is the
        Jacobian with respect to tau_b and Q10. The "±" values are the square
        roots of its diagonal.
    (d) Coefficient of determination: R^2 = 1 - SSE / SST,
        SST = sum_i (tau_i - mean(tau))^2
    (e) Fitted curve: the full temperature range of the data, extended on each
        side by 33.3 % of its width (end points truncated to integers), drawn
        with 200 points.
    (f) Arrhenius regression: x_i = 1000 / (T_i + 273.15), y_i = ln(24 / tau_i),
        simple linear regression y = intercept + slope * x on the fitted data.
            slope = Sxy / Sxx,  intercept = mean(y) - slope * mean(x)
            s^2 = SSE / (n - 2),  SE(slope) = sqrt(s^2 / Sxx),
            SE(intercept) = sqrt(s^2 * (1/n + mean(x)^2 / Sxx))
            Ea = -slope * 1000 * R  (J/mol),  SE(Ea) = SE(slope) * 1000 * R
        with R = 8.314462618 J mol^-1 K^-1.

  ==============================================================================
  3. How to use
  ==============================================================================
    1) On start-up, the window shows a demo calculation with the example data.
    2) Click "Clear" to empty the form.
    3) Enter one "temperature(°C) period(h)" pair per line in the data box
       (separate the values with commas, spaces, tabs or semicolons),
       or use "Load CSV..." to read a file with temperature in column 1 and
       period in column 2.
    4) Optionally tick the temperature-range box and set the range used for
       the fit. Points outside the range are drawn in light green and are not
       used in the calculation. After a calculation, the range can also be
       chosen by dragging across any plot; the fit is recalculated on release.
       Double-click a plot or click "Use all data" to return to all data.
       The axis scales of the plots always match the all-data fit.
    5) Choose the base temperature. Minimum / Maximum / Average / Median are
       computed from the data used in the fit.
    6) Click "Calculate" (or press Ctrl+Enter).
    7) Switch between the Results / Period / Frequency / Arrhenius tabs.
    8) "Save figure..." saves the three plots stacked in one image (PNG, PDF
       or SVG); "Save results..." saves the results text.
    "How it works" opens a window explaining the theory, method and usage.

  ==============================================================================
  4. Privacy, network access and offline use
  ==============================================================================
    - All calculations run locally on your computer.
    - No data is sent anywhere: the program makes no network connections.
    - Nothing is stored: the data you enter and the results exist only while
      the program is running. Files are written only when you use "Save
      figure..." or "Save results...", to the place you choose.
      (matplotlib may keep its own font cache in its configuration folder; it
      contains no data you enter.)
    - Works fully offline. The only way the program touches the internet is
      if you click the GitHub link in the start-up notice or the help window,
      which opens it in your web browser.

  ==============================================================================
  5. File layout
  ==============================================================================
    Constants          : example data, base-temperature choices, help text
    Calculation core   : parsing and fitting (no GUI dependencies)
    Report text        : builds the Results text
    Plotting           : draws each plot onto a matplotlib Axes
    Help window        : renders HELP_CONTENT in a Tk Text widget
    GUI (Q10App)       : layout, start-up notice, buttons, saving
================================================================================
"""

import csv
import math
import re
import statistics
import tkinter as tk
import tkinter.font as tkfont
import webbrowser
from tkinter import filedialog, messagebox, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
from scipy.optimize import curve_fit

# ======================================================================
# Constants
# ======================================================================

DEG = "\u00B0"
R_GAS = 8.314462618       # gas constant, J mol^-1 K^-1
HTML_URL = "https://github.com/dxda6216/q10/blob/main/index.html"  # browser version

DATA_ROWS = 15            # visible lines in the data input box
TITLE_SIZE = 18           # plot-title font size (single-plot tabs)
TITLE_SIZE_COMBINED = 20  # plot-title font size (saved 3-plot figure)

# Demo data, loaded and calculated at start-up. Temperature in °C, period in h.
EXAMPLE_TEMPS = "25.1, 26.9, 27.3, 29.2, 30, 31.7, 32.2, 33.5, 34.1, 35.3, 35.8, 36.5, 37.1, 38.5, 38.9, 39.3, 40.3, 41.1"
EXAMPLE_PERIODS = "33.6, 32.8, 30.4, 29.3, 28.0, 27.9, 26.1, 25.5, 24.9, 24.8, 24.5, 24.4, 24.1, 23.9, 23.6, 21.9, 21.5, 20.1"

EXAMPLE_DATA = "\n".join(
    f"{t.strip()}, {p.strip()}"
    for t, p in zip(EXAMPLE_TEMPS.split(","), EXAMPLE_PERIODS.split(","))
)

# Base-temperature options. Those in FIXED_BASE have a fixed value; Minimum /
# Maximum / Average / Median are computed from the data used in the fit.
CUSTOM_BASE = "Custom value (enter below)"
BASE_CHOICES = [
    f"-273.15{DEG}C (absolute zero)",
    f"0{DEG}C",
    f"4{DEG}C",
    f"25{DEG}C",
    f"30{DEG}C",
    f"37{DEG}C",
    f"100{DEG}C",
    "Minimum",
    "Maximum",
    "Average",
    "Median",
    CUSTOM_BASE,
]
FIXED_BASE = {
    BASE_CHOICES[0]: -273.15,
    BASE_CHOICES[1]: 0.0,
    BASE_CHOICES[2]: 4.0,
    BASE_CHOICES[3]: 25.0,
    BASE_CHOICES[4]: 30.0,
    BASE_CHOICES[5]: 37.0,
    BASE_CHOICES[6]: 100.0,
}

# Texts of the start-up notice (hidden by Clear, editing the data, or loading a CSV).
NOTICE_DEMO = ("You're looking at a demo calculation using example data. To analyse your own "
               "data, click Clear and enter your temperature and period pairs in the form.")
NOTICE_PRIVACY = ("All calculations run on your computer. The data you enter and the results are "
                  "not sent anywhere and are not stored; they are gone when you close the program, "
                  "unless you save them yourself with the Save buttons. The program does not connect "
                  "to the internet, so it works fully offline. A browser version (a single HTML file) "
                  "can be downloaded from:")


# ======================================================================
# Calculation core (no GUI dependencies)
# Flow: parse_pairs() -> analyze() -> fit_q10() and arrhenius_fit()
# ======================================================================

def parse_pairs(text):
    """Parse one 'temperature period' pair per line.

    - separators: comma / space / tab / semicolon (may be mixed)
    - blank lines and lines starting with '#' are ignored
    - one leading non-numeric line is skipped as a header
    - a line with other than 2 values, or an unreadable line, raises
      ValueError with its line number
    Returns (temperatures, periods) as numpy arrays.
    """
    temps, periods = [], []
    header_skipped = False
    for n, line in enumerate(text.splitlines(), start=1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        tokens = [t for t in re.split(r"[,;\s]+", s) if t]
        try:
            vals = [float(t) for t in tokens]
        except ValueError:
            if not temps and not header_skipped:
                header_skipped = True  # header line such as "Temp, Period"
                continue
            raise ValueError(f"Line {n}: cannot read numbers: {s!r}") from None
        if len(vals) != 2:
            raise ValueError(f"Line {n}: expected 2 values (temperature, period), "
                             f"found {len(vals)}: {s!r}")
        temps.append(vals[0])
        periods.append(vals[1])
    return np.array(temps, dtype=float), np.array(periods, dtype=float)


def resolve_base_temperature(selection, custom_value, xa):
    """Turn the selected base temperature T_b into a number.

    Minimum/Maximum/Average/Median are taken from xa, the data USED in the fit
    (only the in-range points when a temperature range is set).
    """
    if selection in FIXED_BASE:
        return FIXED_BASE[selection]
    if selection == "Minimum":
        return float(np.min(xa))
    if selection == "Maximum":
        return float(np.max(xa))
    if selection == "Average":
        return float(statistics.mean(xa))
    if selection == "Median":
        return float(statistics.median(xa))
    return float(custom_value)


def make_q10_model(base_t):
    """Return model (1), tau(T) = tau_b * Q10^((T_b - T)/10), for a given T_b."""
    def q10_model(temperature, tau_bt, q10):
        return tau_bt * q10 ** ((base_t - temperature) * 0.1)
    return q10_model


def initial_guess(xa, ya, base_t):
    """Starting values from a log-linear regression (more robust than a fixed (24, 1)).

    Regress ln(tau) = c + s*T, then use relation (2):
        s = -ln(Q10)/10,   c = ln(tau_b) + (T_b/10)*ln(Q10)
    giving ln(Q10) = -10*s and ln(tau_b) = c - (T_b/10)*ln(Q10).
    """
    try:
        slope, intercept = np.polyfit(xa, np.log(ya), 1)
        ln_q10 = -10.0 * slope
        ln_tau = intercept - (base_t / 10.0) * ln_q10
        p0 = (float(np.exp(ln_tau)), float(np.exp(ln_q10)))
        if all(np.isfinite(p0)) and p0[0] > 0 and p0[1] > 0:
            return p0
    except Exception:
        pass
    return (24.0, 1.0)


def fit_q10(xa, ya, base_t):
    """Non-linear least-squares fit of model (1) to the period data.

    curve_fit minimises SSE = sum (y_i - tau_hat_i)^2 and returns the covariance
    scaled by s^2 = SSE/(n - 2) (absolute_sigma=False).
    Returns (popt=[tau_b, Q10], pcov, R^2, model).
    """
    model = make_q10_model(base_t)
    p0 = initial_guess(xa, ya, base_t)
    popt, pcov = curve_fit(model, xa, ya, p0=p0, maxfev=20000)
    residuals = ya - model(xa, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ya - np.mean(ya)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return popt, pcov, r2, model


def arrhenius_fit(temps, periods):
    """Classical Arrhenius regression: least-squares line of ln(frequency) on 1000/T_K.

        x_i = 1000 / (T_i + 273.15),  y_i = ln(24 / tau_i),  y = intercept + slope * x
    Standard errors use s^2 = SSE / (n - 2):
        SE(slope) = sqrt(s^2 / Sxx),  SE(intercept) = sqrt(s^2 * (1/n + mean(x)^2 / Sxx))
    Activation energy: Ea = -slope * 1000 * R (J/mol), SE(Ea) = SE(slope) * 1000 * R.
    Called with the data used for the Q10 fit, so both fits use the same points.
    """
    t, p = np.asarray(temps, float), np.asarray(periods, float)
    x, y = 1000.0 / (t + 273.15), np.log(24.0 / p)
    n, mx, my = len(x), x.mean(), y.mean()
    sxx = np.sum((x - mx) ** 2)
    syy = np.sum((y - my) ** 2)
    slope = np.sum((x - mx) * (y - my)) / sxx
    intercept = my - slope * mx
    sse = np.sum((y - intercept - slope * x) ** 2)
    s2 = sse / (n - 2)
    se_slope = math.sqrt(s2 / sxx)
    se_int = math.sqrt(s2 * (1.0 / n + mx * mx / sxx))
    return dict(
        slope=slope, intercept=intercept, se_slope=se_slope, se_int=se_int,
        r2=1 - sse / syy if syy > 0 else float("nan"),
        Ea=-slope * 1000 * R_GAS, se_Ea=se_slope * 1000 * R_GAS,  # J/mol
    )


def arrhenius_lnf(arr, temps_c):
    """ln(frequency) predicted by the Arrhenius regression for temperatures in °C."""
    return arr["intercept"] + arr["slope"] * 1000.0 / (np.asarray(temps_c, float) + 273.15)


def analyze(temps, periods, use_range, low, high, base_sel, base_custom):
    """Run the whole calculation and return a dict with everything needed to plot/report.

      temps, periods : all data
      use_range      : restrict the fit to low <= T <= high (°C, inclusive)
      base_sel       : selected base temperature; base_custom: custom value
    Steps:
      1) validate input (equal counts, at least 3 points, periods > 0)
      2) if a range is set, keep only in-range points in xa, ya
      3) resolve the base temperature, fit with fit_q10()
      4) build the fitted-curve points (fcx, fcy)
      5) Arrhenius regression on the same data (arrhenius_fit)
    In the result, xa/ya are the data used for Q10 and x/y are all data.
    """
    x, y = np.asarray(temps, float), np.asarray(periods, float)
    if len(x) != len(y):
        raise ValueError(f"Number of temperatures ({len(x)}) and periods ({len(y)}) differ.")
    if len(x) < 3:
        raise ValueError("Enter at least 3 data points.")
    if np.any(y <= 0):
        raise ValueError("All periods must be positive.")

    if use_range:
        if not low < high:
            raise ValueError("Range low limit must be smaller than the high limit.")
        mask = (x >= low) & (x <= high)
        xa, ya = x[mask], y[mask]
    else:
        xa, ya = x, y
    if len(xa) < 3:
        raise ValueError(f"Only {len(xa)} data point(s) in the selected range; need at least 3.")

    base_t = resolve_base_temperature(base_sel, base_custom, xa)
    popt, pcov, r2, model = fit_q10(xa, ya, base_t)
    perr = np.sqrt(np.diag(pcov))

    # Fitted-curve range: data range extended by 33.3 % of its width on each side.
    margin = (x.max() - x.min()) * 0.333
    fcx = np.linspace(int(x.min() - margin), int(x.max() + margin) + 1, 200)
    fcy = model(fcx, *popt)

    arr = arrhenius_fit(xa, ya)
    # Reference curves fitted to ALL data. The plots take their axis limits from these
    # (plus all data points), so restricting the fit range never changes the axis scales.
    ref = dict(fcy=fcy, arr=arr)
    if use_range:
        try:
            popt_all, _, _, model_all = fit_q10(x, y, base_t)
            ref = dict(fcy=model_all(fcx, *popt_all), arr=arrhenius_fit(x, y))
        except (ValueError, RuntimeError):
            pass  # keep the fitted-range curves as reference

    return dict(
        x=x, y=y, xa=xa, ya=ya, use_range=use_range, low=low, high=high, base_t=base_t,
        popt=popt, perr=perr, r2=r2, fcx=fcx, fcy=fcy, arr=arr, ref=ref,
    )


# ======================================================================
# Report text (the Results tab and "Save results...")
# ======================================================================

def stats_text(arr, name, unit):
    """Minimum, maximum, average and median, one per line."""
    return (
        f"Minimum {name} = {np.min(arr):.3f} {unit}\n"
        f"Maximum {name} = {np.max(arr):.3f} {unit}\n"
        f"Average {name} = {statistics.mean(arr):.3f} {unit}\n"
        f"Median {name}  = {statistics.median(arr):.3f} {unit}\n"
    )


def table_text(x, y):
    """Fixed-width, right-aligned table of index, temperature and period."""
    lines = [f"{'#':>4}  {'Temperature':>12}  {'Period':>10}"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"{i:>4}  {xi:>12.3f}  {yi:>10.3f}")
    return "\n".join(lines)


def summary_line(r):
    """Two-line fit summary shown under the plot title (matplotlib mathtext)."""
    tau, q10 = r["popt"]
    tau_e, q10_e = r["perr"]
    return (f"Est. tau at {r['base_t']:.2f}{DEG}C = {tau:.3f} \u00B1 {tau_e:.3f} h\n"
            f"Q$_{{10}}$ = {q10:.3f} \u00B1 {q10_e:.3f}    R$^2$ = {r['r2']:.4f}")


def report_text(r, show_tab_delimited=False):
    """Assemble the full report.

      INPUT DATA -> (if a range is set) DATA USED FOR Q10 CALCULATION ->
      (if tab-delimited) dataset, fitted curve and Arrhenius regression line ->
      RESULTS -> ARRHENIUS REGRESSION
    The tab-delimited part can be pasted straight into a spreadsheet.
    """
    x, y, xa, ya = r["x"], r["y"], r["xa"], r["ya"]
    tau, q10 = r["popt"]
    tau_e, q10_e = r["perr"]
    arr = r["arr"]
    out = ["INPUT DATA", table_text(x, y), "", stats_text(x, "temperature", DEG + "C"),
           stats_text(y, "period", "h")]
    if r["use_range"]:
        out += ["DATA USED FOR Q10 CALCULATION", table_text(xa, ya), "",
                stats_text(xa, "temperature", DEG + "C"), stats_text(ya, "period", "h")]
    if show_tab_delimited:
        out += ["DATASET", f"Temp ({DEG}C)\tPeriod (hours)"]
        out += [f"{a:.3f}\t{b:.3f}" for a, b in zip(x, y)]
        out += ["", "FITTED CURVE", f"Temp ({DEG}C)\tPeriod (hours)"]
        out += [f"{a:.3f}\t{b:.3f}" for a, b in zip(r["fcx"], r["fcy"])]
        out += ["", "ARRHENIUS REGRESSION LINE", f"Temp ({DEG}C)\t1000/T (1/K)\tln(Frequency)"]
        out += [f"{t:.3f}\t{1000 / (t + 273.15):.5f}\t{arrhenius_lnf(arr, t):.5f}" for t in r["fcx"]]
        out.append("")
    out += [
        "RESULTS",
        f"Estimated period length at {r['base_t']:.2f} {DEG}C = {tau:.3f} \u00B1 {tau_e:.3f} hours",
        f"Q10 (temperature coefficient) = {q10:.3f} \u00B1 {q10_e:.3f}",
        f"R\u00B2 = {r['r2']:.6f}",
        "",
        "ARRHENIUS REGRESSION (ln(Frequency) vs. 1000/T, data used for Q10 calculation)",
        f"Activation energy Ea = {arr['Ea'] / 1000:.3f} \u00B1 {arr['se_Ea'] / 1000:.3f} kJ/mol",
        f"Slope d ln(f) / d(1000/T) = {arr['slope']:.4f} \u00B1 {arr['se_slope']:.4f} (10\u00B3 K)",
        f"Intercept = {arr['intercept']:.4f} \u00B1 {arr['se_int']:.4f}",
        f"R\u00B2 = {arr['r2']:.6f}",
    ]
    return "\n".join(out)


# ======================================================================
# Plotting (each function draws onto a matplotlib Axes)
#   period    : x = T (°C),            y = tau (h)
#   frequency : x = T (°C),            y = 24/tau (day^-1)
#   arrhenius : x = 1000/(T + 273.15), y = ln(24/tau)
# Every drawer takes (ax, r, show_arr); only the Arrhenius plot uses show_arr.
# ======================================================================

def _scatter(ax, r, xs_all, ys_all, xs_used, ys_used):
    """Data points: red = used in the fit, light green = outside the fit range."""
    if not r["use_range"]:
        ax.scatter(xs_used, ys_used, s=35, color="red", label="Data")
    else:
        ax.scatter(xs_all, ys_all, s=30, color="lightgreen", label="Data")
        ax.scatter(xs_used, ys_used, s=35, color="red",
                   label="Data points used for Q$_{10}$ calculation")


def _fix_limits(ax, pairs):
    """Set the axis limits to cover all (x, y) pairs plus a 5 % margin (matplotlib's
    default autoscale margin). The drawers pass all data points and the ALL-DATA
    reference curves, so the scales stay the same when a fit range is chosen;
    range-fit curves that run past these limits are clipped by the axes."""
    xs, ys = [], []
    for px, py in pairs:
        px, py = np.ravel(px).astype(float), np.ravel(py).astype(float)
        ok = np.isfinite(px) & np.isfinite(py)
        xs.append(px[ok]); ys.append(py[ok])
    xs, ys = np.concatenate(xs), np.concatenate(ys)
    if not len(xs):
        return
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    dx, dy = (x1 - x0) or 1.0, (y1 - y0) or 1.0
    ax.set_xlim(x0 - 0.05 * dx, x1 + 0.05 * dx)
    ax.set_ylim(y0 - 0.05 * dy, y1 + 0.05 * dy)


def _fit_label(r):
    return f"Fit   Q$_{{10}}$ = {r['popt'][1]:5.3f}"


def draw_period(ax, r, show_arr=False):
    _scatter(ax, r, r["x"], r["y"], r["xa"], r["ya"])
    ax.plot(r["fcx"], r["fcy"], "--", color="blue", label=_fit_label(r))
    _fix_limits(ax, [(r["x"], r["y"]), (r["fcx"], r["ref"]["fcy"])])
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel("Period (hours)")
    ax.set_title("Temperature vs. Period")
    ax.legend()


def draw_frequency(ax, r, show_arr=False):
    _scatter(ax, r, r["x"], 24 / r["y"], r["xa"], 24 / r["ya"])
    ax.plot(r["fcx"], 24 / r["fcy"], "--", color="blue", label=_fit_label(r))
    _fix_limits(ax, [(r["x"], 24 / r["y"]), (r["fcx"], 24 / r["ref"]["fcy"])])
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel(r"Frequency (day$^{-1}$)")
    ax.set_title("Temperature vs. Frequency")
    ax.legend()


def draw_arrhenius(ax, r, show_arr=False):
    """Arrhenius plot. With show_arr, also draw the Arrhenius regression line
    (solid grey) over the same temperature span as the Q10 curve."""
    inv = lambda t: 1000 / (np.asarray(t) + 273.15)
    ln_f = lambda p: np.log(24 / np.asarray(p))
    _scatter(ax, r, inv(r["x"]), ln_f(r["y"]), inv(r["xa"]), ln_f(r["ya"]))
    ax.plot(inv(r["fcx"]), ln_f(r["fcy"]), "--", color="blue", label=_fit_label(r))
    if show_arr:
        ax.plot(inv(r["fcx"]), arrhenius_lnf(r["arr"], r["fcx"]), "-", color="dimgray", lw=1.5,
                label=f"Arrhenius regression   E$_a$ = {r['arr']['Ea'] / 1000:.2f} kJ/mol")
    pairs = [(inv(r["x"]), ln_f(r["y"])), (inv(r["fcx"]), ln_f(r["ref"]["fcy"]))]
    if show_arr:
        pairs.append((inv(r["fcx"]), arrhenius_lnf(r["ref"]["arr"], r["fcx"])))
    _fix_limits(ax, pairs)
    ax.set_xlabel(r"1000/Temperature (K$^{-1}$)")
    ax.set_ylabel(r"ln(Frequency) (ln(day$^{-1}$))")
    ax.set_title("Arrhenius Plot: 1000/T vs. ln(Frequency)")
    ax.legend()


DRAWERS = [("Period", draw_period), ("Frequency", draw_frequency), ("Arrhenius", draw_arrhenius)]

# x-axis value <-> temperature (°C) for each plot; used for the range band and drag-to-select.
TO_X = {"Period": lambda t: t, "Frequency": lambda t: t,
        "Arrhenius": lambda t: 1000.0 / (t + 273.15)}
FROM_X = {"Period": lambda v: v, "Frequency": lambda v: v,
          "Arrhenius": lambda v: 1000.0 / v - 273.15}


def _estimate_lines(text, size_pt, width_in):
    """Rough count of rendered lines for wrapped text."""
    usable = width_in * 72 * 0.95
    n = 0
    for line in text.split("\n"):
        w = len(line) * 0.55 * size_pt
        n += max(1, math.ceil(w / usable))
    return n


def apply_header(fig):
    """(Re)place the plot title and the fit summary at the top of the figure.

    The title sits at the very top; the fit summary is placed clearly below it, and
    the axes are laid out beneath both. Sizes are computed in inches, so this is
    called again whenever the figure is resized.
    """
    hdr = getattr(fig, "_q10_header", None)
    if hdr is None:
        return
    title, summary, tsize, ssize = hdr
    for t in getattr(fig, "_q10_texts", []):
        try:
            t.remove()
        except Exception:
            pass
    fig._q10_texts = []

    w, h = fig.get_size_inches()
    y_in = 0.10  # top margin
    if title:
        n = _estimate_lines(title, tsize, w)
        fig._q10_texts.append(fig.text(0.5, 1 - y_in / h, title, ha="center", va="top",
                                       fontsize=tsize, wrap=True))
        y_in += n * tsize * 1.25 / 72 + 0.25  # extra gap -> summary sits lower
    fig._q10_texts.append(fig.text(0.5, 1 - y_in / h, summary, ha="center", va="top",
                                   fontsize=ssize))
    y_in += 2 * ssize * 1.3 / 72 + 0.12
    fig.tight_layout(rect=[0, 0, 1, max(0.3, 1 - y_in / h)])


def render_single(fig, drawer, r, title, show_arr=False, name=None):
    """Draw one plot (for a tab) with the title and fit summary; return its Axes.

    With name (the plot's DRAWERS name) and a fit range set, the range is shaded
    (on screen only; the saved figure has no band)."""
    fig.clear()
    ax = fig.add_subplot(111)
    drawer(ax, r, show_arr)
    if name and r["use_range"]:
        ax.axvspan(TO_X[name](r["low"]), TO_X[name](r["high"]), color="blue", alpha=0.07, lw=0, zorder=0)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE, 9)
    apply_header(fig)
    return ax


def make_span_selector(ax, onselect, onmove):
    """Horizontal SpanSelector, compatible with matplotlib 3.3-3.4 (rectprops) and 3.5+ (props)."""
    kw = dict(useblit=True, onmove_callback=onmove)
    style = dict(facecolor="tab:blue", alpha=0.2)
    try:
        return SpanSelector(ax, onselect, "horizontal", props=style, **kw)
    except TypeError:
        return SpanSelector(ax, onselect, "horizontal", rectprops=style, **kw)


def render_combined(fig, r, title, show_arr=False):
    """Figure for saving: title, summary, then the three plots stacked vertically."""
    fig.clear()
    for i, (_, drawer) in enumerate(DRAWERS, start=1):
        ax = fig.add_subplot(3, 1, i)
        drawer(ax, r, show_arr)
        ax.set_box_aspect(3 / 4)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE_COMBINED, 10)
    apply_header(fig)


# ======================================================================
# Help window ("How it works")
# HELP_CONTENT is a list of blocks rendered into a Tk Text widget:
#   ("h2", text)          section heading
#   ("h3", text)          sub-heading
#   ("p", segments)       paragraph
#   ("eq", segments)      displayed equation
#   ("ul"/"ol", [segments, ...])  bulleted / numbered list
#   ("table", [(term, segments), ...])  two-column list
# Segments are a list of plain strings and (style, text) pairs, where style is
# "b" (bold), "sub" (subscript), "sup" (superscript) or "link" (opens a URL).
# ======================================================================

Q10 = ["Q", ("sub", "10")]

HELP_CONTENT = [
    ("h2", "What Q10 means"),
    ("p", Q10 + [" is the factor by which a rate changes for a 10 °C rise in temperature. From rates R",
                 ("sub", "1"), " and R", ("sub", "2"), " measured at temperatures T", ("sub", "1"),
                 " and T", ("sub", "2"), ":"]),
    ("eq", Q10 + [" = (R", ("sub", "2"), " / R", ("sub", "1"), ")", ("sup", "10 / (T₂ − T₁)")]),
    ("p", ["Typical enzyme-catalysed reactions have "] + Q10 + [
        " of about 2 to 3. For a circadian clock, the “rate” is its frequency, 1/period. "
        "Circadian clocks are temperature-compensated, so their period changes only a little with "
        "temperature and "] + Q10 + [" is usually close to 1 (roughly 0.8 to 1.2)."]),
    ("table", [
        ("Q₁₀ > 1", ["Period shortens as temperature rises (the clock runs faster)."]),
        ("Q₁₀ = 1", ["Perfect temperature compensation."]),
        ("Q₁₀ < 1", ["Period lengthens as temperature rises (overcompensation)."]),
    ]),

    ("h2", "The model that is fitted"),
    ("p", ["With τ", ("sub", "b"), " the period at the base temperature T", ("sub", "b"), ":"]),
    ("eq", ["τ(T) = τ", ("sub", "b"), " · "] + Q10 + [("sup", "(Tb − T) / 10"), "        (1)"]),
    ("p", ["The frequency, in cycles per day, is f(T) = 24 / τ(T). Taking logarithms gives a "
           "straight line in T with slope ln("] + Q10 + [")/10:"]),
    ("eq", ["ln f(T) = ln(24 / τ", ("sub", "b"), ") + ((T − T", ("sub", "b"), ") / 10) · ln "] + Q10
     + ["        (2)"]),
    ("h3", "The base temperature"),
    ("p", ["The base temperature only decides at which temperature the fitted period τ", ("sub", "b"),
           " is reported. "] + Q10 + [", its standard error, R² and the fitted curve are the same "
           "whatever base you choose; only τ", ("sub", "b"), " and its error change. A base far "
           "outside your data, such as absolute zero (−273.15 °C), makes τ", ("sub", "b"),
           " a long extrapolation with a large error."]),

    ("h2", "The three plots"),
    ("table", [
        ("Period", ["Period (h) against temperature (°C), with the fitted curve from model (1)."]),
        ("Frequency", ["Frequency, 24/period (day", ("sup", "−1"), "), against temperature."]),
        ("Arrhenius", ["ln(frequency) against 1000/T, where T is in kelvin."]),
    ]),
    ("p", ["The Arrhenius equation k = A · exp(−E", ("sub", "a"), " / RT) gives a straight line of "
           "ln k against 1/T with slope −E", ("sub", "a"), "/R (E", ("sub", "a"),
           " is the activation energy, R = 8.314 J mol", ("sup", "−1"), " K", ("sup", "−1"),
           "). Model (2) is linear in T rather than 1/T, so the fitted curve is slightly bent on "
           "the Arrhenius plot. The two descriptions are approximately related by"]),
    ("eq", ["E", ("sub", "a"), " ≈ R · T", ("sub", "1"), " · T", ("sub", "2"), " · ln("] + Q10
     + [") / 10      (T", ("sub", "1"), ", T", ("sub", "2"), " in kelvin)"]),
    ("p", ["In every plot, red points are used in the fit, light-green points are outside the chosen "
           "temperature range, and the dashed line is the fitted "] + Q10 + [" curve."]),
    ("h3", "Arrhenius regression line"),
    ("p", ["The Arrhenius plot can also show the ordinary straight-line regression of ln(frequency) "
           "on 1000/T (solid grey line), fitted to the same data used for "] + Q10 + [":"]),
    ("eq", ["ln f = intercept + slope · (1000 / T)      E", ("sub", "a"), " = −slope · 1000 · R"]),
    ("p", ["Its E", ("sub", "a"), ", slope, intercept, standard errors and R² are listed in the "
           "Results. The two lines come from different models and are fitted to different residuals "
           "(period for "] + Q10 + [", ln f for the regression), so they nearly coincide within your "
           "data but separate slightly beyond it. Turn the line on or off with ",
           ("b", "Show Arrhenius regression line"), "."]),

    ("h2", "How the fit is calculated"),
    ("ol", [
        [("b", "Objective. "), "Non-linear least squares on the period itself (not its logarithm), "
         "minimising SSE = Σ (τ", ("sub", "i"), " − τ̂(T", ("sub", "i"), "))², using SciPy’s curve_fit."],
        [("b", "Starting values. "), "A straight-line regression of ln(period) on temperature, "
         "converted to τ", ("sub", "b"), " and "] + Q10 + [" through equation (2)."],
        [("b", "Standard errors. "), "From the covariance matrix s² · (JᵀJ)", ("sup", "−1"),
         " with s² = SSE / (n − 2), where J is the Jacobian with respect to τ", ("sub", "b"),
         " and "] + Q10 + [". The ± values are the square roots of its diagonal."],
        [("b", "Goodness of fit. "), "R² = 1 − SSE / SST, where SST = Σ (τ", ("sub", "i"),
         " − mean τ)²."],
        [("b", "Fitted curve. "), "Drawn with 200 points across the data’s temperature range, "
         "extended on each side by a third of its width."],
    ]),

    ("h2", "How to use it"),
    ("ol", [
        ["The program opens with a demo calculation of example data. Click ", ("b", "Clear"),
         " to start with your own."],
        ["Enter one temperature (°C) and period (h) pair per line, separated by commas, spaces, tabs "
         "or semicolons. Or use ", ("b", "Load CSV..."), " with temperature in column 1 and period in "
         "column 2. A header line and lines starting with # are skipped."],
        ["To fit only part of the data, tick ", ("b", "Use only data within a temperature range"),
         " and set Low and High (both ends included). Or, after calculating, simply drag across "
         "any plot: the range you select is filled in and the fit is recalculated at once. The "
         "current range is shaded on the plots, and the axis scales stay the same as with all "
         "data. Double-click a plot or click ", ("b", "Use all data"), " to go back to all data. "
         "While the toolbar’s zoom or pan mode is on, dragging zooms or pans instead."],
        ["Choose the base temperature. Minimum, Maximum, Average and Median are taken from the data "
         "used in the fit."],
        ["Click ", ("b", "Calculate"), ", or press Ctrl+Enter."],
        ["Switch between the Results, Period, Frequency and Arrhenius tabs. Each plot tab has a "
         "toolbar to zoom, pan and save that single plot."],
        [("b", "Save figure..."), " saves all three plots stacked in one PNG, PDF or SVG image. ",
         ("b", "Save results..."), " saves the results as a text file; tick ",
         ("b", "Also show tab-delimited data"), " to include data you can paste into a spreadsheet."],
    ]),

    ("h2", "Privacy and offline use"),
    ("ul", [
        [("b", "Calculations run on your computer. "), "The fit, the plots and the report are all "
         "computed locally by this program."],
        [("b", "Nothing is sent. "), "The program makes no network connections; the data you enter "
         "and the results are never sent anywhere."],
        [("b", "Nothing is stored. "), "Your data and results exist only while the program is "
         "running. Files are written only when you click a Save button, and only to the place you "
         "choose."],
        [("b", "Works offline. "), "No internet connection is needed at any point."],
        [("b", "Browser version. "), "A single-file HTML version of this calculator is available "
         "from ", ("link", HTML_URL), " (use the “Download raw file” button on that page). Clicking "
         "the link opens it in your web browser; that is the only time this program touches the "
         "internet, and only when you click it."],
    ]),
]


def open_url(url):
    """Open a URL in the default web browser (only ever called from a user click)."""
    try:
        webbrowser.open(url, new=2)
    except Exception:
        pass


class HelpWindow(tk.Toplevel):
    """Scrollable, read-only window that renders HELP_CONTENT. Closes with Esc or Close."""

    def __init__(self, master):
        super().__init__(master)
        self.title("How this calculator works")
        self.geometry(f"760x{min(720, master.winfo_screenheight() - 120)}")
        self.transient(master)

        base = tkfont.nametofont("TkDefaultFont")
        fam, size = base.actual("family"), max(10, int(base.actual("size")))
        fonts = {
            "body": (fam, size + 1),
            "b": (fam, size + 1, "bold"),
            "h2": (fam, size + 4, "bold"),
            "h3": (fam, size + 2, "bold"),
            "small": (fam, max(8, size - 1)),
        }

        bar = ttk.Frame(self, padding=(12, 8))
        bar.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(bar, text="How this calculator works", font=fonts["h2"]).pack(side=tk.LEFT)
        ttk.Button(bar, text="Close", command=self.destroy).pack(side=tk.RIGHT)
        ttk.Separator(self).pack(side=tk.TOP, fill=tk.X)

        box = ttk.Frame(self)
        box.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        txt = tk.Text(box, wrap="word", padx=18, pady=10, borderwidth=0, highlightthickness=0,
                      font=fonts["body"], cursor="arrow", spacing1=2, spacing3=2)
        sb = ttk.Scrollbar(box, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        txt.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Text styles used by the renderer.
        txt.tag_configure("b", font=fonts["b"])
        txt.tag_configure("h2", font=fonts["h2"], spacing1=16, spacing3=6)
        txt.tag_configure("h3", font=fonts["h3"], spacing1=10, spacing3=4)
        txt.tag_configure("p", spacing3=8)
        txt.tag_configure("sub", font=fonts["small"], offset=-3)
        txt.tag_configure("sup", font=fonts["small"], offset=5)
        txt.tag_configure("eq", lmargin1=24, lmargin2=24, background="#F1F4F8",
                          spacing1=6, spacing3=10)
        txt.tag_configure("li", lmargin1=12, lmargin2=34, spacing3=4)
        txt.tag_configure("term", font=fonts["b"], lmargin1=12, lmargin2=120)
        txt.tag_configure("link", foreground="#2F4FB5", underline=True)
        txt.tag_bind("link", "<Enter>", lambda e: txt.configure(cursor="hand2"))
        txt.tag_bind("link", "<Leave>", lambda e: txt.configure(cursor="arrow"))

        self._link_n = 0
        self._render(txt, HELP_CONTENT)
        txt.configure(state="disabled")

        self.bind("<Escape>", lambda e: self.destroy())
        self.focus_set()

    def _segments(self, txt, segs, block_tags):
        """Insert a list of segments, combining each segment's style with the block's tags."""
        for seg in segs:
            if isinstance(seg, str):
                txt.insert(tk.END, seg, block_tags)
                continue
            style, text = seg
            if style == "link":
                # Each link gets its own tag so the click knows which URL to open.
                self._link_n += 1
                tag = f"link{self._link_n}"
                txt.tag_bind(tag, "<Button-1>", lambda e, u=text: open_url(u))
                txt.insert(tk.END, text, block_tags + ("link", tag))
            else:
                txt.insert(tk.END, text, block_tags + (style,))

    def _render(self, txt, content):
        """Write every HELP_CONTENT block into the Text widget."""
        for kind, data in content:
            if kind in ("h2", "h3"):
                txt.insert(tk.END, data + "\n", (kind,))
            elif kind in ("p", "eq"):
                self._segments(txt, data, (kind,))
                txt.insert(tk.END, "\n", (kind,))
            elif kind in ("ul", "ol"):
                for i, item in enumerate(data, start=1):
                    txt.insert(tk.END, ("•  " if kind == "ul" else f"{i}.  "), ("li",))
                    self._segments(txt, item, ("li",))
                    txt.insert(tk.END, "\n", ("li",))
                txt.insert(tk.END, "\n")
            elif kind == "table":
                for term, segs in data:
                    txt.insert(tk.END, term + "\t", ("term",))
                    self._segments(txt, segs, ("li",))
                    txt.insert(tk.END, "\n", ("li",))
                txt.insert(tk.END, "\n")


# ======================================================================
# GUI
# User action -> run() -> analyze() -> plots and Results text are redrawn.
# The last result is kept in self.result and reused for redraws and saving.
# ======================================================================

class Q10App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Q10 calculator for circadian period")
        # Tall enough for the data box, but never taller than the screen.
        self.geometry(f"1250x{min(820, self.winfo_screenheight() - 80)}")
        self.minsize(1000, 560)
        self.result = None
        self.help_win = None
        self.selectors = []  # one SpanSelector per plot tab (references keep them alive)

        self.var_range = tk.BooleanVar(value=False)
        self.var_low = tk.StringVar(value="25")
        self.var_high = tk.StringVar(value="40")
        self.var_base = tk.StringVar(value=BASE_CHOICES[5])  # 37°C
        self.var_custom = tk.StringVar(value="37")
        self.var_tab = tk.BooleanVar(value=False)
        self.var_arr = tk.BooleanVar(value=True)   # show the Arrhenius regression line
        self.var_title = tk.StringVar()

        self._build_controls()
        self._build_output()
        self._toggle_states()

        self.bind("<Control-Return>", lambda e: self.run())

        # Start-up: open with the example already calculated and show the notice.
        self.load_example()
        self.run()
        self.status.set("Showing the demo calculation. Click Clear to start with your own data.")
        self._show_notice()

    # ---- layout ------------------------------------------------------
    def _build_controls(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)

        # Header: program name and the "How it works" button.
        head = ttk.Frame(left)
        head.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        ttk.Label(head, text="Q10 calculator", font=("TkDefaultFont", 13, "bold")).pack(side=tk.LEFT)
        ttk.Button(head, text="How it works", command=self.show_help).pack(side=tk.RIGHT)

        # Packing order matters: everything below the data box is packed to the BOTTOM
        # first, so on a short screen only the data box shrinks (it scrolls) while the
        # buttons stay visible.
        f = ttk.LabelFrame(left, text="Plot title", padding=6)
        f.pack(side=tk.TOP, fill=tk.X, pady=(0, 6))
        ttk.Entry(f, textvariable=self.var_title, width=44).pack(fill=tk.X)

        self.status = tk.StringVar(value="Ready.")
        ttk.Label(left, textvariable=self.status, foreground="#555", wraplength=330
                  ).pack(side=tk.BOTTOM, anchor="w", pady=(8, 0))

        row = ttk.Frame(left)
        row.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))
        ttk.Button(row, text="Save figure...", command=self.save_figure).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(row, text="Save results...", command=self.save_results).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(4, 0))

        ttk.Button(left, text="Calculate  (Ctrl+Enter)", command=self.run
                   ).pack(side=tk.BOTTOM, fill=tk.X, ipady=4)

        ttk.Checkbutton(left, text="Also show tab-delimited data in Results",
                        variable=self.var_tab).pack(side=tk.BOTTOM, anchor="w", pady=(0, 6))
        # Arrhenius regression line on/off: only the Arrhenius plot is redrawn.
        ttk.Checkbutton(left, text="Show Arrhenius regression line on the Arrhenius plot",
                        variable=self.var_arr, command=self._redraw_arrhenius
                        ).pack(side=tk.BOTTOM, anchor="w", pady=(0, 2))

        f = ttk.LabelFrame(left, text="Base temperature", padding=6)
        f.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(f, textvariable=self.var_base, values=BASE_CHOICES, state="readonly", width=34)
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._toggle_states())
        row = ttk.Frame(f)
        row.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(row, text="Custom value").pack(side=tk.LEFT)
        self.sp_custom = ttk.Spinbox(row, from_=-273.15, to=200, increment=0.1, width=8,
                                     textvariable=self.var_custom)
        self.sp_custom.pack(side=tk.LEFT, padx=4)
        ttk.Label(row, text=f"{DEG}C").pack(side=tk.LEFT)

        f = ttk.LabelFrame(left, text="Temperature range for the fit", padding=6)
        f.pack(side=tk.BOTTOM, fill=tk.X, pady=(0, 6))
        ttk.Checkbutton(f, text="Use only data within a temperature range",
                        variable=self.var_range, command=self._toggle_states).pack(anchor="w")
        row = ttk.Frame(f)
        row.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(row, text="Low").pack(side=tk.LEFT)
        self.sp_low = ttk.Spinbox(row, from_=-50, to=150, increment=0.1, width=8, textvariable=self.var_low)
        self.sp_low.pack(side=tk.LEFT, padx=(4, 12))
        ttk.Label(row, text="High").pack(side=tk.LEFT)
        self.sp_high = ttk.Spinbox(row, from_=-50, to=150, increment=0.1, width=8, textvariable=self.var_high)
        self.sp_high.pack(side=tk.LEFT, padx=4)
        ttk.Label(row, text=f"{DEG}C").pack(side=tk.LEFT)
        ttk.Button(row, text="Use all data", command=self.use_all_data).pack(side=tk.RIGHT)
        ttk.Label(f, text="Tip: after calculating, drag across a plot to choose the range; "
                          "double-click a plot to use all data.",
                  foreground="#555", wraplength=320, justify="left").pack(anchor="w", pady=(4, 0))

        # Data box: fills the remaining space.
        f = ttk.LabelFrame(left, text="Data (one pair per line)", padding=6)
        f.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(0, 6))
        ttk.Label(f, text=f"temperature({DEG}C)  period(h)  - separated by comma / space / tab"
                  ).pack(side=tk.TOP, anchor="w")
        row = ttk.Frame(f)
        row.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(row, text="Load CSV...", command=self.load_csv).pack(side=tk.LEFT)
        ttk.Button(row, text="Load example", command=self.load_example).pack(side=tk.LEFT, padx=4)
        ttk.Button(row, text="Clear", command=self.clear_data).pack(side=tk.LEFT)
        box = ttk.Frame(f)
        box.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(0, 4))
        self.txt_data = tk.Text(box, width=44, height=DATA_ROWS, wrap="none", undo=True)
        sb = ttk.Scrollbar(box, orient="vertical", command=self.txt_data.yview)
        self.txt_data.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_data.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # Typing in the data box hides the start-up notice.
        self.txt_data.bind("<Key>", lambda e: self._hide_notice(), add="+")

    def _build_output(self):
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Start-up notice (demo + privacy statement + link). Hidden by Clear,
        # editing the data, or loading a CSV.
        self.notice = tk.Frame(right, bg="#EEF2FB", highlightthickness=1,
                               highlightbackground="#C9D3EA", padx=12, pady=8)
        self.notice_labels = [
            tk.Label(self.notice, text=NOTICE_DEMO, bg="#EEF2FB", fg="#18233A",
                     justify="left", anchor="w"),
            tk.Label(self.notice, text=NOTICE_PRIVACY, bg="#EEF2FB", fg="#4B5563",
                     justify="left", anchor="w"),
        ]
        for i, lab in enumerate(self.notice_labels):
            lab.pack(side=tk.TOP, fill=tk.X, pady=(0 if i == 0 else 6, 0))
        # Link row: only the URL is clickable/underlined; the hint after it is plain text.
        link_row = tk.Frame(self.notice, bg="#EEF2FB")
        link_row.pack(side=tk.TOP, fill=tk.X)
        # padx=0 / bd=0 remove the labels' own padding, so the gap is just the one space
        # at the start of the hint text.
        link = tk.Label(link_row, text=HTML_URL, bg="#EEF2FB", fg="#2F4FB5", cursor="hand2",
                        font=("TkDefaultFont", 9, "underline"), padx=0, bd=0)
        link.pack(side=tk.LEFT)
        link.bind("<Button-1>", lambda e: open_url(HTML_URL))
        tk.Label(link_row, text=' (use the "Download raw file" button)', bg="#EEF2FB",
                 fg="#4B5563", font=("TkDefaultFont", 9), padx=0, bd=0).pack(side=tk.LEFT)
        # Re-wrap the notice text whenever its width changes.
        self.notice.bind("<Configure>", lambda e: [lab.configure(wraplength=max(200, e.width - 30))
                                                   for lab in self.notice_labels])

        self.nb = ttk.Notebook(right)
        self.nb.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        tab = ttk.Frame(self.nb)
        self.nb.add(tab, text="Results")
        self.txt_out = tk.Text(tab, wrap="none", font=("Courier", 10))
        ys = ttk.Scrollbar(tab, orient="vertical", command=self.txt_out.yview)
        xs = ttk.Scrollbar(tab, orient="horizontal", command=self.txt_out.xview)
        self.txt_out.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side=tk.RIGHT, fill=tk.Y)
        xs.pack(side=tk.BOTTOM, fill=tk.X)
        self.txt_out.pack(fill=tk.BOTH, expand=True)

        self.figs, self.canvases = [], []
        for name, _ in DRAWERS:
            tab = ttk.Frame(self.nb)
            self.nb.add(tab, text=name)
            fig = Figure(figsize=(6, 4.5), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.mpl_connect("resize_event", lambda ev, f=fig: self._relayout(f))
            canvas.mpl_connect("button_press_event", self._on_plot_press)
            NavigationToolbar2Tk(canvas, tab).update()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self.figs.append(fig)
            self.canvases.append(canvas)

    @staticmethod
    def _relayout(fig):
        try:
            apply_header(fig)
        except Exception:
            pass

    def _toggle_states(self):
        """Enable/disable Low, High and Custom according to the checkbox and selection."""
        rng = "normal" if self.var_range.get() else "disabled"
        self.sp_low.configure(state=rng)
        self.sp_high.configure(state=rng)
        self.sp_custom.configure(state="normal" if self.var_base.get() == CUSTOM_BASE else "disabled")

    # ---- start-up notice and help -----------------------------------
    def _show_notice(self):
        self.notice.pack(side=tk.TOP, fill=tk.X, pady=(0, 8), before=self.nb)

    def _hide_notice(self):
        self.notice.pack_forget()

    def show_help(self):
        """Open the "How it works" window (or bring the open one to the front)."""
        if self.help_win is not None and self.help_win.winfo_exists():
            self.help_win.lift()
            self.help_win.focus_set()
            return
        self.help_win = HelpWindow(self)

    # ---- data helpers ------------------------------------------------
    def _set_text(self, widget, text):
        widget.delete("1.0", tk.END)
        widget.insert("1.0", text)

    def load_example(self):
        self._set_text(self.txt_data, EXAMPLE_DATA)
        self.var_title.set("Example")

    def clear_data(self):
        """Clear the data box, the plot title, the Results text and all three plots."""
        self._hide_notice()
        self._set_text(self.txt_data, "")
        self.var_title.set("")
        self._set_text(self.txt_out, "")
        for fig, canvas in zip(self.figs, self.canvases):
            fig._q10_header = None  # stop the resize handler from redrawing the header
            fig.clear()
            canvas.draw()
        self.selectors = []
        self.result = None
        self.status.set("Cleared.")

    def load_csv(self):
        """Read CSV/TSV/text; keep rows whose first two columns are numbers.
        The rows are written into the data box so they can be checked before calculating."""
        path = filedialog.askopenfilename(
            title="Open CSV / TSV (column 1 = temperature, column 2 = period)",
            filetypes=[("Delimited text", "*.csv *.tsv *.txt"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, newline="", encoding="utf-8-sig") as fh:
                sample = fh.read(4096)
                fh.seek(0)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters=",;\t ")
                except csv.Error:
                    dialect = csv.excel
                temps, periods = [], []
                for row in csv.reader(fh, dialect):
                    row = [c.strip() for c in row if c.strip()]
                    if len(row) < 2:
                        continue
                    try:
                        t, p = float(row[0]), float(row[1])
                    except ValueError:
                        continue  # header line
                    temps.append(t)
                    periods.append(p)
            if not temps:
                raise ValueError("No numeric rows found (need two numeric columns).")
        except Exception as exc:
            messagebox.showerror("Load failed", str(exc))
            return
        self._hide_notice()
        self._set_text(self.txt_data, "\n".join(f"{t:g}, {p:g}" for t, p in zip(temps, periods)))
        self.status.set(f"Loaded {len(temps)} rows from {path.split('/')[-1]}")

    # ---- actions -----------------------------------------------------
    def run(self, select_tab=True):
        """Calculate: read inputs -> analyze() -> redraw plots and Results.

        Returns True on success. select_tab=False keeps the current tab (used when the
        range is chosen by dragging on a plot)."""
        try:
            temps, periods = parse_pairs(self.txt_data.get("1.0", tk.END))
            use_range = self.var_range.get()
            low = float(self.var_low.get()) if use_range else 0.0
            high = float(self.var_high.get()) if use_range else 0.0
            custom = float(self.var_custom.get()) if self.var_base.get() == CUSTOM_BASE else 0.0
            result = analyze(temps, periods, use_range, low, high, self.var_base.get(), custom)
        except (ValueError, RuntimeError) as exc:
            messagebox.showerror("Cannot calculate", str(exc))
            self.status.set("Error: " + str(exc))
            return False

        self.result = result
        self.selectors = [None] * len(DRAWERS)
        for i in range(len(DRAWERS)):
            self._render_plot(i)
        self._set_text(self.txt_out, report_text(result, self.var_tab.get()))
        tau, q10 = result["popt"]
        msg = (f"Q10 = {q10:.3f}, tau = {tau:.3f} h at {result['base_t']:.2f}{DEG}C, "
               f"R\u00B2 = {result['r2']:.4f}")
        if result["use_range"]:
            msg += f"  (fit range {result['low']:g}–{result['high']:g} {DEG}C, {len(result['xa'])} points)"
        self.status.set(msg)
        if select_tab:
            self.nb.select(1)
        return True

    def _render_plot(self, i):
        """Draw plot tab i (with the range band) and attach a fresh drag selector."""
        name, drawer = DRAWERS[i]
        ax = render_single(self.figs[i], drawer, self.result, self.var_title.get().strip(),
                           self.var_arr.get(), name=name)
        self.canvases[i].draw()
        self.selectors[i] = make_span_selector(
            ax, lambda a, b, n=name, ax=ax: self._on_span(n, ax, a, b),
            lambda a, b, n=name: self._on_span_move(n, a, b))

    def _redraw_arrhenius(self):
        """Redraw only the Arrhenius tab after the regression-line checkbox changes."""
        if self.result is None:
            return
        self._render_plot([name for name, _ in DRAWERS].index("Arrhenius"))

    # ---- drag-to-select temperature range ---------------------------
    # Dragging horizontally across a plot selects a temperature range (matplotlib
    # SpanSelector). While dragging, the status line shows the range; on release it is
    # written into Low/High, the range checkbox is ticked and the calculation runs again
    # (the current tab is kept). The x value is converted back to temperature (for the
    # Arrhenius plot T = 1000/x - 273.15); Low is rounded down and High up to 0.1 °C so
    # that points at the edges are included. If the new range cannot be used, the
    # previous settings and result are kept. A double-click switches the range off.
    # While the toolbar's zoom or pan mode is active, dragging zooms/pans instead.
    MIN_DRAG_PX = 4  # shorter drags count as clicks and are ignored

    @staticmethod
    def _span_temps(name, a, b):
        t1, t2 = FROM_X[name](a), FROM_X[name](b)
        return min(t1, t2), max(t1, t2)

    def _on_span_move(self, name, a, b):
        lo, hi = self._span_temps(name, a, b)
        self.status.set(f"Selected: {lo:.1f}–{hi:.1f} {DEG}C - release to recalculate")

    def _on_span(self, name, ax, a, b):
        px = ax.transData.transform([(a, 0), (b, 0)])[:, 0]
        if abs(px[1] - px[0]) < self.MIN_DRAG_PX:
            return
        lo, hi = self._span_temps(name, a, b)
        lo = math.floor(lo * 10 + 1e-6) / 10
        hi = math.ceil(hi * 10 - 1e-6) / 10
        # Defer: redrawing the figure inside the selector's own callback is unsafe.
        self.after_idle(lambda: self._apply_range(lo, hi))

    def _apply_range(self, lo, hi):
        """Set the fit range to [lo, hi] and recalculate; restore the previous settings if it fails."""
        prev = (self.var_range.get(), self.var_low.get(), self.var_high.get())
        self.var_range.set(True)
        self.var_low.set(f"{lo:.1f}")
        self.var_high.set(f"{hi:.1f}")
        self._toggle_states()
        if not self.run(select_tab=False):
            self.var_range.set(prev[0])
            self.var_low.set(prev[1])
            self.var_high.set(prev[2])
            self._toggle_states()

    def _on_plot_press(self, event):
        """Double-click inside a plot -> use all data (ignored while zoom/pan is active)."""
        if event.dblclick and event.inaxes is not None and not event.canvas.widgetlock.locked():
            self.after_idle(self.use_all_data)

    def use_all_data(self):
        """Switch the fit range off and recalculate with all data."""
        if not self.var_range.get() or self.result is None:
            return
        self.var_range.set(False)
        self._toggle_states()
        self.run(select_tab=False)

    def _need_result(self):
        if self.result is None:
            messagebox.showinfo("Nothing to save", "Run a calculation first.")
            return False
        return True

    def save_figure(self):
        """Save the three plots stacked in one image (PNG / PDF / SVG)."""
        if not self._need_result():
            return
        path = filedialog.asksaveasfilename(
            title="Save all three plots", defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not path:
            return
        fig = Figure(figsize=(8, 18), dpi=150)
        render_combined(fig, self.result, self.var_title.get().strip(), self.var_arr.get())
        fig.savefig(path)
        self.status.set(f"Saved {path}")

    def save_results(self):
        if not self._need_result():
            return
        path = filedialog.asksaveasfilename(
            title="Save results text", defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")])
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(report_text(self.result, self.var_tab.get()) + "\n")
        self.status.set(f"Saved {path}")


def main():
    # Crisp text on high-DPI Windows displays (no effect elsewhere).
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    Q10App().mainloop()


if __name__ == "__main__":
    main()
