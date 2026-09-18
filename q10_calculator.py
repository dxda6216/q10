#!/usr/bin/env python3
"""
Q10 calculator for circadian period - standalone desktop app.

Fits   period(T) = tau_bt * Q10 ** ((base_T - T) / 10)   to temperature/period
data and shows the fit as Period-vs-Temperature, Frequency-vs-Temperature and
Arrhenius plots.

Requirements:  Python 3.9+, numpy, scipy, matplotlib   (tkinter ships with Python)
Run:           python q10_calculator.py
"""

import csv
import math
import re
import statistics
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

DEG = "\u00B0"

DATA_ROWS = 15          # visible lines in the data input box
TITLE_SIZE = 18          # plot-title font size (single-plot tabs)
TITLE_SIZE_COMBINED = 20  # plot-title font size (saved 3-plot figure)

EXAMPLE_TEMPS = "25.1, 26.9, 27.3, 29.2, 30, 31.7, 32.2, 33.5, 34.1, 35.3, 35.8, 36.5, 37.1, 38.5, 38.9, 39.3, 40.3, 41.1"
EXAMPLE_PERIODS = "33.6, 32.8, 30.4, 29.3, 28.0, 27.9, 26.1, 25.5, 24.9, 24.8, 24.5, 24.4, 24.1, 23.9, 23.6, 21.9, 21.5, 20.1"

EXAMPLE_DATA = "\n".join(
    f"{t.strip()}, {p.strip()}"
    for t, p in zip(EXAMPLE_TEMPS.split(","), EXAMPLE_PERIODS.split(","))
)

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


# ======================================================================
# Calculation core (no GUI dependencies)
# ======================================================================

def parse_numbers(text):
    """Parse numbers separated by commas, spaces, semicolons or newlines."""
    tokens = [t for t in re.split(r"[,\s;]+", text.strip()) if t]
    try:
        return np.array([float(t) for t in tokens], dtype=float)
    except ValueError as exc:
        raise ValueError(f"Could not read a number: {exc}") from None


def parse_pairs(text):
    """Parse one 'temperature period' pair per line.

    Separators may be commas, spaces, tabs, semicolons, or any mixture.
    Blank lines and lines starting with '#' are ignored; a non-numeric first
    line is treated as a header and skipped.
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
    def q10_model(temperature, tau_bt, q10):
        return tau_bt * q10 ** ((base_t - temperature) * 0.1)
    return q10_model


def initial_guess(xa, ya, base_t):
    """Start from a log-linear regression (more robust than a fixed (24, 1))."""
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
    model = make_q10_model(base_t)
    p0 = initial_guess(xa, ya, base_t)
    popt, pcov = curve_fit(model, xa, ya, p0=p0, maxfev=20000)
    residuals = ya - model(xa, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((ya - np.mean(ya)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return popt, pcov, r2, model


def stats_text(arr, name, unit):
    return (
        f"Minimum {name} = {np.min(arr):.3f} {unit}\n"
        f"Maximum {name} = {np.max(arr):.3f} {unit}\n"
        f"Average {name} = {statistics.mean(arr):.3f} {unit}\n"
        f"Median {name}  = {statistics.median(arr):.3f} {unit}\n"
    )


def table_text(x, y):
    lines = [f"{'#':>4}  {'Temperature':>12}  {'Period':>10}"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"{i:>4}  {xi:>12.3f}  {yi:>10.3f}")
    return "\n".join(lines)


def analyze(temps, periods, use_range, low, high, base_sel, base_custom):
    """Run the whole calculation and return a dict with everything needed to plot/report."""
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

    margin = (x.max() - x.min()) * 0.333
    fcx = np.linspace(int(x.min() - margin), int(x.max() + margin) + 1, 200)
    fcy = model(fcx, *popt)

    return dict(
        x=x, y=y, xa=xa, ya=ya, use_range=use_range, base_t=base_t,
        popt=popt, perr=perr, r2=r2, fcx=fcx, fcy=fcy,
    )


def summary_line(r):
    tau, q10 = r["popt"]
    tau_e, q10_e = r["perr"]
    return (f"Est. tau at {r['base_t']:.2f}{DEG}C = {tau:.3f} \u00B1 {tau_e:.3f} h\n"
            f"Q$_{{10}}$ = {q10:.3f} \u00B1 {q10_e:.3f}    R$^2$ = {r['r2']:.4f}")


def report_text(r, show_tab_delimited=False):
    x, y, xa, ya = r["x"], r["y"], r["xa"], r["ya"]
    tau, q10 = r["popt"]
    tau_e, q10_e = r["perr"]
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
        out.append("")
    out += [
        "RESULTS",
        f"Estimated period length at {r['base_t']:.2f} {DEG}C = {tau:.3f} \u00B1 {tau_e:.3f} hours",
        f"Q10 (temperature coefficient) = {q10:.3f} \u00B1 {q10_e:.3f}",
        f"R\u00B2 = {r['r2']:.6f}",
    ]
    return "\n".join(out)


# ======================================================================
# Plotting (each function draws onto a matplotlib Axes)
# ======================================================================

def _scatter(ax, r, xs_all, ys_all, xs_used, ys_used):
    if not r["use_range"]:
        ax.scatter(xs_used, ys_used, s=35, color="red", label="Data")
    else:
        ax.scatter(xs_all, ys_all, s=30, color="lightgreen", label="Data")
        ax.scatter(xs_used, ys_used, s=35, color="red",
                   label="Data points used for Q$_{10}$ calculation")


def _fit_label(r):
    return f"Fit   Q$_{{10}}$ = {r['popt'][1]:5.3f}"


def draw_period(ax, r):
    _scatter(ax, r, r["x"], r["y"], r["xa"], r["ya"])
    ax.plot(r["fcx"], r["fcy"], "--", color="blue", label=_fit_label(r))
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel("Period (hours)")
    ax.set_title("Temperature vs. Period")
    ax.legend()


def draw_frequency(ax, r):
    _scatter(ax, r, r["x"], 24 / r["y"], r["xa"], 24 / r["ya"])
    ax.plot(r["fcx"], 24 / r["fcy"], "--", color="blue", label=_fit_label(r))
    ax.set_xlabel(f"Temperature ({DEG}C)")
    ax.set_ylabel(r"Frequency (day$^{-1}$)")
    ax.set_title("Temperature vs. Frequency")
    ax.legend()


def draw_arrhenius(ax, r):
    inv = lambda t: 1000 / (np.asarray(t) + 273.15)
    ln_f = lambda p: np.log(24 / np.asarray(p))
    _scatter(ax, r, inv(r["x"]), ln_f(r["y"]), inv(r["xa"]), ln_f(r["ya"]))
    ax.plot(inv(r["fcx"]), ln_f(r["fcy"]), "--", color="blue", label=_fit_label(r))
    ax.set_xlabel(r"1000/Temperature (K$^{-1}$)")
    ax.set_ylabel(r"ln(Frequency) (ln(day$^{-1}$))")
    ax.set_title("Arrhenius Plot: 1000/T vs. ln(Frequency)")
    ax.legend()


DRAWERS = [("Period", draw_period), ("Frequency", draw_frequency), ("Arrhenius", draw_arrhenius)]


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


def render_single(fig, drawer, r, title):
    fig.clear()
    ax = fig.add_subplot(111)
    drawer(ax, r)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE, 9)
    apply_header(fig)


def render_combined(fig, r, title):
    fig.clear()
    for i, (_, drawer) in enumerate(DRAWERS, start=1):
        ax = fig.add_subplot(3, 1, i)
        drawer(ax, r)
        ax.set_box_aspect(3 / 4)
    fig._q10_header = (title, summary_line(r), TITLE_SIZE_COMBINED, 10)
    apply_header(fig)


# ======================================================================
# GUI
# ======================================================================

class Q10App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Q10 calculator for circadian period")
        # Tall enough for the data box, but never taller than the screen.
        self.geometry(f"1250x{min(820, self.winfo_screenheight() - 80)}")
        self.minsize(1000, 560)
        self.result = None

        self.var_range = tk.BooleanVar(value=False)
        self.var_low = tk.StringVar(value="25")
        self.var_high = tk.StringVar(value="40")
        self.var_base = tk.StringVar(value=BASE_CHOICES[5])  # 37°C
        self.var_custom = tk.StringVar(value="37")
        self.var_tab = tk.BooleanVar(value=False)
        self.var_title = tk.StringVar()

        self._build_controls()
        self._build_output()
        self._toggle_states()

        self.bind("<Control-Return>", lambda e: self.run())

    # ---- layout ------------------------------------------------------
    def _build_controls(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)

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

        # Data box: fills the remaining space, 25 lines tall by default.
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

    def _build_output(self):
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.nb = ttk.Notebook(right)
        self.nb.pack(fill=tk.BOTH, expand=True)

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
        rng = "normal" if self.var_range.get() else "disabled"
        self.sp_low.configure(state=rng)
        self.sp_high.configure(state=rng)
        self.sp_custom.configure(state="normal" if self.var_base.get() == CUSTOM_BASE else "disabled")

    # ---- data helpers ------------------------------------------------
    def _set_text(self, widget, text):
        widget.delete("1.0", tk.END)
        widget.insert("1.0", text)

    def load_example(self):
        self._set_text(self.txt_data, EXAMPLE_DATA)
        self.var_title.set("Example")

    def clear_data(self):
        """Clear the data box, the plot title, the Results text and all three plots."""
        self._set_text(self.txt_data, "")
        self.var_title.set("")
        self._set_text(self.txt_out, "")
        for fig, canvas in zip(self.figs, self.canvases):
            fig._q10_header = None  # stop the resize handler from redrawing the header
            fig.clear()
            canvas.draw()
        self.result = None
        self.status.set("Cleared.")

    def load_csv(self):
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
        self._set_text(self.txt_data, "\n".join(f"{t:g}, {p:g}" for t, p in zip(temps, periods)))
        self.status.set(f"Loaded {len(temps)} rows from {path.split('/')[-1]}")

    # ---- actions -----------------------------------------------------
    def run(self):
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
            return

        self.result = result
        title = self.var_title.get().strip()
        for fig, canvas, (_, drawer) in zip(self.figs, self.canvases, DRAWERS):
            render_single(fig, drawer, result, title)
            canvas.draw()
        self._set_text(self.txt_out, report_text(result, self.var_tab.get()))
        tau, q10 = result["popt"]
        self.status.set(f"Q10 = {q10:.3f}, tau = {tau:.3f} h at {result['base_t']:.2f}{DEG}C, "
                        f"R\u00B2 = {result['r2']:.4f}")
        self.nb.select(1)

    def _need_result(self):
        if self.result is None:
            messagebox.showinfo("Nothing to save", "Run a calculation first.")
            return False
        return True

    def save_figure(self):
        if not self._need_result():
            return
        path = filedialog.asksaveasfilename(
            title="Save all three plots", defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not path:
            return
        fig = Figure(figsize=(8, 18), dpi=150)
        render_combined(fig, self.result, self.var_title.get().strip())
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
