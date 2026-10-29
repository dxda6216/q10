# Q10 calculator for circadian period

A small tool that estimates the temperature coefficient **Q<sub>10</sub>** of a circadian clock from period measurements at several temperatures. It fits

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Ctau%28T%29%20%3D%20%5Ctau_b%20%5Ccdot%20Q_%7B10%7D%5E%7B%28T_b%20-%20T%29%2F10%7D">
  <img alt="\tau(T) = \tau_b \cdot Q_{10}^{(T_b - T)/10}" src="https://latex.codecogs.com/svg.image?%5Ctau%28T%29%20%3D%20%5Ctau_b%20%5Ccdot%20Q_%7B10%7D%5E%7B%28T_b%20-%20T%29%2F10%7D">
</picture>
</p>

to your data and shows the result as Period, Frequency and Arrhenius plots. It also reports an Arrhenius regression line and the activation energy E<sub>a</sub>.

<p align="center">
<picture>
<img src="https://raw.githubusercontent.com/dxda6216/q10/refs/heads/main/images/q10_calculator_web_app_screenshot.png" alt="Q10 calculator web app screenshot" width="750">
</picture>
</p>

The same calculator comes in three versions. They use the same model, the same fitting method and the same inputs, and give the same numbers.

| | Browser (HTML) | Desktop (Python) | Google Colab |
|---|---|---|---|
| File | [`index.html`](https://github.com/dxda6216/q10/blob/main/index.html) | [`q10_calculator.py`](https://github.com/dxda6216/q10/blob/main/q10_calculator.py) | [`circadian_period_q10.ipynb`](https://github.com/dxda6216/q10/blob/main/circadian_period_q10.ipynb) |
| Runs on | Your computer, in the browser | Your computer | Google’s servers |
| Installation | None | Python with numpy, SciPy, matplotlib | None (Google account) |
| Internet | Not needed once downloaded | Not needed | Required |
| Your data leaves your computer | No | No | Yes, it is processed by Google Colab |
| Best for | Quick use, sharing with others, sensitive data | Offline work, PDF figures | Running without downloading anything |

---

## Contents

- [Quick start](#quick-start)
- [Background](#background)
- [How the fit is calculated](#how-the-fit-is-calculated)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Example data and expected results](#example-data-and-expected-results)
- [Error messages](#error-messages)
- [Privacy and data handling](#privacy-and-data-handling)
- [Version-specific instructions](#version-specific-instructions)
- [Differences between the versions](#differences-between-the-versions)
- [Use of AI](#use-of-ai)
- [License](#license)

---

## Quick start

**▶ Open the web app: https://dxda6216.github.io/q10/index.html**

The easiest way is to use the web application (browser version) available at [`https://dxda6216.github.io/q10/index.html`](https://dxda6216.github.io/q10/index.html). Although this app is hosted on GitHub Pages, the Q<sub>10</sub> calculations are performed within your web browser (on your computer), and neither the data you enter nor the analysis results are ever sent or stored anywhere.

### Browser version ### 

Open [`https://github.com/dxda6216/q10/blob/main/index.html`](https://github.com/dxda6216/q10/blob/main/index.html), click **Download raw file** (You can save the downloaded HTML file anywhere on your computer and rename it as you like, but please keep the .html file extension). Open the saved file in your web browser. It opens with a demo calculation. Click **Clear**, enter your data, then click **Calculate**.

Once the HTML file is downloaded into your computer, it can be opened in a web browser without an internet connection (offline), and the Q<sub>10</sub> calculations are performed within the web browser on your computer.

### Python version ###

Download [`q10_calculator.py`](https://github.com/dxda6216/q10/blob/main/q10_calculator.py) into your computer and execute it from the command line.

The .py version needs Python 3.9 or later and three packages: numpy, SciPy and matplotlib (version 3.3 or later). It also needs Tkinter, which comes with Python on most systems.

You can install the three packages with:
```
pip install numpy scipy matplotlib
```
and run:
```
python q10_calculator.py
```

### Colab version ###

Open [`https://github.com/dxda6216/q10/blob/main/circadian_period_q10.ipynb`](https://github.com/dxda6216/q10/blob/main/circadian_period_q10.ipynb) and click **Download raw file** to download the notebook file (`circadian_period_q10.ipynb`) into your local computer.

In [Google Colab](https://colab.research.google.com) (a Google account is required), choose **File → Upload notebook**, and select the saved notebook file (`circadian_period_q10.ipynb`).

Replace **Temperatures** and **Periods** in the form with your values, then choose **Runtime → Run all** (Ctrl+F9).

---

## Background

### What Q<sub>10</sub> means

Q<sub>10</sub> is the factor by which a rate changes for a 10 °C rise in temperature \([see Van\'t Hoff\'s classical Q<sub>10</sub> table](https://dxda6216.github.io/q10/images/vant_hoff_q10.png)\). From rates R<sub>1</sub> and R<sub>2</sub> measured at temperatures T<sub>1</sub> and T<sub>2</sub>:

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DQ_%7B10%7D%20%3D%20%5Cleft%28%5Cfrac%7BR_2%7D%7BR_1%7D%5Cright%29%5E%7B10%2F%28T_2%20-%20T_1%29%7D">
  <img alt="Q_{10} = \left(\frac{R_2}{R_1}\right)^{10/(T_2 - T_1)}" src="https://latex.codecogs.com/svg.image?Q_%7B10%7D%20%3D%20%5Cleft%28%5Cfrac%7BR_2%7D%7BR_1%7D%5Cright%29%5E%7B10%2F%28T_2%20-%20T_1%29%7D">
</picture>
</p>

Typical enzyme-catalysed reactions have Q<sub>10</sub> of about 2 to 3. For a circadian clock, the “rate” is its frequency, 1/period. Circadian clocks are *temperature-compensated*: their period changes only a little with temperature, so Q<sub>10</sub> is usually close to 1 (roughly 0.8 to 1.2).

| Q<sub>10</sub> | Meaning |
|---|---|
| > 1 | Period shortens as temperature rises (the clock runs faster). |
| = 1 | Perfect temperature compensation. |
| < 1 | Period lengthens as temperature rises (overcompensation). |

### The model that is fitted

With τ<sub>b</sub> the period at the base temperature T<sub>b</sub>:

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Ctau%28T%29%20%3D%20%5Ctau_b%20%5Ccdot%20Q_%7B10%7D%5E%7B%28T_b%20-%20T%29%2F10%7D%20%5Cqquad%20%281%29">
  <img alt="\tau(T) = \tau_b \cdot Q_{10}^{(T_b - T)/10} \qquad (1)" src="https://latex.codecogs.com/svg.image?%5Ctau%28T%29%20%3D%20%5Ctau_b%20%5Ccdot%20Q_%7B10%7D%5E%7B%28T_b%20-%20T%29%2F10%7D%20%5Cqquad%20%281%29">
</picture>
</p>

The frequency, in cycles per day, is f(T) = 24 / τ(T). Taking logarithms gives a straight line in T with slope ln(Q<sub>10</sub>)/10:

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cln%20f%28T%29%20%3D%20%5Cln%5Cfrac%7B24%7D%7B%5Ctau_b%7D%20%2B%20%5Cfrac%7BT%20-%20T_b%7D%7B10%7D%5C%2C%5Cln%20Q_%7B10%7D%20%5Cqquad%20%282%29">
  <img alt="\ln f(T) = \ln\frac{24}{\tau_b} + \frac{T - T_b}{10}\,\ln Q_{10} \qquad (2)" src="https://latex.codecogs.com/svg.image?%5Cln%20f%28T%29%20%3D%20%5Cln%5Cfrac%7B24%7D%7B%5Ctau_b%7D%20%2B%20%5Cfrac%7BT%20-%20T_b%7D%7B10%7D%5C%2C%5Cln%20Q_%7B10%7D%20%5Cqquad%20%282%29">
</picture>
</p>

The fitting yields an approximate solution for the coeficients τ<sub>b</sub> and Q<sub>10</sub>.

### The base temperature

The base temperature only decides at which temperature the fitted period τ<sub>b</sub> is reported. Q<sub>10</sub>, its standard error, R² and the fitted curve are the same whatever base you choose; only τ<sub>b</sub> and its error change. A base far outside your data, such as absolute zero (−273.15 °C), makes τ<sub>b</sub> a long extrapolation with a large error.

### The Arrhenius plot

The Arrhenius equation k = A · exp(−E<sub>a</sub> / RT) gives a straight line of ln k against 1/T, with slope −E<sub>a</sub>/R. E<sub>a</sub> is the activation energy and R = 8.314 J mol<sup>−1</sup> K<sup>−1</sup>. The calculator plots ln(frequency) against 1000/T, with T in kelvin.

Model (2) is linear in T, not in 1/T, so the fitted Q<sub>10</sub> curve is slightly bent on the Arrhenius plot. The two descriptions are approximately related by

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE_a%20%5Capprox%20%5Cfrac%7BR%5C%2CT_1%20T_2%20%5Cln%20Q_%7B10%7D%7D%7B10%7D%20%5Cqquad%20%28T_1%2C%20T_2%20%5Ctext%7B%20in%20kelvin%7D%29">
  <img alt="E_a \approx \frac{R\,T_1 T_2 \ln Q_{10}}{10} \qquad (T_1, T_2 \text{ in kelvin})" src="https://latex.codecogs.com/svg.image?E_a%20%5Capprox%20%5Cfrac%7BR%5C%2CT_1%20T_2%20%5Cln%20Q_%7B10%7D%7D%7B10%7D%20%5Cqquad%20%28T_1%2C%20T_2%20%5Ctext%7B%20in%20kelvin%7D%29">
</picture>
</p>

### The Arrhenius regression line

Optionally (on by default), the Arrhenius plot also shows the ordinary straight-line regression of ln(frequency) on 1000/T. It is fitted to the same data points used for Q<sub>10</sub>:

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cln%20f%20%3D%20%5Ctext%7Bintercept%7D%20%2B%20%5Ctext%7Bslope%7D%5Ccdot%5Cfrac%7B1000%7D%7BT%7D%2C%20%5Cqquad%20E_a%20%3D%20-%5Ctext%7Bslope%7D%5Ccdot%201000%20%5Ccdot%20R">
  <img alt="\ln f = \text{intercept} + \text{slope}\cdot\frac{1000}{T}, \qquad E_a = -\text{slope}\cdot 1000 \cdot R" src="https://latex.codecogs.com/svg.image?%5Cln%20f%20%3D%20%5Ctext%7Bintercept%7D%20%2B%20%5Ctext%7Bslope%7D%5Ccdot%5Cfrac%7B1000%7D%7BT%7D%2C%20%5Cqquad%20E_a%20%3D%20-%5Ctext%7Bslope%7D%5Ccdot%201000%20%5Ccdot%20R">
</picture>
</p>

The regression line and the Q<sub>10</sub> curve nearly coincide within the range of your data, but they are not identical, for two reasons:

1. The models differ in form: the Q<sub>10</sub> model is linear in T, the Arrhenius model is linear in 1/T.
2. They minimise different residuals: the Q<sub>10</sub> fit uses the period itself, the regression uses ln(frequency).

Outside the data range, the two lines separate slightly.

---

## How the fit is calculated

1. **Objective.** Non-linear least squares on the period itself (not its logarithm), minimising SSE = Σ (τ<sub>i</sub> − τ̂(T<sub>i</sub>))². The Python and Colab versions use SciPy’s `curve_fit`. The HTML version uses its own Levenberg–Marquardt implementation in JavaScript, which reaches the same solution.
2. **Starting values.** A straight-line regression of ln(period) on temperature, converted to τ<sub>b</sub> and Q<sub>10</sub> through equation (2).
3. **Standard errors.** From the covariance matrix s² · (JᵀJ)<sup>−1</sup>, with s² = SSE / (n − 2) and J the Jacobian with respect to τ<sub>b</sub> and Q<sub>10</sub>. The ± values are the square roots of its diagonal, the same definition as `curve_fit` with `absolute_sigma=False`.
4. **Goodness of fit.** R² = 1 − SSE / SST, where SST = Σ (τ<sub>i</sub> − mean τ)².
5. **Fitted curve.** Drawn with 200 points across the data’s temperature range, extended on each side by 33.3 % of its width. The end points are truncated to whole degrees.
6. **Arrhenius regression.** With x<sub>i</sub> = 1000 / (T<sub>i</sub> + 273.15) and y<sub>i</sub> = ln(24 / τ<sub>i</sub>):
   - slope = S<sub>xy</sub> / S<sub>xx</sub>
   - intercept = mean(y) − slope · mean(x)
   - SE(slope) = √(s² / S<sub>xx</sub>)
   - SE(intercept) = √(s² · (1/n + mean(x)² / S<sub>xx</sub>))
   - s² = SSE / (n − 2)
   - E<sub>a</sub> = −slope · 1000 · R, and SE(E<sub>a</sub>) = SE(slope) · 1000 · R, with R = 8.314462618 J mol<sup>−1</sup> K<sup>−1</sup>

---

## Inputs

All three versions take the same inputs. The names in brackets are the Colab form fields.

| Input | Description |
|---|---|
| Plot title (`Data_description`) | Optional title shown above the plots. |
| Data (`Temperatures`, `Periods`) | Temperatures in °C and periods in hours, one period per temperature. At least 3 points are needed. All periods must be positive. |
| Temperature range (`Set_temperature_range`, `Range_low_limit`, `Range_high_limit`) | Optionally use only points with Low ≤ T ≤ High for the fit. Both ends are included. Points outside the range are still plotted, in light green. |
| Base temperature (`Select_base_temperature`, `Base_temperature`) | See the options below. |
| Show Arrhenius regression line (`Show_Arrhenius_regression_line`) | Draws the regression line and its E<sub>a</sub> on the Arrhenius plot. The E<sub>a</sub> results are always reported. |
| Tab-delimited data (`Display_tab_delimited_data`) | Adds the dataset, the fitted curve and the Arrhenius regression line as tab-separated columns. These can be pasted straight into a spreadsheet. |

The base temperature options are:

- **Fixed values:** −273.15 °C (absolute zero), 0, 4, 25, 30, 37 or 100 °C.
- **Taken from the data used in the fit:** Minimum, Maximum, Average or Median. When a range is set, only the in-range points count.
- **Custom value:** any temperature you enter (the slider in Colab).

**Entering data.** In the HTML and Python versions, type one pair per line, for example `25.1, 33.6`. Values can be separated by commas, spaces, tabs or semicolons. Lines starting with `#` and one header line are skipped. Both versions can also load a CSV/TSV file with temperature in column 1 and period in column 2. In Colab, enter the temperatures and the periods as two comma-separated lists, in the same order.

---

## Outputs

### Results

The results contain:

- **Input data:** the input table, with the minimum, maximum, average and median of the temperatures and of the periods.
- **Data used for Q10 calculation:** the same table and statistics for the in-range points. This appears only when a range is set.
- **Tab-delimited data:** included only if you choose it.
- **RESULTS:**
  - the estimated period at the base temperature: τ<sub>b</sub> ± SE (hours)
  - Q<sub>10</sub> ± SE
  - R² of the Q<sub>10</sub> fit
- **ARRHENIUS REGRESSION:**
  - E<sub>a</sub> ± SE (kJ/mol)
  - slope ± SE
  - intercept ± SE
  - R² of the regression

### Plots

| Plot | x-axis | y-axis |
|---|---|---|
| Period | Temperature (°C) | Period (h) |
| Frequency | Temperature (°C) | Frequency, 24/period (day<sup>−1</sup>) |
| Arrhenius | 1000/T (K<sup>−1</sup>) | ln(frequency) |

In every plot, points used in the fit are red (orange in the HTML version). Points outside the chosen range are light green. The dashed line is the fitted Q<sub>10</sub> curve, and the solid grey line is the Arrhenius regression. The saved figure stacks all three plots under the title and a summary of τ<sub>b</sub>, Q<sub>10</sub> and R².

---

## Example data and expected results

All three versions open with the following example data. You can use it to check that your copy works.

| Temperature (°C) | 25.1 | 26.9 | 27.3 | 29.2 | 30 | 31.7 | 32.2 | 33.5 | 34.1 | 35.3 | 35.8 | 36.5 | 37.1 | 38.5 | 38.9 | 39.3 | 40.3 | 41.1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Period (h) | 33.6 | 32.8 | 30.4 | 29.3 | 28.0 | 27.9 | 26.1 | 25.5 | 24.9 | 24.8 | 24.5 | 24.4 | 24.1 | 23.9 | 23.6 | 21.9 | 21.5 | 20.1 |

**All 18 points, base temperature 37 °C:**

```
Estimated period length at 37.00 °C = 23.662 ± 0.221 hours
Q10 (temperature coefficient) = 1.327 ± 0.019
R² = 0.958656

Activation energy Ea = 21.859 ± 1.178 kJ/mol
Slope d ln(f) / d(1000/T) = -2.6290 ± 0.1416 (10³ K)
Intercept = 8.4911 ± 0.4612
R² = 0.955620
```

With a base temperature of 30 °C (the Colab default), the first line becomes 28.840 ± 0.232 hours. Everything else is unchanged.

**Range 28–40 °C (13 points used):** Q<sub>10</sub> = 1.261 ± 0.026, R² = 0.9197, E<sub>a</sub> = 18.051 ± 1.646 kJ/mol.

---

## Error messages

When the input can’t be used, no results are produced and a message explains why. The HTML version shows it in the status line, the Python version in a dialog box, and Colab in the cell output starting with “Cannot calculate:”.

| Message | What to do |
|---|---|
| Number of temperatures (*n*) and periods (*m*) differ. | Give exactly one period for each temperature. |
| Enter at least 3 data points. | Add more data. |
| All periods must be positive. | Check for zeros, negative values or typos. |
| Range low limit must be smaller than the high limit. | Swap or correct the range limits. |
| Only *n* data point(s) in the selected range; need at least 3. | Widen the range. |
| Line *n*: … (HTML and Python only) | That line doesn’t contain exactly two numbers; fix or delete it. |

---

## Privacy and data handling

**HTML and Python versions.** All calculations run on your own computer.

- **Nothing is sent.** The data you enter and the results are not sent anywhere.
- **Nothing is stored.** Data and results disappear when you close or reload the page, or quit the program. Files are written only when you click a Save button, and only to the place you choose.
- **Network use.** The Python version makes no network connections. The HTML version contacts only Google Fonts for its typefaces. That request contains none of your data, but like any web request it lets Google see your IP address and browser type. Offline, the page falls back to your system fonts and everything else works unchanged.
- **The GitHub link.** Both versions show a link to the download page. It opens only if you click it.

**Colab version.** The code runs in a Colab runtime on Google’s servers, so your data is processed by Google Colab.

- **What is saved with the notebook.** Form values become part of the notebook’s code, and results and plots are kept in its output. If you save the notebook (for example to Google Drive), they are saved with it.
- **Before sharing.** Remove your data from the form and choose **Edit → Clear all outputs**.
- **What the code itself does.** It makes no other network requests and writes no files.
- **Keeping data local.** To keep your data on your own computer, use the HTML or Python version.

---

## Version-specific instructions

### Browser version (`index.html`)

**Getting it.** Open [`https://github.com/dxda6216/q10/blob/main/index.html`](https://github.com/dxda6216/q10/blob/main/index.html) and click **Download raw file**. Then open the saved file in any modern browser. No installation or server is needed.

**Using it.**

1. The page opens with a demo calculation. A notice above the results says so and summarises how your data is handled. Click **Clear** to start; editing the data or loading a CSV also hides the notice.
2. Enter your data or use **Load CSV…**, set the options, and click **Calculate** (or press Ctrl+Enter / ⌘+Enter).
3. Switch between the **Period**, **Frequency**, **Arrhenius** and **Report** tabs.
4. **Save figure (PNG / SVG)** saves the three plots stacked in one image. **Save report (TXT)** saves the results text.

**How it works** (top left) opens an explanation of the theory, method, usage and privacy. The page follows your system’s light or dark mode and works on phones.

### Python version (`q10_calculator.py`)

**Requirements**

- [Python](https://www.python.org/) 3.9 or later
- numpy
- SciPy
- matplotlib 3.3 or later
- Tkinter, which is usually included with Python

Install the packages with:

```
pip install numpy scipy matplotlib
```

Whether Tkinter is included depends on how Python was installed:

| Python installation | Tkinter |
|---|---|
| Windows, macOS (python.org installers) | Included |
| macOS with Homebrew Python | `brew install python-tk` |
| Debian/Ubuntu | `sudo apt install python3-tk` (Fedora: `python3-tkinter`) |
| Anaconda | Everything is already included |
| Miniconda | `conda install numpy scipy matplotlib` |

To check that Tkinter works, run `python -m tkinter`. A small test window should appear.

**Anaconda notes**

With [Anaconda](https://anaconda.org/), everything the program needs is already installed, so no extra installation is required. The full Anaconda distribution includes Python, numpy, SciPy, matplotlib and Tkinter (conda’s tk package).

With Miniconda, it installs only a minimal set of packages, so add the three libraries:
  ```
  conda install numpy scipy matplotlib
  ```
Tkinter comes with conda’s Python, so it doesn’t need to be installed separately.

- **Separate environment (optional).** To keep this program apart from other projects:
  ```
  conda create -n q10 python=3.12 numpy scipy matplotlib
  conda activate q10
  ```
  Avoid mixing `pip install` and `conda install` in the same environment.
- **Spyder.** Running the program inside Spyder can make the window freeze or not appear. Run it from the Anaconda Prompt or a terminal instead, or set Run → Configuration per file → “Execute in an external system terminal”.
- **Jupyter.** This is a desktop window app and is not meant to run inside a notebook. For notebooks, use the `.ipynb` version.
- **Linux.** With older conda builds of Tk, text in the window can look jagged. This is cosmetic; `conda install -c conda-forge tk` may improve it.

**Using it.** Run:
```
python q10_calculator.py
```
The window opens with the demo calculation and a notice above the results; **Clear** removes both. Enter your data or use **Load CSV…**, set the options, and click **Calculate** (or press Ctrl+Enter). Each plot tab has a toolbar to zoom, pan and save that single plot. **Save figure…** saves the three plots stacked in one PNG, PDF or SVG file, and **Save results…** saves the results text. **How it works** opens an explanation window.

### Colab version (`circadian_period_q10.ipynb`)

**Getting it.** Open [`https://github.com/dxda6216/q10/blob/main/circadian_period_q10.ipynb`](https://github.com/dxda6216/q10/blob/main/circadian_period_q10.ipynb) and click **Download raw file**.

**Opening it.** In [Google Colab](https://colab.research.google.com), choose **File → Upload notebook**, and select the saved notebook file (`circadian_period_q10.ipynb`). A Google account is required.

**Using it.**

1. The notebook opens with the example data in the form and its result below. A notice at the top of the output says when the demo is shown.
2. Replace **Temperatures** and **Periods** with your values, as comma-separated lists in the same order.
3. Set the other form fields as needed, then choose **Runtime → Run all** (Ctrl+F9) or click the cell’s play button.
4. To save the figure, right-click it and choose **Save image as…**. To save the results, copy the printed text. Set **Display_tab_delimited_data** to “Yes” to get columns you can paste into a spreadsheet.

The **How this calculator works** section at the bottom of the notebook explains the theory, method, usage and privacy. The notebook also runs in Jupyter, where you edit the values directly in the code because the form fields are Colab-only.

---

## Differences between the versions

The calculations and results are the same in all three versions. Only the interface differs:

| | HTML | Python | Colab |
|---|---|---|---|
| Default base temperature | 37 °C | 37 °C | 30 °C |
| Data entry | Text box or CSV file | Text box or CSV file | Two comma-separated lists |
| Saving the figure | PNG, SVG | PNG, PDF, SVG | Right-click the image |
| Saving the results | TXT file | TXT file | Copy the output |
| Zoom and pan plots | No | Yes (plot toolbar) | No |
| Explanation | “How it works” pop-up | “How it works” window | Text section in the notebook |
| Demo notice disappears | Clear, editing data or loading a CSV | Clear, editing data or loading a CSV | When the data differs from the example |

---

## Use of AI

I wrote the core code for Q<sub>10</sub> calculation and graph generation in Python. The graphical user interface was coded with the assistance of generative AI, based on my design. Additionally, the HTML version was generated by AI based on the Python version. I reviewed the AI-generated code and made manual corrections as needed. Furthermore, this README was created with the assistance of AI, including the design and proofreading of the document. The user manual (PDF) was entirely AI-generated based on the README. The following is what the AI stated: \(Note that, in the quote, \"I,\" \"my\" and \"myself\" refer to the AI chatbot, and \"you\" and \"yours\" refer to the author of this app.\)

> I didn’t copy or reproduce code from any protected software. Everything I changed or added was written for your project, building on your original code. There are some third-party components the files depend on, but none of them were copied into your code.
> 
> **What I wrote, and what it’s based on**
> 
> The HTML version’s fitting code: I wrote the JavaScript Levenberg–Marquardt fit myself. It’s based on the published, standard algorithm, not ported from SciPy, MINPACK or any other library’s source code. The same goes for the other mathematics: the least-squares and standard-error formulas, the Arrhenius regression, and the “nice number” (1-2-5) rule for axis ticks. These are textbook methods.
>
> The Python and Colab changes: these extend your own code. The example data and the overall structure come from your files.
> 
> **Third-party components used, not copied**
> 
> numpy, SciPy, matplotlib, pandas and Tkinter: your original code already used these, and the updated versions still import them. They’re open source: numpy, SciPy and pandas use the BSD licence, matplotlib its own BSD-style licence, and Python and Tkinter the PSF licence. Importing them doesn’t put any licence obligations on your code. Obligations would apply only if you redistributed the libraries themselves.
> 
> Google Fonts (HTML version): the typefaces Instrument Sans and JetBrains Mono are under the SIL Open Font License. The page loads them from Google rather than including them in the file.
>
> CodeCogs (README): the formula images are generated and hosted by the CodeCogs service, not stored in your repository. If you rely on it publicly, check CodeCogs’ terms of use, which I haven’t verified.
> 
> Colab form syntax (#@title, #@param): this is a Colab feature that your original notebook already used.
> 
> **Tools that aren’t in your files**
> 
> To test the code, I used a browser automation tool (Playwright), ImageMagick, and SciPy’s linregress for comparison. None of these are part of the delivered files.
> 
> **One caveat**
> 
> I can’t give an absolute guarantee: my writing comes from patterns learned in training, and very common code idioms can resemble existing code. But nothing here was taken from a specific source, and the code consists of standard techniques written for your files.
> 
> If you publish the repository, it’s worth adding a LICENSE file (for example MIT or BSD-3) so others know how they may use your code. I’m not a lawyer, so for formal licensing questions it’s best to consult one.

## License

The Q10 calculator is released under the MIT License. The licence text is in the file [`LICENSE`](https://github.com/dxda6216/q10/blob/main/LICENSE) in the repository [`https://github.com/dxda6216/q10`](https://github.com/dxda6216/q10).
