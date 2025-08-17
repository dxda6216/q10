# Circadian period Q<sub>10</sub> calculator
This is a simple code to calculate Q<sub>10</sub> values for circadian period length by using SciPy Optimize non-linear least squares fit on Google Colab.<br />

<img width="836" alt="q10cc" src="https://github.com/dxda6216/q10/blob/main/images/q10cc.jpg"><br />

This program is for calculating Q<sub>10</sub> values and plotting data only. Your data entered to this calculator will not be recorded or stored in any location (please see the code).

This calculator fits an exponential equation:

<img src="https://latex.codecogs.com/svg.image?\tau&space;_{t}=\frac{\tau&space;_{bt}}{{Q_{10}}^{^({\frac{t-bt}{10}})}">

to a set of data points, and estimates a period length at base temperature *bt* ˚C (τ<sub>*bt*</sub>) and temperature coefficient (Q<sub>10</sub>).

For mammalian cells data, it fits an exponential equation:

<img src="https://latex.codecogs.com/svg.image?\tau&space;_{t}=\frac{\tau&space;_{37}}{{Q_{10}}^{^({\frac{t-37}{10}})}">

to a set of data points, and estimates a period length at 37˚C (τ<sub>37</sub>) and temperature coefficient (Q<sub>10</sub>).

For cyanobacteria data, it fits an exponential equation:

<img src="https://latex.codecogs.com/svg.image?\tau&space;_{t}=\frac{\tau&space;_{30}}{{Q_{10}}^{^({\frac{t-30}{10}})}">

to a set of data points, and estimates a period length at 30˚C (τ<sub>30</sub>) and temperature coefficient (Q<sub>10</sub>).

## Usage

<b><ins>Open the calculator (Google Colab notebook)</ins></b>

Click on the link below:<br/>

https://colab.research.google.com/github/dxda6216/q10/blob/main/circadian_period_q10.ipynb<br/>

*For mammalian cells data, please use:*<br/>
https://colab.research.google.com/github/dxda6216/q10/blob/main/circadian_period_q10_for_mammalian_cells.ipynb<br/>

*For cyanobacteria data, please use:*<br/>
https://colab.research.google.com/github/dxda6216/q10/blob/main/circadian_period_q10_for_cyano.ipynb<br/>

Note: To use the Colab, you will need to login to your Google account.

<b><ins>Sample Description</ins></b>

To the "<i>Data_description</i>" input form, enter a description of your sample. This field can be blank.

<pre><b>Data_description:</b> <ins>Enter a description of your sample. </ins></pre>

<b><ins>Temperature and Period Data</ins></b>

Enter your data into the "<i>Temperatures</i>" and "<i>Periods</i>" input forms.
For example, if you have the following dataset

<pre><ins>Temp (˚C)</ins>   <ins>Period (hr)</ins>
  32.2        28.5
  35.3        27.1
  37.0        24.2
  39.6        22.6
  41.8        20.7</pre>
enter your data as follows:
<pre><b>Temperatures:</b> <ins>32.2, 35.3, 37.0, 39.6, 41.8     </ins>
<b>Periods:</b> <ins>28.5, 27.1, 24.2, 22.6, 20.7          </ins></pre>
The values should be delimited by a single comma. Leading and trailing spaces will be ignored.

<b><ins>Data Plotting</ins></b>

Scales of X-axis and Y-axis, ticks, etc will automatically be adjusted.

<b><ins>Exporting Dataset and Fitted Curve</ins></b>

To output the dataset and fitted curve in tab-delimited format, select "<b>Yes</b>" in the "<i>Display_tab_delimited_data</i>" dropdown, and 2-column datasets will be displayed. To export the data to Excel, copy the 2-column dataset and paste it onto an Excel worksheet (right click -> Paste Options: Paste Special... -> Paste as Unicode Text -> click "OK").

If you have period data obtained only at 2 different temperatures (only 2 data points), please see:
https://github.com/dxda6216/q10_two_temperatures
