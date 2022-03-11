# Circadian period Q<sub>10</sub> calculator
This is a simple code to calculate Q<sub>10</sub> values for circadian period length by using SciPy Optimize non-linear least squares fit on Google Colab.

**Open in Google Colab:**<br/>
https://colab.research.google.com/github/dxda6216/q10/blob/main/circadian_period_q10.ipynb

<img width="700" alt="q10ss" src="https://user-images.githubusercontent.com/101025597/157565673-5f3826cb-7808-4b62-8908-efda81ab1b93.png">

## Usage

This program is for calculating Q<sub>10</sub> values and plotting data only. Your data entered to this calculator will not be recorded or stored in any location (please see the code).

This calculator fits an exponential equation:<br />
![q10eq_s2](https://user-images.githubusercontent.com/101025597/157140650-b4d4112d-d1ea-4ffb-98b8-dbb24e2e45a7.png) <br />
to a set of data points, and estimates a period length at 37˚C and Q<sub>10</sub> value.

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

<b>(Optional)</b> To manually adjust the figure size, axis scales & ticks, etc, change the "<i>Set_graphical_parameters_of_figure</i>" dropdown field from "<b>No</b>" to "<b>Yes</b>", and modify arguments in the "<i>Figure_graphical_parameters</i>" input field. Those are for plotting only and don't affect the result of Q<sub>10</sub> calculation. The arguments (int or float values) are in the following order:
<pre><b>Figure_graphical_parameters:</b> <ins>8, 6, 27, 45, 28, 45, 2, 14, 45, 14, 45, 2, 28, 43     </ins>

 0th argument : Figure size - width
 1st argument : Figure size - height
 2nd argument : X-axis scale - min
 3rd argument : X-axis scale - max
 4th argument : X-axis ticks - starting
 5th argument : X-axis ticks - ending
 6th argument : X-axis ticks - interval
 7th argument : Y-axis scale - min
 8th argument : Y-axis scale - max
 9th argument : Y-axis ticks - starting
10th argument : Y-axis ticks - ending
11th argument : Y-axis ticks - interval
12th argument : Fitted curve - starting
13th argument : Fitted curve - ending</pre>

<b><ins>Exporting Dataset and Fitted Curve</ins></b>

To output the dataset and fitted curve in tab-delimited format, select "<b>Yes</b>" in the "<i>Display_tab_delimited_data</i>" dropdown, and 2-column datasets will be displayed. To export the data to Excel, copy the dataset and paste it onto an Excel worksheet (right click →  Paste Options: Paste Special... → Paste as Unicode Text → click "OK").

If you have period data obtained only at 2 different temperatures (only 2 data points), please see:
https://github.com/dxda6216/q10_two_temperatures
