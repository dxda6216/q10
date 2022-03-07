# Circadian period Q<sub>10</sub> calculator
This is a simple code to calculate Q<sub>10</sub> values for circadian period length by using SciPy Optimize non-linear least squares fit on Google Colab.

https://colab.research.google.com/github/dxda6216/q10/blob/main/circadian_period_q10.ipynb

![q10c_ss_80p](https://user-images.githubusercontent.com/101025597/157003902-737a4096-d741-494e-b4b9-74103bb56a45.png)

## Usage

This calculator is calculatig Q<sub>10</sub> values and plotting data only. Your data entered to this calculator will not be recorded or saved anywhere.

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

<b>(Optional)</b> To manually adjust the figure size, axis scales & ticks, etc, chenge the "<i>Set_graphical_parameters_of_firure</i>" dropdown field from "<b>No</b>" to "<b>Yes</b>", and modify arguments in the "<i>Figure_graphical_parameters</i>" input field. Those are for plotting only and don't affect the result of Q<sub>10</sub> calculation. The arguments (int or float values) are in the following order:
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
