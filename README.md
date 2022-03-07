# Circadian period Q<sub>10</sub> calculator
This is a simple code to calculate Q10 values for circadian period length by using SciPy Optimize non-linear least squares fit on Google Colab.

## Usage

This calculator is calculatig Q<sub>10</sub> values and plotting data only. Your data entered to this calculator will not be recorded or saved anywhere.

Enter your data into the "Temperatures" and "Periods" input forms.
For exampel, if you have the following dataset

<pre>Temp (˚C)   Period (hr)
  32.0        28.5
  37.0        24.2
  42.0        20.8</pre>
enter your data as follows:
<pre>Temperatures: 32.0, 37.0, 42.0         
Periods: 28.5, 24.2, 20.8              
The values should be delimited by a single comma
(leading and trailing spaces will be ignored).</pre>

To adjust figure size, scale of axis, etc. of the data plot, chenge the "Set_graphical_parameters_of_firure" dropdown field from "No" to "Yes", and modified arguments in the "Figure_graphical_parameters" input field. Those are for plotting only and don't affect the result of Q10 calculation. The arguments (int or float values) are in the following order:
<pre>Figure_graphical_parameters: 8, 6, 27, 45, 28, 45, 2, 14, 45, 14, 45, 2, 28, 43   

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
