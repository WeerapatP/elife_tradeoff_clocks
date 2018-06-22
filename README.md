# Biophysical clocks face a trade-off between internal and external noise resistance

This repository provides the code to do the simulations we have in the paper (Biophysical clocks face a trade-off between internal and external noise resistance). Folders are divided according to the figures that the code generates.

## Fig1

To use the [kaiABC.py](https://github.com/WeerapatP/elife_tradeoff_clocks/blob/master/Fig1/kaiABC.py "kaiABC.py"), try
```
python kaiABC.py [ATP_min] [ATP_max] [eta] [sfactor] [stoch_flag] [env_noise]
```

For example,
```
python kaiABC.py 0.2 0.8 1 1 0 0
```

The code will print out the MI value. The details of the values we use in the figure 1 can be found in the Appendix 1.

## Fig2

Files beginning with 'plot' will run their corresponding file beginning with 'simulate' multiple times to generate a plot similar to Figure 2 in the paper. For example, try running [plot_MI_vs_intnoise_millar.m](https://github.com/WeerapatP/elife_tradeoff_clocks/blob/master/Fig2/plot_MI_vs_intnoise_millar.m "plot_MI_vs_intnoise_millar.m") in MATLAB. This will run [simulate_millar.m](https://github.com/WeerapatP/elife_tradeoff_clocks/blob/master/Fig2/simulate_millar.m "simulate_millar.m") with different values of internal noise and parameters (which can be found in the Appendix 2 of the paper).

## Fig4, 6, Appendix 5 Fig 1 and 2
Similar to Fig2 folder, the files beginning with 'plot' will run their corresponding file beginning with 'simulate' multiple times to generate the plots. For example, [plot_MI_vs_days.m](https://github.com/WeerapatP/elife_tradeoff_clocks/blob/master/Fig4%2C%206%2C%20Appendix%205%20Fig%201%20and%202/plot_MI_vs_days.m "plot_MI_vs_days.m") will run [simulate_MI_vs_day.m](https://github.com/WeerapatP/elife_tradeoff_clocks/blob/master/Fig4%2C%206%2C%20Appendix%205%20Fig%201%20and%202/simulate_MI_vs_day.m "simulate_MI_vs_day.m") multiple times to generate the subfigure in figure 6. 

## Fig5 panel c and e:
First run the Matlab code “run_this.m” to generate the raw datas. Then the first step of making the figure is done by the Matlab code “priliminary_generation_of_plot.m”. With the results from the Matlab code, save the figure 1 as “priliminary_c.fig”, and save figure 2 as “priliminary_e.fig”. Then run the Matlab code “Final_generation_of_figure_5_c_and_e.m” to collect data from these plots and generate Figure 5 panel e and c. The final plots will be saved under the name “c.pdf” and “e.pdf”

## Fig5 panel d
First, run the Matlab code “run_this.m” to generate the raw data. Then generate the plot by running the Matlab code “code_for_plotting.m” which gives the plot for Figure 5 panel d. The plot can be saved into pdf format and the axis can be scaled as the reader wishes.