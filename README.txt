README

This is a list of what each script accomplishes and which answers
they apply to for Prelim 1.
    1d.jl: Is used to answer an initial portion of question 1d.
        This code initially is used to convert <n> to a specific volume basis.
        This code then is used to find unknown variables W1 and W2.

    one_d.m: Is used to continue answering and finally plotting 1d.
        The variable values are directly copied from the results of 1d.jl
        A K_Inhibitor value was determined by doing a visual best fit on a semilogx plot
        The calculated data was then plotted against the actual data on a semilogx plot.

    2cplot.jl: This code is used to answer and plot for question 2c.
    This code finds the outer limits of S that lead to a three fixed-point system
    These S values are then utilized to produce a plot of all stable points of X vs. S
    A range from 0.0 to 2.6 is shown and the trend continues with no other known unstable points.

    timesol.jl: This code is used to answer and plot for question 2d.
        This produces a time-solved plot of the ODEs with S values of 0.02
        10.0, and 10^5. This produced three overlayed lines for X vs. time.

    bifurcationPoints.jl: Produces a plot of X over time. This code was used to change
    the S value to visualize at what S values does bifurcation occur. This code played
    a crucial role for which S values to utilize in Saddle.jl and Hopf.jl.

    Hopf.jl: This answers part 1 of question 2e. and provides a plot. A steady state
    at a value below the Hopf bifurcation point was selected and then modified for
    Cell B (+25% values) and Cell C (-25% values). The S value was then changed to 100
    causing the system to pass through the bifurcation point. This transition is plotted.

    Saddle.jl: This answers part 2 of question 2e. and provides a plot. A steady state
    at a value above the Saddle bifurcation point was selected and then modified for
    Cell B (+25% values) and Cell C (-25% values). The S value was then changed to 100
    causing the system to pass through the bifurcation point. This transition is plotted
