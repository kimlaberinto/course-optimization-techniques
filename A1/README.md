# ECE4860T14 Optimization - Assignment 1 (Winter 2021)

## Running and making all the plots

Run the following code with julia already installed in the root folder. The following line will install any necessary libraries and generate all the plots (in SVG), and all the output text files used to make the report. All these generated plots and outputs are stored in assets/

```
julia --project=. ./src/makeplots.jl
```

## LaTeX Code for the Report

The LaTeX code is found in lab_report/

Note that the SVG plots are inserted in the LaTeX report using the svg package. The lab report with the inserted SVG files can be generated with the following command while inside the A1/lab_report/ folder

```
pdflatex --shell-escape A1_LaberintoKim.tex
```

The --shell-escape flag is necessary for the svg package.


