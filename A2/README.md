# ECE4860T14 Optimization - Assignment 2 (Winter 2021)

## Running and making all the plots

Run the following code with julia already installed in the root folder. The following line will install any necessary libraries and generate all the plots (in PNG), and all the output text files used to make the report. All these generated plots and outputs are stored in the assets folder.

```
julia --project=. ./src/main.jl
```

## LaTeX Code for the Report

The LaTeX code is found in lab_report/

The lab report with the inserted PNG files can be generated with the following command.

```
cd ./report/
pdflatex A2_LaberintoKim.tex
```


