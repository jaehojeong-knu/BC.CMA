# Bias-corrected estimation in causal mediation analysis (BC.CMA)

This repository provides reproducible R scripts for the manuscript **“Bias-corrected estimation in causal mediation analysis”** by *Jaeho Jeong, Jongho Im, and Young Min Kim*.

---

## Overview

The project introduces two likelihood-based bias-correction methods (BC1 and BC2) for causal mediation analysis (CMA).  
These methods address **transformation-induced finite-sample bias** in estimating the natural direct effect (NDE), natural indirect effect (NIE), and mediation proportion (MP).

---

## Project Structure

```
BC.CMA/
├── BC.CMA.Rproj
├── renv.lock                 # Reproducible environment snapshot
├── README.md
│
├── functions/                # Utility functions commonly used
│   ├── 00.basic_functions.R
│   ├── 01.CC_functions.R
│   ├── 02.CB_functions.R
│   ├── 03.BC_functions.R
│   └── 04.BB_functions.R
│
├── Sec4.Simulation/          # Section 4: Monte Carlo simulations
│   ├── CC_Simulation.R       
│   ├── CB_Simulation.R       
│   ├── BC_Simulation.R       
│   ├── BB_Simulation.R       
│   ├── Simulation_data/ 
│   └── Simulation_result/   
│
└── Sec5.Data_application/    # Section 5: Real-data applications
    ├── 01.BC_Data_application.R
    ├── 02.BB_Data_application.R
    └── Data_application_result/
```

---

## Quick Start
1. Clone  
   `git clone https://github.com/jaehojeong-knu/BC.CMA.git`

2. Open  
   Open `BC.CMA.Rproj` in RStudio (renv auto-activates).

3. Restore  
   ```r
   renv::restore()
   
---

## Setup & Reproducibility

This project uses **[`renv`](https://rstudio.github.io/renv/)** to ensure an identical R environment.

### 1. Clone the repository
```bash
git clone https://github.com/jaehojeong-knu/BC.CMA.git
cd BC.CMA
```

### 2. Restore the R environment
```r
install.packages("renv")
renv::restore()
```

### 3. Required packages
All dependencies are automatically handled by `renv`.  
Major packages include:
```r
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(psych)
library(boot)
library(medflex)
library(mediation)
```

---

## How to Run

### Simulation studies (Section 4)
Run each scenario script individually:
```r
source("Sec4.Simulation/CC_Simulation.R")  # Subsection 4.1
source("Sec4.Simulation/CB_Simulation.R")  # Subsection 4.2
source("Sec4.Simulation/BC_Simulation.R")  # Subsection 4.3
source("Sec4.Simulation/BB_Simulation.R")  # Subsection 4.4
```
Results are automatically saved under:
```
./Sec4.Simulation/Simulation_result/
```

### Real-data applications (Section 5)
```r
source("Sec5.Data_application/01.BC_Data_application.R")  # Subsection 5.1
source("Sec5.Data_application/02.BB_Data_application.R")  # Subsection 5.2
```
Results are automatically saved under:
```
./Sec5.Data_application/Data_application_result/
```


---

## License

This project is licensed under the **MIT License**.  
You are free to use, modify, and distribute this code with appropriate citation.

---

## Contact

**Jaeho Jeong**  
Ph.D. Candidate, Department of Statistics  
Kyungpook National University (KNU), South Korea  
E-mail: jaeho.jeong@knu.ac.kr  
