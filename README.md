# Bias-corrected estimation in causal mediation analysis

This repository provides reproducible R scripts for the manuscript **“Bias-corrected estimation in causal mediation analysis”** by *Jaeho Jeong, Jongho Im, and Young Min Kim*.

---


## Overview


It contains all R codes used in the manuscript, including Monte Carlo simulations (Section 4) and real-data applications (Section 5). It introduces two likelihood-based bias-correction methods (BC1 and BC2) for causal mediation analysis (CMA). These methods address **transformation-induced finite-sample bias** in estimating the natural direct effect (NDE), natural indirect effect (NIE), and mediation proportion (MP).

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

## Setup & Reproducibility

This project uses **[`renv`](https://rstudio.github.io/renv/)** to ensure a fully reproducible R environment.

---

### **1. Download the Repository**

#### **Option 1 : Clone via Git**
```bash
git clone https://github.com/jaehojeong-knu/BC.CMA.git
cd BC.CMA
```

#### **Option 2: Download as ZIP**
1. On the GitHub page, click **“Code” → “Download ZIP.”**  
2. Unzip the file and keep all contents together under a single directory, for example:
   ```
   BC.CMA-main/
   ├── BC.CMA.Rproj
   ├── functions/
   ├── Sec4.Simulation/
   └── Sec5.Data_application/
   ```
3. Open **`BC.CMA.Rproj`** in RStudio.  
   This automatically sets the project root so that all relative paths work correctly.

<hr style="border:0.5px solid #d3d3d3;">

### **2. Restore Required Packages**

All dependencies are managed through `renv`.  
To recreate the same environment used in the manuscript, run:
```r
install.packages("renv")
renv::restore()
```

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

<hr style="border:0.5px solid #d3d3d3;">

### **3. Run Simulation Studies (Section 4)**

Each simulation scenario can be executed independently:
```r
source("Sec4.Simulation/CC_Simulation.R")  # Subsection 4.1
source("Sec4.Simulation/CB_Simulation.R")  # Subsection 4.2
source("Sec4.Simulation/BC_Simulation.R")  # Subsection 4.3
source("Sec4.Simulation/BB_Simulation.R")  # Subsection 4.4
```

Simulation results are automatically saved under:
```
./Sec4.Simulation/Simulation_result/
```

<hr style="border:0.5px solid #d3d3d3;">

### **4. Run Real-Data Applications (Section 5)**

Run the following scripts for real-data analyses:
```r
source("Sec5.Data_application/01.BC_Data_application.R")  # Subsection 5.1
source("Sec5.Data_application/02.BB_Data_application.R")  # Subsection 5.2
```

Outputs are automatically saved under:
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
