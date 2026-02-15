# 📦 monalisa_tutorial_ui

User-interface-based educational tutorial pipeline for **MONALISA** non-Cartesian MRI reconstruction.

This repository provides a structured, step-by-step workflow for:

- Coil sensitivity estimation  
- Data preparation (Mitosius wrapping)  
- Gridded reconstruction  
- Optional compressed sensing reconstruction  

It is designed for:

- 🎓 New lab members  
- 🧠 Researchers learning MONALISA  
- 🧪 Educational demonstrations  
- 🔬 Rapid reconstruction testing  

---

# 🔗 Official MONALISA Documentation

This repository is built on top of **MONALISA**.

📘 Full documentation:  
👉 https://mattechlab.github.io/monalisa/

Please refer to the official documentation for:

- API reference  
- Theory and algorithms  
- Advanced configuration  
- Supported data formats  
- Reconstruction backends  

This repository focuses on **guided usage via UI**, not core framework development.

---

# 🚀 Quick Start

## 1️⃣ Install Dependencies

- MONALISA  
- Pulseq  
- mapVBVD  

Example:

```matlab
addpath(genpath('/path/to/monalisa'));
addpath(genpath('/path/to/pulseq'));
addpath(genpath('/path/to/mapVBVD'));
```

## 2️⃣ Run the Educational Pipeline

```matlab
run main_script.m
```

The script will guide you through:

	1.	Optional RMS data inspection

	2.	Coil sensitivity estimation

	3.	Saving coil maps

	4.	Mitosius preparation
    
	5.	Final reconstruction (gridded or CS)