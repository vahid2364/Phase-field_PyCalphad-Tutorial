# **Phase-field_PyCalphad-Tutorial**
### **Dynamic Coupling of Microstructure Modeling with the Phase-field Method and CALPHAD Using PyCalphad**

---

## 🚀 **Overview**
Welcome to the **Phase-field_PyCalphad-Tutorial**, a comprehensive resource for exploring the dynamic coupling of microstructure evolution with thermodynamic modeling. This tutorial provides a practical approach to combining the **phase-field method** with the **PyCalphad library**, enabling predictive modeling of microstructural transformations in materials.

Developed as part of the *Texas A&M University MSEN 210 Course - Thermodynamics of Materials*, this tutorial serves as a stepping stone for students and researchers interested in computational materials science.

---

## 📚 **Introduction to Microstructure Thermodynamics**
Understanding the thermodynamics of microstructures is a cornerstone of materials science. This tutorial bridges the theoretical foundations with practical implementation using cutting-edge computational tools:
- **Phase-field Method**: Simulates the evolution of microstructures during processes such as solidification, phase transformations, and grain growth.
- **PyCalphad Library**: A Python-based library for thermodynamic calculations based on the CALPHAD approach, seamlessly integrated with phase-field simulations.

This repository is ideal for students, educators, and researchers aiming to deepen their understanding of microstructure thermodynamics.

---

## 🔧 **Instructions**
Follow the steps below to get started with the tutorial:

### 1️⃣ **Explore the Jupyter Notebook**
   - Start with the provided **Jupyter notebook**, which introduces the concepts and walks you through the dynamic coupling process.
   - The notebook includes visualizations and explanations to guide your understanding.

### 2️⃣ **Review the PowerPoint Slides**
   - Supplement your learning with the provided **PowerPoint presentation**, which summarizes key concepts and methodologies.

### 3️⃣ **Run the Phase-field Code**
   - The phase-field code is designed to run within the Jupyter notebook or from the terminal. Here’s how to execute it:
     - Ensure you have **gfortran** or **Intel Fortran** pre-installed on your system.
     - Open a terminal and navigate to the directory containing the code.
     - Execute the provided **Makefile**:
       ```bash
       make
       ```
     - This will compile the code and create an executable file called `main`.
     - Run the executable by typing:
       ```bash
       ./main
       ```
     - Wait 20-30 seconds for the simulation to complete.

### 4️⃣ **Visualize Microstructure Results**
   - The resulting microstructure files are stored in the `microstructure` folder in `.plt` format.
   - You can visualize these results using:
     - **Jupyter Notebook**: Directly load and analyze the `.plt` files.
     - **Tecplot**: A third-party (paid) visualization software.
     - **ParaView**: A free and open-source visualization tool.

---

## 💻 **Requirements**
- **Programming Environment**: Python 3.x, Jupyter Notebook
- **Libraries**:
  - `PyCalphad`
  - `Numpy`
  - `Matplotlib`
- **Compilers**: gfortran or Intel Fortran

---


## 🙋 Support

If you encounter any issues or have questions, feel free to reach out:
	•	Email: attari.v@tamu.edu
	•	Website: [vahid2364.tamu.edu](https://vahid2364.github.io)


## 🛠️ **Setup**
1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo/Phase-field_PyCalphad-Tutorial.git
   cd Phase-field_PyCalphad-Tutorial

## 📜 Acknowledgments

This tutorial was developed as part of a collaborative teaching effort at Texas A&M University. Special thanks to the Arroyave Research Group for their guidance and contributions to this work.

## 🎓 License

This project is licensed under the MIT License. See the LICENSE file for details.


### Key Improvements:
1. **Structured Overview**: Emphasized the significance of the project and its educational purpose.
2. **Expanded Instructions**: Detailed each step with clear terminal commands.
3. **Requirements Section**: Listed dependencies for clarity.
4. **Support and Acknowledgments**: Added a professional support section and acknowledgment of contributors.
5. **Engaging Tone**: Made the language more inviting and professional.

Feel free to modify links and additional content specific to your repository!
