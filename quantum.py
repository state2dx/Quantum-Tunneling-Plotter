import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Global variables to store wavefunctions
wave_original = None
wave_inside_barrier = None
wave_after_barrier = None
barrier_start = None
barrier_end = None


# Function to create Gaussian wavefunction (defining the format)
def gaussian_wave(x, a=1.0, x0=0.0, k=5.0):
    return np.exp(-a * (x - x0)**2) * np.cos(k * x)


# Function to calculate the energy for the particle
def calculate_energy(n, width):
    hbar = 1.0  # Reduced Planck's constant (simplified units)
    m = 1.0  # Particle mass (simplified units)
    L = width  # Width of the potential well
    energy = (n**2 * np.pi**2 * hbar**2) / (2 * m * L**2)
    return energy

#Function to update the plot
def update_plot():
    global wave_original, wave_inside_barrier, wave_after_barrier, barrier_start, barrier_end

    try:
        wave_expr = wave_var.get() #Fetches the wavefunction from the dialog box
        width = float(width_var.get()) #Fetches the width from the dialog box
        potential = float(potential_var.get()) #Fetches the potential from the dialog box
        dropdown_selection = energy_var.get() #Fetches the energy from the dropdown menu

        n = int(dropdown_selection.split()[0])
        energy = calculate_energy(n, width)
        display_var.set(f"Energy (E) = {energy:.3f} eV")

        x = np.linspace(-10, 10, 1000)
        dx = x[1] - x[0]  # Position spacing

        #Splitting into three regions
        # Define barrier boundaries
        barrier_start = 0
        barrier_end = width

        # Region 1: Before the barrier
        wave_original = eval(wave_expr, {"x": x, "np": np})

        # Region 2: Inside the barrier
        kappa = np.sqrt(2.0 * (potential - energy)) #k for region 2 (V0 potential)
        psi_at_start = wave_original[np.abs(x - barrier_start).argmin()]
        wave_inside_barrier = np.zeros_like(x)
        for i in range(len(x)):
            if barrier_start <= x[i] <= barrier_end:
                wave_inside_barrier[i] = psi_at_start * np.exp(-kappa * (x[i] - barrier_start))

        # Region 3: After the barrier
        wave_after_barrier = np.zeros_like(x)
        if energy > 0:
            k_free = np.sqrt(2.0 * energy) #k outside the barrier region 3 (0 potential)
            psi_at_end = wave_inside_barrier[np.abs(x - barrier_end).argmin()]
            a = 0.5
            x0 = barrier_end
            for i in range(len(x)):
                if x[i] > barrier_end:
                    wave_after_barrier[i] = np.exp(-a * (x[i] - x0) ** 2) * (
                        psi_at_end * np.cos(k_free * (x[i] - barrier_end))
                    )

        # Combine all regions into a single wavefunction for plotting
        wave_tunnel = np.copy(wave_original)
        wave_tunnel[(x >= barrier_start) & (x <= barrier_end)] = wave_inside_barrier[(x >= barrier_start) & (x <= barrier_end)]
        wave_tunnel[x > barrier_end] = wave_after_barrier[x > barrier_end]

        probability_density = np.abs(wave_tunnel)**2
        probability_density /= np.sum(probability_density) * dx

        # Plotting
        ax1.clear()
        ax2.clear()

        ax1.plot(x, wave_tunnel, label="Wavefunction (ψ)", color="blue", linewidth=1.5)
        ax1.axvline(barrier_start, color='red', linestyle='--', label="Barrier edges")
        ax1.axvline(barrier_end, color='red', linestyle='--')
        ax1.legend(fontsize=10)
        ax1.set_title("Wavefunction", fontsize=14, fontweight='bold')
        ax1.set_xlabel("Position (a.u.)", fontsize=12)
        ax1.set_ylabel("Amplitude", fontsize=12)
        ax1.grid(True, linestyle='--', alpha=0.6)

        ax2.plot(x, probability_density, label="Probability Density (|ψ|²)", color="green", linewidth=1.5)
        ax2.axvline(barrier_start, color='red', linestyle='--')
        ax2.axvline(barrier_end, color='red', linestyle='--')
        ax2.legend(fontsize=10)
        ax2.set_title("Probability Density", fontsize=14, fontweight='bold')
        ax2.set_xlabel("Position (a.u.)", fontsize=12)
        ax2.set_ylabel("Density", fontsize=12)
        ax2.grid(True, linestyle='--', alpha=0.6)

        canvas.draw()
        result_display.insert(tk.END, "Plot updated successfully!\n")
        result_display.yview(tk.END)

    except Exception as e:
        display_var.set(f"Error: {str(e)}")
        result_display.insert(tk.END, f"Error in update_plot: {str(e)}\n")


def calculate_expectation(calc_type):
    try:
        x_min, x_max = -10, 10
        x_points = 1000
        x = np.linspace(x_min, x_max, x_points)
        dx = x[1] - x[0]  

        # Get the wavefunction
        wave_expr = wave_var.get()
        wave = eval(wave_expr, {"x": x, "np": np})  # Evaluate wavefunction

        if np.all(wave == 0):
            result_display.insert(tk.END, "Wavefunction is zero!\n")
            return

        # Normalize the wavefunction
        probability_density = np.abs(wave)**2
        normalization_factor = np.sum(probability_density) * dx  # Integral of |ψ(x)|^2
        wave /= np.sqrt(normalization_factor)  # Normalize wavefunction

        # Recalculate the normalized probability density
        probability_density = np.abs(wave)**2

        if calc_type == "position":
            # Calculate Position Expectation: <x> = ∫ x * |ψ(x)|^2 dx
            expectation_position = np.sum(x * probability_density) * dx
            result_display.insert(tk.END, f"Position Expectation: {expectation_position:.3f} a.u.\n")
        
        elif calc_type == "momentum":
            # Fourier transform to momentum space (k-space)
            k = np.fft.fftfreq(x_points, d=dx) * 2 * np.pi  # Momentum space wavevector
            wave_k = np.fft.fft(wave)  # Fourier transform of wavefunction

            # Normalize in k-space
            prob_k = np.abs(wave_k)**2
            normalization_factor_k = np.sum(prob_k) * (k[1] - k[0])  # Integral in k-space
            prob_k /= normalization_factor_k  # Normalize probability density in momentum space

            # Expectation value of momentum <p> = ∫ k * |ψ(k)|^2 dk
            expectation_momentum = np.sum(k * prob_k) * (k[1] - k[0])  # Multiply by spacing in k-space
            result_display.insert(tk.END, f"Momentum Expectation: {expectation_momentum:.3f} ħk\n")
        
        # Update the display
        result_display.yview(tk.END)

    except Exception as e:
        result_display.insert(tk.END, f"Error in calculate_expectation: {str(e)}\n")
        result_display.yview(tk.END)

def plot_fourier_transform():
    global wave_original, wave_inside_barrier, wave_after_barrier, barrier_start, barrier_end

    try:
        x = np.linspace(-10, 10, 1000)
        dx = x[1] - x[0]  # Position spacing
        k = np.fft.fftfreq(len(x), d=dx) * 2 * np.pi  # Fourier space wavevector

        # Fourier transforms for each region
        fft_region1 = np.fft.fft(wave_original)
        fft_region2 = np.fft.fft(wave_inside_barrier)
        fft_region3 = np.fft.fft(wave_after_barrier)

        # Combine results
        wave_k = fft_region1 + fft_region2 + fft_region3
        prob_k = np.abs(wave_k) ** 2
        prob_k /= np.sum(prob_k) * (k[1] - k[0])  # Normalize in k-space

        # Plotting
        ax1.clear()
        ax2.clear()

        # Plot real and imaginary parts of ψ(k)
        ax1.plot(k, wave_k.real, label="Re(ψ(k))", color="blue", linewidth=1.5)
        ax1.plot(k, wave_k.imag, label="Im(ψ(k))", color="orange", linewidth=1.5)
        ax1.legend(fontsize=10)
        ax1.set_title("Fourier Transform of Wavefunction", fontsize=14, fontweight='bold')
        ax1.set_xlabel("Wavevector (k)", fontsize=12)
        ax1.set_ylabel("Amplitude", fontsize=12)
        ax1.grid(True, linestyle='--', alpha=0.6)

        # Plot probability density in k-space
        ax2.plot(k, prob_k, label="Probability Density (|ψ(k)|²)", color="green", linewidth=1.5)
        ax2.legend(fontsize=10)
        ax2.set_title("Probability Density in k-space", fontsize=14, fontweight='bold')
        ax2.set_xlabel("Wavevector (k)", fontsize=12)
        ax2.set_ylabel("Density", fontsize=12)
        ax2.grid(True, linestyle='--', alpha=0.6)

        canvas.draw()
        result_display.insert(tk.END, "Fourier transform plotted successfully!\n")
        result_display.yview(tk.END)

    except Exception as e:
        result_display.insert(tk.END, f"Error in plot_fourier_transform: {str(e)}\n")
        result_display.yview(tk.END)


def fourier_trans():
    global wave_original, wave_inside_barrier, wave_after_barrier, barrier_start, barrier_end

    try:
        x = np.linspace(-10, 10, 1000)
        dx = x[1] - x[0]  # Position spacing
        k = np.fft.fftfreq(len(x), d=dx) * 2 * np.pi  # Fourier space wavevector

        # Fourier transforms for each region
        fft_region1 = np.fft.fft(wave_original)
        fft_region2 = np.fft.fft(wave_inside_barrier)
        fft_region3 = np.fft.fft(wave_after_barrier)

        # Combine results
        wave_k = fft_region1 + fft_region2 + fft_region3

        # Display the Fourier transform results in the text box
        result_display.insert(tk.END, f"Fourier Transform (k-space values):\n{k}\n")
        result_display.insert(tk.END, f"Real part of Fourier Transform:\n{wave_k.real}\n")
        result_display.insert(tk.END, f"Imaginary part of Fourier Transform:\n{wave_k.imag}\n")
        result_display.insert(tk.END, f"Probability Density (|ψ(k)|²):\n{np.abs(wave_k)**2}\n")
        result_display.yview(tk.END)

    except Exception as e:
        result_display.insert(tk.END, f"Error in fourier_trans: {str(e)}\n")
        result_display.yview(tk.END)

root = tk.Tk()
root.title("Quantum Tunneling Simulation")
root.geometry("1200x800")
root.configure(bg="#2C3E50")

#Defining the plotting area
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
fig.tight_layout(pad=4.0, h_pad=3.0, rect=[0, 0, 1, 0.96])
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

#Defining defaults
wave_var = tk.StringVar(value="np.exp(-0.5 * x**2) * np.cos(5 * x)")
width_var = tk.StringVar(value="2.0")
potential_var = tk.StringVar(value="5.0")
display_var = tk.StringVar()
energy_dropdown_values = [f"{n} ({calculate_energy(n, 2.0):.3f} eV)" for n in range(1, 11)]
energy_var = tk.StringVar(value=energy_dropdown_values[0])

#Defining Control Panel area
control_panel = tk.Frame(root, bg="#34495E", padx=20, pady=20)
control_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

#Dialog box, Button and Dropdown box creationo
tk.Label(control_panel, text="Wavefunction:", font=("Arial", 14, "bold"), bg="#34495E", fg="white").grid(row=0, column=0, sticky="e")
tk.Entry(control_panel, textvariable=wave_var, width=30, font=("Arial", 14), fg="black", bg="#ECF0F1").grid(row=0, column=1)
tk.Label(control_panel, text="Energy Level (n):", font=("Arial", 14, "bold"), bg="#34495E", fg="white").grid(row=1, column=0, sticky="e")
energy_dropdown = ttk.Combobox(
    control_panel, textvariable=energy_var, 
    values=energy_dropdown_values, width=30, font=("Arial", 14), state="readonly"
)
energy_dropdown.grid(row=1, column=1)

tk.Label(control_panel, text="Barrier Width (a.u.):", font=("Arial", 14, "bold"), bg="#34495E", fg="white").grid(row=2, column=0, sticky="e")
tk.Entry(control_panel, textvariable=width_var, width=10, font=("Arial", 14), fg="black", bg="#ECF0F1").grid(row=2, column=1)

tk.Label(control_panel, text="Potential (eV):", font=("Arial", 14, "bold"), bg="#34495E", fg="white").grid(row=3, column=0, sticky="e")
tk.Entry(control_panel, textvariable=potential_var, width=10, font=("Arial", 14), fg="black", bg="#ECF0F1").grid(row=3, column=1)

tk.Button(control_panel, text="Update Plot", command=update_plot, font=("Arial", 14, "bold"), bg="#1ABC9C", fg="white").grid(row=4, column=0, columnspan=2, pady=10)

tk.Button(control_panel, text="Calc Position", command=lambda: calculate_expectation("position"), font=("Arial", 14), bg="#5DADE2", fg="black").grid(row=5, column=0, columnspan=2, pady=5)
tk.Button(control_panel, text="Calc Momentum", command=lambda: calculate_expectation("momentum"), font=("Arial", 14), bg="#5DADE2", fg="black").grid(row=6, column=0, columnspan=2, pady=5)

result_display_frame = tk.Frame(root, bg="#2C3E50")
result_display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

result_display = tk.Text(result_display_frame, height=20, width=50, font=("Arial", 12), bg="#ECF0F1", fg="black", wrap=tk.WORD)
result_display.pack(fill=tk.BOTH, expand=True)
result_display.insert(tk.END, "Welcome to the Quantum Tunneling Simulation!\n")

tk.Button(control_panel, text="Plot Fourier Transform", command=plot_fourier_transform, font=("Arial", 14), bg="#E67E22", fg="black").grid(row=9, column=0, columnspan=2, pady=10)
tk.Button(control_panel, text="Calculate Fourier Transform", command=fourier_trans, font=("Arial", 14), bg="#E67E22", fg= "black").grid(row=10, column=0, columnspan=2, pady=10)

root.mainloop()
